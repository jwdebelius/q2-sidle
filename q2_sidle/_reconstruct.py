import itertools as it
import os

import biom
import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd

from qiime2 import Metadata, Artifact
from qiime2.plugin import ValidationError
from q2_sidle._utils import (_setup_dask_client, 
                             _convert_generator_to_seq_block,
                             _convert_generator_to_delayed_seq_block, 
                             _convert_seq_block_to_dna_fasta_format,
                             _check_manifest,
                             degen_reps,
                             _read_manifest_files
                             )

def reconstruct_counts(manifest: Metadata,
    count_degenerates: bool=True,
    per_nucleotide_error: float=0.005,
    max_mismatch: int=2,
    min_abund: float=1e-5,
    region_normalize: bool=True,
    block_size : int=10000,
    debug: bool=False, 
    n_workers: int=0,
    client_address: str=None,
    ) -> (biom.Table, Metadata, pd.DataFrame):
    """
    Reconstructs regional alignments into a full length 16s sequence

    Parameters
    ----------
    manifest: Metadata
         A manifest file describing the relationship between regions and their
        alignment mapping. The manifest must have at least three columns
        (`kmer-map`, `alignment-map` and `frequency-table`) each of which
        contains a unique filepath. The regions should be sorted and the 
        region labels must match between the kmer map and the alignment map.
    count_degenerates: bool
        Whether degenerate sequences should be counted as unique kmers during
        reconstruction or only unique sequences should be counted.
    per_nucleotide_error: float
        The assumed per-nucleotide error rate from amplifican and sequencing.
        The value used here is based on benchmarking from the original smurf
        paper
    max_mismatch: int
        The maximum number of nucleotides that can differ between a kmer and
        a sequence for it to be considered a match
    min_abund: float
        The minimum relative abundance for a feature to be retained during 
        reconstruction. The default from the original smurf algorithm is '
          '1e-10, the number here may depend on the total number of '
          'sequences in your sample.
    debug: bool
        Whether the function should be run in debug mode (without a client)
        or not. `debug` superceeds all options
    n_workers: int, optional
        The number of jobs to initiate. When `n_workers` is 0, the cluster 
        will be able to access all avalaibel resources.
    client_address: str
        The IP address for an existing dask client/cluster

    Returns
    -------
    biom.Table
        The reconstructed, region normalized counts
    Metadata
        A summary of the statitics for the regional map describing the number
        of regions mapped to each reference sequence and the number of kmers.
        The kmer mapping estimate can account for degeneracy when the 
        `count_degenerates` flag is used or can ignore degenrate sequences in
        mapping
    DataFrame
        A map between the final kmer name and the original database sequence.
        Useful for reconstructing taxonomy and trees.
    """
    
    # Sets up the client
    _setup_dask_client(debug=debug, cluster_config=None,  
                       n_workers=n_workers, address=client_address)

    # Checks the manifest
    _check_manifest(manifest)
    
    # Gets the order of the regions
    region_order = \
        manifest.get_column('region-order').to_series().astype(int)
    region_order = region_order - region_order.min()
    region_names = {i: r for r, i in region_order.items()}
    num_regions = len(region_order)

    # Imports the alignment maps and gets the kmers that were aligned
    align_map = pd.concat(
        axis=0, 
        sort=False, 
        objs=_read_manifest_files(manifest, 'alignment-map', 
                                  'FeatureData[KmerAlignment]', pd.DataFrame)
        )
    align_map.reset_index(drop=True, inplace=True)
    align_map[['mismatch', 'length']] = \
        align_map[['mismatch', 'length']].astype(int)
    align_map.replace({'region': region_order}, inplace=True)
    kmers = _get_unique_kmers(align_map['kmer'])

    print('Regional Alignments Loaded')

    # ### Untangles the database to get the unique regional mapping

    # Filters database down to kmers which are present in the sequences 
    # because otherwise we're trying to untangle a huge amount of data and it 
    # just gets memory intensive and slow
    kmer_map = dd.concat(
        axis=0,
        dfs=_read_manifest_files(manifest, 'kmer-map', 
                                 'FeatureData[KmerMap]', dd.DataFrame)
    )
    kmer_map = kmer_map.reset_index().set_index('db-seq').reset_index()
    kmer_map['region'] = kmer_map['region'].replace(region_order)
    kmer_map  = kmer_map.loc[kmer_map['db-seq'].isin(kmers)]
    
    print('Regional Kmers Loaded')

    # Builds database mapping bettween the kmer, the original database
    # sequence and the original name
    db_map = _untangle_database_ids(kmer_map.copy(),
                                    num_regions=num_regions,
                                    block_size=block_size / 10,
                                    )
    print('Database map assembled')    

    ### Summarizes the database
    kmer_map = kmer_map.compute()
    kmer_map.set_index('db-seq', inplace=True)
    kmer_map['clean_name'] = db_map
    kmer_map.reset_index(inplace=True)

    db_summary = _count_mapping(kmer_map.reset_index(), 
                                count_degenerates, 
                                kmer='seq-name')
    if region_normalize == False:
        print('normalize!')
        db_summary['num-regions'] = 1

    ### Constructs the regional alignment
    align_mat = _construct_align_mat(align_map,
                                    sequence_map=db_map.to_dict(),
                                    seq_summary=db_summary,
                                    nucleotide_error=per_nucleotide_error, 
                                    max_mismatch=max_mismatch, 
                                    blocksize=block_size,
                                    )
    print('Alignment map constructed')
    
    ### Solves the relative abundance
    counts = pd.concat(
        axis=1, 
        sort=False, 
        objs=_read_manifest_files(manifest, 'frequency-table', 
                                 'FeatureTable[Frequency]', pd.DataFrame)).T
    counts.fillna(0, inplace=True)
    print('counts loaded')

    # We have to account for the fact that some of hte ASVs may have been 
    # discarded because they didn't meet the match parameters we've set or 
    # because they're not in the database.
    keep_asvs = list(set(align_mat['asv'].values) & (set(counts.index)))
    unaligned_counts = counts.copy().drop(keep_asvs).sum(axis=1)
    counts = counts.loc[keep_asvs]
    align_mat = align_mat.loc[align_mat['asv'].isin(keep_asvs)]
    keep_kmers = align_mat['clean_name'].unique()
    db_summary = db_summary.loc[keep_kmers]

    # Normalizes the alignment table
    n_table = counts / counts.sum(axis=0)
    print('counts normalized')

    # Performs the maximum liklihood reconstruction on a per-sample basis. 
    # Im not sure if this could be refined to optimize the alogirthm
    # to allow multiple samples ot be solved together, but... eh?
    sample_ids = n_table.columns.values
    n_samples = len(sample_ids)

    rel_abund = _solve_iterative_noisy(align_mat=align_mat, 
                                       table=n_table,
                                       min_abund=min_abund,
                                       seq_summary=db_summary)

    db_summary = db_summary.loc[rel_abund.index]
    print('Relative abundance calculated')

    # Puts together the regional normalized counts
    count_table = _scale_relative_abundance(align_mat=align_mat,
                                            relative=rel_abund,
                                            counts=counts,
                                            seq_summary=db_summary).fillna(0)
    count_table = biom.Table(count_table.values,
                             sample_ids=count_table.columns,
                             observation_ids=count_table.index)
    count_table.filter(lambda v, id_, md: v.sum() > 0, axis='observation')
    summary = db_summary.loc[count_table.ids(axis='observation')]
    summary['mapped-asvs'] = \
        align_mat.groupby('clean_name')['asv'].apply(lambda x: '|'.join(x))
    summary.index.set_names('feature-id', inplace=True)
    summary = Metadata(summary)

    # Puts together a kmer-based mapping which describes the regions covered
    kmer_map.sort_values(['db-seq', 'region'], ascending=True, inplace=True)
    kmer_map.drop_duplicates(['db-seq', 'region'], inplace=True)
    kmer_map['region'] = kmer_map['region'].astype(int)

    first_fwd = kmer_map.drop_duplicates('db-seq', keep='first')
    first_fwd = first_fwd.set_index('db-seq')[['fwd-primer']]
    last_rev = kmer_map.drop_duplicates('db-seq', keep='last')
    last_rev = last_rev.set_index('db-seq')[['fwd-primer', 'kmer-length']]
    
    mapping = pd.concat(axis=1, sort=False, objs=[
        db_map[db_map.isin(count_table.ids(axis='observation'))],
        first_fwd.add_prefix('first-'),
        last_rev.add_prefix('last-'),
        ])
    mapping.index.set_names('db-seq', inplace=True)
    mapping.dropna(subset=['clean_name'], inplace=True)
    mapping.sort_index(inplace=True)

    return count_table, summary, mapping


def _construct_align_mat(match, sequence_map, seq_summary, 
    nucleotide_error=0.005, max_mismatch=2, kmer_name='kmer', asv_name='asv',
    seq_name='clean_name', miss_col='mismatch', blocksize=5000):
    """
    Constructs a long form matrix describing alignmnet

    Parameters
    ----------
    matches: DataFrame
        A list of long-form dataframes describing the kmer-region matches for
        each region. The output of `align_kmers`
    sequence_map : dict
        Maps the regional database kmer names back to the original database
        identifiers
    nucleotide_error: float, optional
        The assumed per-nucleotide error rate.
    max_mismatch : int, optional
        The maximum number of nucelotides which can differ between an ASV
        and the asigned kmer
    kmer_name: str, optional
        The column name in `mismatch` that identifies the reference sequences
        the kmers came from
    seq_name: str, optional
        The label for the original sequence in the database a kmer comes from
        (assuming that we're collapsing degenerate sequences back to orginal)
        fits. 
    asv_name : str, optional
        The column in `mismatch` which identifies the ASV identifer for each
        sequence that was mapped to a kmer.
    miss_col : str, optional
        The column in `mismatch` that describes the mismatch between the kmers

    Returns
    -------
    DataFrame
        A mapping of the match probability and error rate between matched
        kmers for iterative processing and reconstruction.

    """

    # Expands the kmer identifiers and then maps back to the database
    align_mat = dd.from_pandas(
        _expand_duplicate_sequences(match, id_col=kmer_name),
        chunksize=blocksize
    )

    align_mat['db-seq'] = \
        align_mat[kmer_name].apply(lambda x: x.split('@')[0],
                                   meta=('kmer_name', 'str'))
    align_mat[seq_name] = align_mat['db-seq'].replace(sequence_map)

    # Collapses into ASV-ref seq combos where for degenerates, the closet
    # match is retained as the sequence.
    asv_ref = align_mat.groupby([seq_name, asv_name, 'region'])
    asv_ref = asv_ref[[miss_col, 'length']].min().reset_index()

    ## E(ih) : Pr(read = i | kmer = h)
    # Calculates the probability of an error as the probability that bases
    # are correct (1 - error rate)^(right bases)  and the probability the 
    # mismatches are wrong (error_rate / 3)^mismatch (where we divide by 3)
    # for any of the other bases to subsitute.
    asv_ref['match_prob'] = \
        (np.power(1 - nucleotide_error, asv_ref['length'] - asv_ref[miss_col])
         * (np.power(nucleotide_error / 3, asv_ref[miss_col]))).astype(float)
    asv_ref['perfect_prob'] = \
        np.power(1-nucleotide_error, asv_ref['length'])
    asv_ref['max_error_prob'] = (
        np.power(1 - nucleotide_error, asv_ref['length'] - max_mismatch) * 
        np.power(nucleotide_error / 3, max_mismatch)
    )

    # Calculates the maximum error rate to get the mismatch. 0.1 multiplier comes
    # directly from the SMURF code
    # https://github.com/NoamShental/SMURF/blob/master/mFiles/solve_iterative_noisy.m#L53
    asv_ref['error_thresh'] = \
        (asv_ref['perfect_prob'] - 0.1*asv_ref['max_error_prob'] - 
            np.spacing(1))
    asv_ref['perf_match'] = \
        (asv_ref['match_prob'] > asv_ref['error_thresh'])

    # # ### Computes Mij: Pr(kmer = j | read = j)
    # Normalizes the alignment probability based on the pre-region 
    # amplification
    region_amp = seq_summary['num-regions']

    asv_ref = asv_ref.set_index(seq_name)
    asv_ref['region-norm'] = region_amp.astype(int)

    asv_ref['norm'] = (asv_ref['match_prob'] / asv_ref['region-norm'])

    asv_ref  = asv_ref.reset_index().compute()
    asv_ref.drop(columns='region-norm', inplace=True)

    return asv_ref


def _count_mapping(long_, count_degen, kmer='kmer', region='region', 
    clean_name='clean_name', db_seq='db-seq'):
    """
    Counts the number of sequences which have been mapped 
    """
    # If we're counting degenerates, we'll keep the kmers whcih contain
    # degeneracies. Otherwise, we work with the db_seqs which don't.
    if count_degen:
        long_.drop_duplicates([kmer, region, clean_name], 
                              inplace=True)
    else:
        long_.drop_duplicates([db_seq, region, clean_name], 
                              inplace=True)

    # Gets the regional counts
    regional_seqs = \
        long_.groupby([clean_name, region]).count()[[kmer]]
    regional_seqs.reset_index(inplace=True)

    # Then, we count the data.
    counts = pd.DataFrame(
        data=[regional_seqs.groupby('clean_name')[region].count(),
              regional_seqs.groupby('clean_name')[kmer].sum(),
              regional_seqs.groupby('clean_name')[kmer].mean(),
              regional_seqs.groupby('clean_name')[kmer].std().fillna(0),
              # regional_seqs.groupby('clean_name')['']
              ],
        index=['num-regions', 
               'total-kmers-mapped', 
               'mean-kmer-per-region',
               'stdv-kmer-per-region'
               ],
        ).T
    return counts


def _detangle_names(long_):
    """
    Splits up the mapping into clearer clusters

    Parameters
    ----------
    long_: DataFrame
        A mapping between the database sequence (`db-seq`), the matched 
        sequence  (`clean_name`), and a column (`counter`) containing an 
        integer describing how many sequences are mapped to the 
        database sequence in a sorted order.
        It is important that this be a long, defense association: condensed/
        upper triangle matrices will not work here.

    Returns
    -------
    Series
        A mapping between the databset sequence and a unique group name.
    """

    long_.sort_values(['db-seq', 'clean_name'], ascending=True, inplace=True)
    long_['counter'] = 1
    long_['counter'] = long_.groupby(['db-seq'])['counter'].cumsum() - 1
    
    # gets the horizontal mapping that pivots the dataframe into a square
    # distance/mapping matrix
    map_ = \
        long_.set_index(['db-seq', 'clean_name'])['counter'].unstack().notna()
    diagnonal = (np.tril(np.ones(map_.shape)) == np.triu(np.ones(map_.shape)))
    map_ = (map_ | diagnonal)
    
    # Clusters the sequences
    map_.sort_values(list(map_.columns), inplace=True)
    map_ = map_[map_.index]
        
    # Calculates the number of connections within the cluster
    num_connections = (long_.groupby('db-seq')['counter'].max() + 1)
    long_.set_index('clean_name', inplace=True)
    long_['co-connection'] = num_connections[long_.index]
    co_min = long_.groupby('db-seq')['co-connection'].min()
    connection_penalty = (num_connections - co_min + 1)

    # Builds a mask that penalizes for the number of connections that should
    # be present
    penalty_v = (np.ones((len(connection_penalty), 1)) * 
                 connection_penalty.loc[map_.index].values)
    penalty_mask = (penalty_v == penalty_v.T)
    
    overlap = penalty_mask & map_
    overlap = overlap.mask(~overlap, np.nan)
    overlap = overlap.unstack().dropna().reset_index()

    # # In some cases, there are sequences that have the same number of 
    # # connections but don't  connect to the same sequences. So,for pairs 
    # # of unique sequences which are initally mapped, we want to check the
    # # number of sequences that overlap
    shared = \
        overlap.groupby('db-seq')['clean_name'].apply(lambda x: set(x.values))
    def check_shared(x):
        seqs0 = shared.loc[x['clean_name']]
        seqs1 = shared.loc[x['db-seq']]
        if len(seqs0) != len(seqs1):
            return 0
        else:
            return len(seqs0 & seqs1) - 1

    check_diff = overlap['clean_name'] != overlap['db-seq']

    if check_diff.any():
        overlap['seq_count'] = overlap.loc[check_diff].apply(check_shared, axis=1)
    else:
        overlap['seq_count'] = np.nan
    overlap['matches'] = (overlap['seq_count'] == overlap[0]) | ~check_diff

    overlap.sort_values(['db-seq', 'clean_name'], 
                        inplace=True, 
                        ascending=True)
    overlap.drop_duplicates(['db-seq', 'clean_name'], 
                            inplace=True)
    new_name = \
        overlap.groupby('db-seq')['clean_name'].apply(lambda x: "|".join(x))

    return new_name


def _expand_duplicate_sequences(df, id_col, delim='|'):
    """
    Expands delimited IDs into rows with unique identifiers

    Paramters
    ---------
    df : Dataframe
        The dataframe with the duplicated ids
    id_col : str
        The column in `df` containing the ID information
    data_cols : list, optional
        The columns from the dataframe ot be retained
    delim : str
        The delimiting string

    Returns
    -------
    DataFrame
        An expanded dataframe.
    """
    data_cols = df.drop(columns=[id_col]).columns
    wide_exp = pd.concat(axis=1, objs=[
        df[id_col].apply(lambda x: pd.Series(x.split(delim))),
        df[data_cols]
        ])
    # Converts to long form and drops the duplicated anything without an ID
    long_exp = wide_exp.melt(
        id_vars=data_cols,
        value_vars=wide_exp.drop(columns=data_cols).columns,
        value_name=id_col
    ).drop(columns=['variable']).dropna(subset=[id_col])
    long_exp = long_exp[np.hstack([[id_col], data_cols.values])]
    long_exp.sort_values(by=list(long_exp.columns), inplace=True)
    long_exp.reset_index(inplace=True, drop=True)

    return long_exp


def _get_clean(df):
    clean_kmers = \
        df.groupby(['db-seq', 'region'])['value'].apply(
            lambda x: "|".join(x.values))
    return clean_kmers.reset_index()


def _get_shared_seqs(df):
    wide = df.set_index(['db-seq', 'region'])['kmer'].apply(
            lambda x: pd.Series(x.split("|"))).reset_index()
    long_ = wide.melt(id_vars=['db-seq', 'region']).dropna()
    long_['value'] = long_['value'].apply(lambda x: x.split("@")[0])
    long_.drop_duplicates(['db-seq', 'region', 'value'], inplace=True)
    
    return long_.drop(columns=['variable'])


def _get_unique_kmers(series):
    kmers = np.hstack([[a.split("@")[0] for a in kmer.split('|')] 
                       for kmer in series])
    return np.sort(np.unique(kmers))


def _scale_relative_abundance(align_mat, relative, counts, seq_summary,
    seq_name='clean_name', asv_name='asv'):
    """
    Scales the relative abundance data to give sequence count data
    Parameters
    ----------
    align_mat: DataFrame
        A mapping of the match probability and error rate between matched
        kmers for iterative processing and reconstruction.
    relative : DataFrame
        a sequence x sample table giving the relative frequency of each
        sequence in each sample.
    counts : DataFrame
        The relationship between each ASV and the number of counts observed
        for that sequence
    sequence_summary : DataFrame
        A summary of the kmers which were mapped to each of the refernece
        sequences. This describes the number regions covered, the total
        number of kmers mapped, the average and standard deviation in the 
        number per region.
    seq_name: str, optional
        The label for the original sequence in the database a kmer comes from
        (assuming that we're collapsing degenerate sequences back to orginal)
        fits. 
    asv_name : str, optional
        The column in `mismatch` which identifies the ASV identifer for each
        sequence that was mapped to a kmer.
    Returns
    -------
    DataFrame
        A table describing the region-scaled counts mapped to each reference 
        sequence where the total counts represents the sum of counts across 
        all sequencing regions. 
    """
    align_mat = align_mat.loc[align_mat[seq_name].isin(relative.index)]
    align = align_mat.pivot_table(values='norm', index=asv_name, 
                                  columns=seq_name, fill_value=0
                                  )

    # And then we get indexing because I like to have the indexing
    # references
    align_seqs = align.columns
    align_asvs = align.index
    align = align.values

    relative = relative.loc[align_seqs]
    counts = counts.loc[align_asvs]
    samples = relative.columns
    spacer =  np.spacing(1)

    # A simple per-sample function to solve the counts because frequency is fun
    # and pandas rounding may be what's killing this
    @dask.delayed
    def _solve_sample(align, freq, count):
        # Probability a count shows up in a particular ASV/reference pairing
        p_r_and_j = freq.values.T * align
        # Normalized for the number of asvs mapped to that reference
        p_r_given_j = p_r_and_j / np.atleast_2d(p_r_and_j.sum(axis=1) + spacer).T
        # sums the counts
        count_j_given_r = count.values * p_r_given_j

        return np.atleast_2d(count_j_given_r.sum(axis=0))

    scaled_counts = pd.DataFrame(
        data=np.vstack(dask.compute(*[
            _solve_sample(align, relative[[sample]], counts[[sample]])
            for sample in samples])),
        index=samples,
        columns=align_seqs
        )
    counts_r = scaled_counts / seq_summary['num-regions']


    return counts_r.round(0).T


def _solve_iterative_noisy(align_mat, table, seq_summary, tolerance=1e-7,
    min_abund=1e-10, num_iter=1e5, 
    seq_name='clean_name', asv_name='asv', threads=1):
    """
    Maps ASV abundance to reference sequences
    Parameters
    ----------
    align_mat: DataFrame
        A mapping of the match probability and error rate between matched
        kmers for iterative processing and reconstruction.
    table: DataFrame
        A table combining all ASV sequences per sample (coutns for all regions) 
        where samples are columns and ASVs are rows (biom transposition) which 
        has been normalized per sample.
    sequence_summary : DataFrame
        A summary of the kmers which were mapped to each of the refernece
        sequences. This describes the number regions covered, the total
        number of kmers mapped, the average and standard deviation in the 
        number per region.
    tolerance: float, optional
        The error tolerance in alignmnet
    min_abund : float, optional
        The minimum relative abundance to be considered "present". Bacteria
        below this threshhold will be discarded as having a 0-value
    num_iter: int, optional
        The number of iterations for fitting
    seq_name: str, optional
        The label for the original sequence in the database a kmer comes from
        (assuming that we're collapsing degenerate sequences back to orginal)
        fits. 
    asv_name : str, optional
        The column in `mismatch` which identifies the ASV identifer for each
        sequence that was mapped to a kmer.
    Returns 
    -------
    DataFrame
        a sequence x sample table giving the relative frequency of each
        sequence in each sample.
    Raises
    ------
    ValueError
        If there are not the same number of matches and tables
    ValueError
        If there are not the same number of read lengths and matches if 
        a list of read lengths was specified.
    """

    # Gets the sparse alignment matrix
    # Gets the alignment matrix because its cheaper outside the loop
    # We get a sparse matrix because hopefully it works better.
    align_mat = align_mat.loc[align_mat['asv'].isin(table.index)]
    align = align_mat.pivot_table(values='norm', index=asv_name, 
                                  columns=seq_name, fill_value=0
                                  )
    align = align.loc[table.index]
    # And then we get indexing because I like to have the indexing
    # references
    align_seqs = align.columns
    align_asvs = align.index

    recon = []
    for sample, col_ in table.iteritems():
        freq_ = dask.delayed(_solve_ml_em_iterative_1_sample)(
            align=align.loc[col_ > 0].values,
            abund=dask.delayed(col_[col_ > 0].values),
            sample=sample,
            align_kmers=align_seqs,
            num_iter=num_iter,
            tolerance=tolerance,
            min_abund=min_abund,
            )
        df = dask.delayed(pd.DataFrame)(freq_, columns=[sample])
        recon.append(df)
    # Puts together the table as optimiziation
    recon = pd.concat(axis=1, sort=False, objs=dask.compute(*recon))
    recon.dropna(how='all', axis=0, inplace=True)
    recon.fillna(0, inplace=True)

    # # Region normalization
    recon = recon / seq_summary.loc[recon.index, ['num-regions']].values
    recon = recon / recon.sum(axis=0)
    recon.sort_index(axis='index', inplace=True)
    recon.sort_index(axis='columns', inplace=True)

    return recon


def _solve_ml_em_iterative_1_sample(align, abund, align_kmers, sample, 
    num_iter=10000,
    tolerance=1e-7,  min_abund=1e-10, 
    kmer_name='kmer', asv_name='asv'):
    """
    Iterative expected maximization of coutns and alignment
    
    Python implementation of `ml_em_iterative.m` from SMURF.
    See http://bjlkeng.github.io/posts/the-expectation-maximization-algorithm/
    as at least an original estimation of what's happening  
    Parameters
    ----------
    align : 
    n : int
        The number of regions explored
    n_iter: int
        the number of iterations to try before giving up
    tol : float
        The error tolerance???
    Returns
    --------
    DataFrame
        The assigned relative abundance of each read
    # TO DO: can this be solved for multipl erads at once?
    """

    # Our starting bact_freq estimate is here.
    bact_freq = np.dot(align.T, abund) / np.sum(np.dot(align.T, abund))

    # Solves the mixture. I'm not sure if this can be reimplemented as a
    # gaussian mixed models problem if we have multiple samples? Because like
    # that might actuallky make it easier? Maybe???
    for i in np.arange(0, num_iter):
        ### Expectation
        # assign theta estimate for each bacteria
        theta_i = np.dot(align, bact_freq)
        
        ### Maximization
        # Adjusts the abundance and bacterial frequency
        r_weighted = abund / (theta_i + np.spacing(1))
        bact_factor = np.dot(r_weighted.T, align)
        
        # Computes errro 
        error = np.dot(np.absolute(1 - bact_factor), bact_freq)
        
        # And then we get the new frequency... 
        bact_freq = bact_freq * bact_factor
        
        if error < 1e-7:
            break

        high_enough = bact_freq > min_abund
        # if sum(~high_enough) > 0:
            # print(np.arange(0, len(bact_freq))[~high_enough])
        align = align[:, high_enough]
        bact_freq = bact_freq[high_enough]
        align_kmers = align_kmers[high_enough]

    # And then we do hard threshholding
    bact_freq[bact_freq <= min_abund] = 0
    bact_freq = bact_freq / bact_freq.sum()

    return pd.Series(bact_freq, align_kmers, name=sample)


def _sort_untidy(df, clean_seqs):
    if df['tidy']:
        return df['shared-set']
    else:
        return df['shared-set'] - clean_seqs
    

def _tidy_sequence_set(clean_kmers, clean_seqs):
    """
    Iteratively cleans kmers by correcting overlapping sets
    of sequences
    """
    
    # Gets the sequences that are still shared with other sequences
    shared_seqs = clean_kmers.loc[~clean_kmers['tidy']].copy()
    shared_seqs = shared_seqs.groupby('db-seq')['shared-set'].apply(
        lambda x: set.intersection(*x), 
        meta=('shared-set', 'object')
    )
    shared_seqs = shared_seqs.reset_index()
    shared_seqs.columns = ['db-seq', 'shared-set']
    shared_seqs['kmer'] = shared_seqs['shared-set'].apply(
        lambda x: '|'.join(sorted(x)), 
        meta=('shared-set', 'str')
    )
    # Updates the mapped sequence
    clean_check = shared_seqs.groupby('kmer').apply(
        lambda x: (set.union(*x['shared-set'].values)) == set(x['db-seq']),
        meta=(None, 'bool')).compute()
    if clean_check.any():
        cleaned = np.hstack([
            x.split("|") for x in clean_check.index[clean_check].values
        ])
        clean_seqs = clean_seqs | set(cleaned)
    else:
        clean_seqs = clean_seqs

    clean_kmers['tidy'] = clean_kmers['db-seq'].isin(clean_seqs)
    clean_kmers['shared-set'] = \
        clean_kmers.apply(_sort_untidy, 
                          clean_seqs=clean_seqs,
                          axis=1, 
                          meta=(None, 'object')
                         )
    return clean_kmers, clean_seqs


def _untangle_database_ids(region_db, num_regions, block_size=500):
    """
    Untangles regional sequence identifier to map sequences to each other

    Parameters
    ----------
    region_db: DataFrame
        A sequence array as a dataframe containing the trimmed region
    kmers: list
        A list of kmers to keep for the analysis
    block_size : int, optional
        The number of sequence groups whcih should be combined for mapping

    Returns
    -------
    DataFrame
        a map between the regional database identifiers and the sequences
        based on the combinations of regions used.
    DataFrame
        A summary of the kmers which were mapped to each of the refernece
        sequences. This describes the number regions covered, the total
        number of kmers mapped, the average and standard deviation in the 
        number per region.
    """
    # Cleans up the kmers by pulling off shared labels
    shared = [dask.delayed(_get_shared_seqs)(df) 
              for df in region_db.to_delayed()]
    clean_kmers = dd.from_delayed([dask.delayed(_get_clean)(df) 
                                   for df in shared])
    clean_kmers['tidy'] = False
    # Formats for an id set
    clean_kmers['shared-set'] = clean_kmers['value'].apply(
        lambda x: set(x.split("|")), meta=('value', 'object'))

    # Sets up a while loop. Because somehow it makes sense here?
    untidy = True
    clean_seqs = set([])
    last_tidy = len(clean_kmers)

    while untidy:
        clean_kmers, clean_seqs = _tidy_sequence_set(clean_kmers, clean_seqs)
        tidy_seqs = ~clean_kmers['tidy'].compute()
        untidy = tidy_seqs.any()
        if (last_tidy == tidy_seqs.sum()):
            break
        last_tidy = tidy_seqs.sum()

    # Maps the tidy kmers into a single kmer
    db_map1 = clean_kmers.loc[clean_kmers['tidy']
                             ].groupby('db-seq')['shared-set'].apply(
        lambda x: '|'.join(sorted(set.intersection(*x))), 
        meta=('shared-set', str)
    )
    db_map1.name = 'clean_name'
    # If anything else remains untidy, we need to continue mapping...
    if untidy:
        # Finds the sequences that need to be mapped and pulls out the data
        bad_ = ~clean_kmers['tidy']
        to_map = clean_kmers.loc[bad_].compute()
        unique_ids = to_map['db-seq'].unique()

        to_map = to_map.groupby('db-seq')['shared-set'].apply(
                lambda x: pd.Series(sorted(set.intersection(*x.values))),
                ).reset_index()
        to_map.columns = ['db-seq', 'counter', 'clean_name']

        db_map2 = _detangle_names(to_map) 
    else:
        db_map2 = pd.Series([], 
                            index=pd.Index([], name='db-seq'),  
                            name='clean_name')
    
    db_map = pd.concat(axis=0, objs=[db_map1.compute(), db_map2])
    
    return db_map

