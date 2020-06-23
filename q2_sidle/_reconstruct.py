import itertools as it
import os

import biom
import dask
import numpy as np
import pandas as pd
import sparse as sp

from qiime2 import Metadata, Artifact
from qiime2.plugin import ValidationError
from q2_sidle._utils import (_setup_dask_client, 
                             _convert_generator_to_seq_block,
                             _convert_generator_to_delayed_seq_block, 
                             _convert_seq_block_to_dna_fasta_format,
                             degen_reps,
                             )

def reconstruct_counts(manifest: Metadata,
    count_degenerates: bool=True,
    per_nucleotide_error: float=0.005,
    max_mismatch: int=2,
    min_abund: float=1e-10,
    debug: bool=False, 
    n_workers: int=0
    ) -> (biom.Table, Metadata, pd.DataFrame):
    """
    Reconstructs regional alignments into a full length 16s sequence

    Parameters
    ----------
    kmer_map : DataFrame
        Okay, technically if this works, this is a list of the per-region 
        dataframes that describe the relationship between kmer names and the 
        original database sequence names
    regional_alignment : DataFrame
        Again, a list of dataframes that map the regional kmers ot an ASV 
        label
    count_table : biom.Table
        Once again, a list of 
    """
    # # Sets up the client
    # _setup_dask_client(debug=debug, cluster_config=None,  
    #                    n_workers=n_workers)

    # Checks the manifest
    _check_manifest(manifest)
    
    # Gets the order of the regions
    region_order = \
        manifest.get_column('region-order').to_series().astype(int).to_dict()
    region_names = {i: r for r, i in region_order.items()}

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

    ### Untangles the database to get the unique regional mapping

    # Filters database down to kmers which are present in the sequences 
    # because otherwise we're trying to untangle a huge amount of data and it 
    # just gets memory intensive and slow
    kmer_map = pd.concat(
        axis=0, 
        sort=False,
        objs=_read_manifest_files(manifest, 'kmer-map', 
                                 'FeatureData[KmerMap]', pd.DataFrame)
        )
    kmer_map.reset_index(inplace=True)
    kmer_map['region'] = kmer_map['region'].replace(region_order)
    
    print('Regional Kmers Loaded')

    # Builds database mapping bettween the kmer, the original database
    # sequence and the original name
    db_map = _untangle_database_ids(kmer_map.copy(),
                                    kmers=kmers,
                                    block_size=100,
                                    )
    print('Database map assembled')    

    ### Summarizes the database
    kmer_map['clean_name'] = kmer_map['db-seq'].replace(db_map.to_dict())
    db_summary = _count_mapping(kmer_map.reset_index(), 
                                count_degenerates, 
                                kmer='seq-name')

    ### Constructs the regional alignment
    align_mat = _construct_align_mat(align_map,
                                    sequence_map=db_map.to_dict(),
                                    seq_summary=db_summary,
                                    nucleotide_error=per_nucleotide_error, 
                                    max_mismatch=max_mismatch, 
                                    ).compute()
    print('Alignment map constructed')
    
    ### Solves the relative abundance
    counts = pd.concat(
        axis=1, 
        sort=False, 
        objs=_read_manifest_files(manifest, 'frequency-table', 
                                 'FeatureTable[Frequency]', pd.DataFrame)).T
    counts.fillna(0, inplace=True)

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

    # Performs the maximum liklihood reconstruction on a per-sample basis. 
    # Im not sure if this could be refined to optimize the alogirthm
    # to allow multiple samples ot be solved together, but... eh?
    rel_abund = _solve_iterative_noisy(align_mat=align_mat, 
                                       table=n_table,
                                       min_abund=min_abund,
                                       seq_summary=db_summary,
                                       ).fillna(0)
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


def _build_id_set(tangle, block_size=500):
    """
    Identifies a set of sequences that exist together

    Parameters
    ----------
    tangle: pd.DataFrame
        A dataframe which describes the sequence (index) and the linked
        sequences across all regions which still have partners
    block_size
        The number of sequences to be linked to the analysis block
    Returns
    -------
    list
        list of arrays of groups of sequences which are linked together for 
        untangling

    """
    tangle.sort_values(list(tangle.columns)[::-1], inplace=True)

    def _check_ids(id_):

        last_set = [id_]
        new_set = set([])
        stopper = 0

        while stopper < 2:
            check = tangle.isin(last_set).any(axis=1).values
            new_set = \
                np.sort(tangle.loc[check].melt()['value'].dropna().unique())
            stopper += 1.*(set(new_set) == set(last_set))
            last_set = new_set

        return last_set

    # And now we untangle until we get the blocks of sequences.
    # The while loops are a bit dirty, but the idea is that we dont know how
    # many sequences are going to get linked or how many sets need to build
    # the blocks. A for-loop might cause issues... 
    all_ids = tangle['db-seq'].unique()
    skip = set([])
    id_set = []
    while len(all_ids) > 0:
        block = [] 
        # For each remainig id 
        for id_ in all_ids:
            if (id_ in skip):
                continue
            seq_group = _check_ids(id_)
            skip = skip.union(set(seq_group))
            block = np.hstack([block, seq_group])
            if (len(block) >= block_size):
                break
        id_set.append(block)
        all_ids = all_ids[~np.isin(all_ids, list(skip))]

    return id_set


def _check_manifest(manifest):
    """
    Makes sure that everything is in the manifest

    Parameters
    ---------
    Manifest : qiime2.Metadata
        A manifest file describing the relationship between regions and their
        alignment mapping. The manifest must have at least three columns
        (`kmer-map`, `alignment-map` and `frequency-table`) each of which
        contains a unique filepath. 
    """
    manifest = manifest.to_dataframe()
    cols = manifest.columns
    if not (('kmer-map' in cols) & ('alignment-map' in cols) & 
            ('frequency-table' in cols) & ('region-order' in cols)):
        raise ValidationError('The manifest must contain the columns '
                              'kmer-map, alignment-map and frequency-table.\n'
                              'Please check the manifest and make sure all'
                              ' column names are spelled correctly')
    manifest = manifest[['kmer-map', 'alignment-map', 'frequency-table']]
    if pd.isnull(manifest).any().any():
        raise ValidationError('All regions must have a kmer-map, '
                              'alignment-map and frequency-table. Please '
                              'check and make sure that you have provided '
                              'all the files you need')
    if (len(manifest.values.flatten()) != 
            len(np.unique(manifest.values.flatten()))):
        raise ValidationError('All paths in the manifest must be unique.'
                             ' Please check your filepaths')
    if not np.all([os.path.exists(fp_) for fp_ in manifest.values.flatten()]):
        raise ValidationError('All the paths in the manifest must exist.'
                             ' Please check your filepaths')
    if not np.all([os.path.isfile(fp_) for fp_ in manifest.values.flatten()]):
        raise ValidationError('All the paths in the manifest must be files.'
                             ' Please check your filepaths')


@dask.delayed
def _construct_align_mat(match, sequence_map, seq_summary, 
    nucleotide_error=0.005, max_mismatch=2, kmer_name='kmer', asv_name='asv',
    seq_name='clean_name', miss_col='mismatch'):
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
    align_mat = _expand_duplicate_sequences(match, id_col=kmer_name)
    align_mat['db-seq'] = \
        align_mat[kmer_name].apply(lambda x: x.split('@')[0])
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

    # Calculates the maximum error rate to get the mismatch
    asv_ref['error_thresh'] = \
        (asv_ref['perfect_prob'] - 0.1*asv_ref['max_error_prob'] - 
            np.spacing(1))
    asv_ref['perf_match'] = \
        (asv_ref['match_prob'] > asv_ref['error_thresh'])

    # # ### Computes Mij: Pr(kmer = j | read = j)
    # Normalizes the alignment probability based on the pre-region 
    # amplification
    region_amp = seq_summary['num-regions']
    region_norm = asv_ref[seq_name].replace(region_amp)
    asv_ref['norm'] = (asv_ref['match_prob'] / region_norm)

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


@dask.delayed
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
    long_.reset_index(inplace=True)
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

    overlap['seq_count'] = overlap.loc[check_diff].apply(check_shared, axis=1)
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


def _get_unique_kmers(series):
    kmers = np.hstack([[a.split("@")[0] for a in kmer.split('|')] 
                       for kmer in series])
    return np.sort(np.unique(kmers))


@dask.delayed
def _map_id_set(id_set, tangle):
    """
    Untangles a set of sequence ids back into a comprehensive identifier
    Parmaters
    ---------
    id_set : array-like
        The list of ids to be explored in tangle. Possibly generated using
        `_build_id_set`.
    tangle : DataFrame
        An expanded dataframe mapping each database sequence `db_seq` to
        the region (`region`) and all accompanying sequences
    centroids : set
        The set of possible sequences which may be centroid sequences
    
    Returns
    -------
    DataFrame
        The mapping between the oringal database sequence (`db_seq`) and the 
        name that will appear in the intermediate dataset (`clean_name`) and 
        the number of regions the sequence covers (`num_regions`)
    """

    id_set = np.sort(np.unique(id_set))
    tangle_filter = tangle.drop(columns=['region']).isin(id_set).any(axis=1)
    exp1 = tangle.loc[tangle_filter].copy()
    exp1 = exp1.set_index('db-seq').loc[(id_set)].reset_index()
    exp1.dropna(axis=0, how='all', inplace=True)

    # Minimizes the region in an attempt to make this make more sense for
    # sparse. 
    exp1.replace({r: i for i, r in 
                  enumerate(sorted(exp1['region'].dropna().unique()))},
                  inplace=True)
    num_regions = len(exp1['region'].unique())
    seq_coords = {seq: i for i, seq in enumerate(sorted(id_set))}
    num_seqs = len(id_set)

    # Identifies sequences which may be missing and compares the coverage of
    # thsoe sequences, with the idea that sequences which intersect should 
    # overlap along the full lenght of a sequence, rather than having a shared
    # central sequence, and then missing the sequence on either end. We're
    # assuming the coverage behavior based on the fact that our databases 
    # should represent some level of clustering, and so if the database
    # doesn't cover the region, we still treat the sequence as unqiue
    miss_seq = exp1.drop_duplicates(['region', 'db-seq']).copy()
    miss_seq = miss_seq.pivot(index='db-seq', 
                              columns='region', 
                              values='0').isna()

    miss_seq = miss_seq.loc[id_set]
    cum_coverage = (~miss_seq).cumsum(axis=1) * (~miss_seq)
    cum_coverage = cum_coverage.loc[id_set]

    horiz = np.dstack([cum_coverage.values] * num_seqs).swapaxes(1, 2)
    coverage_overlap = sp.as_coo(
        ((horiz <= horiz.swapaxes(0, 1)).all(axis=2) | 
         (horiz >= horiz.swapaxes(0, 1)).all(axis=2))
        )#
    # Finds where sequences are missing, which allows the possible mapping
    any_missing = sp.as_coo(np.dstack([
        miss_seq[[region]].values | miss_seq[[region]].values.T
        for region in miss_seq.columns
        ]))
    

    # Finds the regions and coordinates where the sequences overlap based 
    # on the shared sequences
    long_mat = exp1.drop_duplicates(['region', 'db-seq']).melt(
        id_vars=['db-seq', 'region'],
        value_vars=exp1.drop(columns=['db-seq', 'region']),
        # value_name='',
        ).dropna().copy()
    long_mat.sort_values(['db-seq', 'region', 'value'], inplace=True)
    long_mat.drop_duplicates(['db-seq', 'region', 'value'], inplace=True)

    # Gets the alignment matrix
    num_coords = len(id_set)
    #  Builds a sparse matrix of the matches
    match_array = sp.COO(
        long_mat[['db-seq', 'value', 'region']].replace(seq_coords).values.T,
        data=True,
        shape=(num_coords, num_coords, num_regions)
    )

    # We find matching coordinates where we look for places that things are
    # missing or matched given that there is at least one hard match between 
    # the sequences in any region

    ### Finds anything that overlaps
    possible_match = (match_array | any_missing).all(axis=2)
    match_coords = possible_match & coverage_overlap

    # ...And we map the matches back to the database sequences
    long_match = pd.DataFrame(match_coords.coords, 
                              index=['db-seq', 'clean_name']).T
    long_match.replace({i: seq for seq, i in seq_coords.items()},
                       inplace=True)
    long_match.rename(columns={0: 'db-seq', 1: 'clean_name'}, inplace=True)

    return long_match


@dask.delayed
def _map_singletons(unique):
    """
    Maps the sequences which do not share any their region with any partners

    Parameters
    ----------
    regional_data: DataFrame
        the data which links sequences (kmers) to regions

    Returns
    -------
    DataFrame
        The mapping between the oringal database sequence (`db_seq`) and the 
        name that will appear in the intermediate dataset (`clean_name`) and 
        the number of regions the sequence covers (`num_regions`)
    """
    unique.drop_duplicates(['db-seq', 'region'], inplace=True)
    seq_map = pd.Series(
        unique['db-seq'].drop_duplicates().values, 
        index=pd.Index(unique['db-seq'].drop_duplicates(), name='db-seq'),
        name='clean_name'
    )
    return seq_map.reset_index()


def _read_manifest_files(manifest, dataset, semantic_type=None, view=None):
    """
    Extracts files from the manifest and turns them into a list of objects
    for analysis
    """
    paths = manifest.get_column(dataset).to_series()
    artifacts = [Artifact.load(path) for path in paths]
    if semantic_type is not None:
        type_check = np.array([str(a.type) == semantic_type 
                               for a in artifacts])
        if not np.all(type_check):
            err_ = '\n'.join([
                'Not all %s Artifacts are of the %s semantic type.' 
                    % (dataset.replace('-', ' '),   semantic_type),
                'Please review semantic types for these regions:',
                '\n'.join(paths.index[type_check == False])
                ])
            raise TypeError(err_)
    if view is not None:
        return [a.view(view) for a in artifacts]
    else:
        return artifacts


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

    # Gets the sequence positions
    align_mat = align_mat.loc[align_mat[seq_name].isin(relative.index)].copy()
    align_seqs = align_mat[seq_name].values
    align_asvs = align_mat[asv_name].values
    samples = relative.columns

    # Calculates the probability of a count showing up in a particular 
    # reference x ASV pair
    p_r_and_j = pd.concat(axis=1, objs=[
        align_mat.set_index(seq_name)[asv_name],
        (relative.loc[align_seqs].T * align_mat['norm'].values).T
        ])
    p_r_and_j.index.set_names(seq_name, inplace=True)
    p_r_and_j.reset_index(inplace=True)
    
    # Calculates the probability of a count mapping to the reference given 
    # that it mapped to the jth ASV
    j_ind = p_r_and_j.groupby('asv')[samples].sum() + np.spacing(1)
    p_r_given_j = \
        (p_r_and_j.set_index([seq_name, asv_name]) /
         j_ind.loc[align_asvs].values)

    # Maps the counts to the sequence given the probability
    counts_r_given_j = p_r_given_j * counts.loc[align_asvs].values
    counts_r_given_j.reset_index(inplace=True)

    # And then sums the counts over the sequence (so we get the total counts)
    counts_r = (counts_r_given_j.groupby(seq_name)[samples].sum().T / 
                seq_summary['num-regions']).T

    return counts_r.round(0)


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
    # And then we get indexing because I like to have the indexing
    # references
    align_seqs = align.columns
    align_asvs = align.index
    align = sp.as_coo(align.values)

    recon = []
    for sample, col_ in table.loc[align_asvs].iteritems():
        abund = dask.delayed(sp.as_coo)(col_.values)
        freq_ = _solve_ml_em_iterative_1_sample(align=align,
                                                abund=abund,
                                                sample=sample,
                                                align_kmers=align_seqs,
                                                num_iter=num_iter,
                                                tolerance=tolerance,
                                                min_abund=min_abund,
                                                )
        df = dask.delayed(pd.DataFrame)(freq_, columns=[sample])
        recon.append(df)
    recon = pd.concat(axis=1, sort=False, objs=dask.compute(*recon))
    recon.columns = table.columns

    # Tidies the table
    def tidy_table(df):
        df = df /  seq_summary.loc[df.index, ['num-regions']].values
        df = df / df.sum(axis=0)
        df.fillna(0, inplace=True)
        df = df.loc[(df > 0).any(axis=1)]
        df.sort_index(axis='index', inplace=True)
        df.index.set_names(seq_name, inplace=True)
        df.sort_index(axis='columns', inplace=True)
        return df

    tidy = tidy_table(recon)

    return tidy


@dask.delayed
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
    bact_freq = sp.dot(align.T, abund) / np.sum(sp.dot(align.T, abund))

    # Solves the mixture. I'm not sure if this can be reimplemented as a
    # gaussian mixed models problem if we have multiple samples? Because like
    # that might actuallky make it easier? Maybe???
    for i in np.arange(0, num_iter):
        ### Expectation
        # assign theta estimate for each bacteria
        theta_i = sp.dot(align, bact_freq)
        
        ### Maximization
        # Adjusts the abundance and bacterial frequency
        r_weighted = abund / (theta_i + np.spacing(1))
        bact_factor = sp.dot(r_weighted.T, align)
        
        # Computes errro 
        error = np.dot(np.absolute(1 - bact_factor.todense()), bact_freq)
        
        # And then we get the new frequency... 
        bact_freq = bact_freq * bact_factor
        
        if error < 1e-7:
            break

        high_enough = bact_freq.todense() > min_abund
        # if sum(~high_enough) > 0:
            # print(np.arange(0, len(bact_freq))[~high_enough])
        align = align[:, high_enough]
        bact_freq = bact_freq[high_enough]
        align_kmers = align_kmers[high_enough]

    # And then we do hard threshholding
    bact_freq = bact_freq.todense()
    bact_freq[bact_freq <= min_abund] = 0
    bact_freq = bact_freq / bact_freq.sum()

    return pd.Series(bact_freq, align_kmers, name=sample)


def _untangle_database_ids(region_db, kmers=None, block_size=2000):
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

    # Gets list of sequences that should be untangled together
    region_db['region'] = region_db['region'].astype(int)
    region_db['region'] = region_db['region'] - region_db['region'].min()
    region_db.sort_values(['db-seq', 'region', 'seq-name'], inplace=True)
    region_db.drop_duplicates(['db-seq', 'region', 'kmer'], inplace=True)

    if kmers is None:
        kmers = region_db['db-seq'].unique()

    def _split_seq(x):
        return pd.Series(np.unique([y.split('@')[0] for y in x.split('|')]))

    # Gets the sequences that are shared across the kmers in the regions
    shared = \
        region_db.set_index(['db-seq', 'region'])['kmer'].apply(_split_seq)
    shared.reset_index(inplace=True)
    shared = shared.melt(id_vars=['db-seq', 'region']).dropna()
    shared.sort_values(['db-seq', 'region', 'value'], inplace=True)
    shared.drop_duplicates(['db-seq', 'region', 'value'], inplace=True)
    shared['variable'] = 1.
    shared['variable'] = \
        shared.groupby(['db-seq', 'region'])['variable'].cumsum() - 1
    shared['variable'] = shared['variable'].astype(int).astype(str)
    shared = shared.set_index(['db-seq', 'region', 'variable']).unstack()
    shared.columns = shared.columns.droplevel(None)
    shared.reset_index(inplace=True)

    in_kmers = shared.drop(columns=['region']).isin(kmers).any(axis=1)
    shared = shared.loc[in_kmers].copy()

    # Finds all the sequences that are linked across all regions and checks
    # sequences converage
    def find_partners(g):
        shared_seqs = g.values.flatten()[~pd.isnull(g.values.flatten())]
        return pd.Series(np.sort(np.unique(shared_seqs)))
    partners = \
        shared.drop(columns=['region']).groupby('db-seq').apply(find_partners)
    partners = partners.unstack()
    num_seqs = len(partners)

    print('finding all shared sequences. There are %i sequences.' % num_seqs)

    # Splits out sequences that are never linked to other sequences. This way,
    # we dont spend extra time findling sequences, etc.
    share_coverage = \
        ~partners.drop(columns=[0]).isin(partners.index.values).any(axis=1)

    # Handles the sequences which do not share any regions with anything else
    unique_seqs = share_coverage.index[share_coverage]
    unique = region_db.loc[region_db['db-seq'].isin(unique_seqs)].copy()
    unique_map = _map_singletons(unique)
    
    # Then, we need to deal with sequences that *are* linked to others. 
    still_tangled = partners.loc[~share_coverage].copy()
    print('found unique sequences. %i/%i seqs mapped, %i/%i remaining' 
          % (len(unique_seqs), num_seqs, len(still_tangled), num_seqs))

    seq_blocks = _build_id_set(still_tangled.reset_index(), 
                               block_size=block_size)

    # Builds the new sequence map
    seq_map = [_map_id_set(id_set, shared.loc[shared['db-seq'].isin(id_set)])
               for id_set in seq_blocks]
    seq_map.append(unique_map)

    # And now we get hte inital mapping to all the sequences.
    initial_map = pd.concat(axis=0, objs=dask.compute(*seq_map))

    # Tides up the mapping into something we'll use later to make pivoting easier
    initial_map.sort_values(['db-seq', 'clean_name'], inplace=True)
    initial_map['counter'] = 1
    initial_map['counter'] = \
        initial_map.groupby('db-seq')['counter'].cumsum() - 1
    initial_names = \
        initial_map.groupby('db-seq')['clean_name'].apply(lambda x: "|".join(x))

    # Checks for correct associations. We assume there should be the same 
    # number of instances of a sequence as there are sequences in a cluster
    count_check = initial_names.value_counts().reset_index()
    count_check['seq_count'] = \
        count_check['index'].apply(lambda x: x.count("|") + 1)
    count_check['okay'] = \
        (count_check['seq_count'] == count_check['clean_name'])

    @dask.delayed
    def _pass_okay_clusters(okay_seqs):
        return okay_seqs

    if count_check['okay'].all():
        seq_names = initial_names
    else:
        okay_ids = count_check.loc[count_check['okay'], 'index'].values
        okay_seqs = initial_names[initial_names.isin(okay_ids)].copy()
        tangled_seqs = \
            initial_names.index[~initial_names.isin(okay_ids)].copy()

        print('found initial sequence blocks. %i/%i seqs mapped, '
              '%i/%i remaining'  % (num_seqs - len(tangled_seqs), 
                                    num_seqs, len(tangled_seqs), num_seqs))
    
        more_tangles = \
            initial_map.loc[initial_map.isin(tangled_seqs).any(axis=1)].copy()
        wide_ = more_tangles.pivot(index='db-seq', 
                                   columns='counter', 
                                   values='clean_name')

        # Repartitions the entangled sequences and then untangles them
        recluster = _build_id_set(wide_.reset_index(), block_size)
        plaiting = [_detangle_names(
            initial_map.loc[initial_map['db-seq'].isin(cluster)].reset_index()
            ) for cluster in recluster]
        plaiting.append(_pass_okay_clusters(okay_seqs))
        seq_names = pd.concat(dask.compute(*plaiting))


    seq_names.sort_values(inplace=True)

    print('found final sequence blocks. %i/%i seqs mapped,0/%i remaining' 
          % (len(seq_names), num_seqs, num_seqs))

    return seq_names

