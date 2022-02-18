import copy
import itertools as it
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import biom
import dask
from dask.delayed import Delayed
import dask.dataframe as dd
import numpy as np
import pandas as pd

from qiime2 import Metadata, Artifact
from qiime2.plugin import ValidationError
from q2_sidle._utils import (_setup_dask_client, 
                             degen_reps,
                             _check_regions,
                             )
from q2_sidle._build_database import _check_regions

def reconstruct_counts(
    region: str,
    regional_alignment: Delayed,
    regional_table: biom.Table,
    database_map: pd.Series,
    database_summary: Metadata,
    per_nucleotide_error: float=0.005,
    min_abund: float=1e-5,
    region_normalize: str='average',
    block_size : int=10000,
    min_counts: int=1000,
    debug: bool=False, 
    n_workers: int=1,
    client_address: str=None,
    ) -> (biom.Table):
    """
    Reconstructs regional alignments into a full length 16s sequence

    Parameters
    ----------
    # region: 
    count_degenerates: bool
        Whether degenerate sequences should be counted as unique kmers during
        reconstruction or only unique sequences should be counted.
    per_nucleotide_error: float
        The assumed per-nucleotide error rate from amplifican and sequencing.
        The value used here is based on benchmarking from the original smurf
        paper
    min_abund: float
        The minimum relative abundance for a feature to be retained during 
        reconstruction. The default from the original smurf algorithm is '
          '1e-10, the number here may depend on the total number of '
          'sequences in your sample.
    region_normalize: str ('average', 'weighted', 'unweighted')
        Controls the behavior of count normalization by region. When 
        `region_normalize="average"`, the depth will be the average depth
        over all regions. When "weight" is used, the depth will correspond to
        all sequneces with a coverage normalization factor. So, there are 2 regions
        with 1000 sequences each, "average" will give a total count of 1000 and 
        "weight" will give a total count of 2000. "unweighted" will ignore
        the regional normalization in the abundance calculation and is preferable
        for meta analysis.
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

    region_order, region_names, num_regions = _check_regions(region)

    ### Constructs the regional alignment map
    database_summary = database_summary.to_dataframe()
    database_summary.index.set_names('clean_name', inplace=True)
    if region_normalize == 'unweighted':
        database_summary['num-regions'] = 1

    align_mat = _construct_align_mat(matches=regional_alignment,
                                     sequence_map=database_map,
                                     seq_summary=database_summary,
                                     nucleotide_error=per_nucleotide_error, 
                                     blocksize=block_size,
                                     )
    print('Alignment map constructed')
    
    ### Solves the relative abundance
    counts = regional_table[0]
    if len(regional_table) > 1:
        for table_ in regional_table[1:]:
            counts = counts.merge(table_)

    counts = pd.DataFrame(
        counts.matrix_data.toarray(),
        index=counts.ids(axis='observation'),
        columns=counts.ids(axis='sample'),
    )
    counts.fillna(0, inplace=True)

    # We have to account for the fact that some of hte ASVs may have been 
    # discarded because they didn't meet the match parameters we've set or 
    # because they're not in the database.
    keep_asvs = list(set(align_mat['asv'].values) & (set(counts.index)))
    unaligned_counts = counts.copy().drop(keep_asvs).sum(axis=1)
    counts = counts.loc[keep_asvs]
    keep_samples = counts.sum(axis=0) > min_counts
    if keep_samples.all() == False:
        raise ValueError('There are {samples:d} samples with fewer than '
                         '{depth:d} total sequences. Please check your '
                         'minimum counts and make sure your representative '
                         'sequences are aligned with the database.'.format(
                            samples=(keep_samples == False).sum(), 
                            depth=min_counts)
                         )
    # elif not keep_samples.all():
    #     warnings.warn("There are %i samples with fewer than %i total"
    #                   " reads. These samples will be discarded."
    #                   % ((keep_samples==False).sum(), min_counts),
    #                   UserWarning)
    counts = counts[keep_samples.index.values[keep_samples.values]]
    counts = counts.loc[counts.sum(axis=1) > 0]
    keep_asvs = counts.index

    align_mat = align_mat.loc[align_mat['asv'].isin(keep_asvs)]

    # Normalizes the alignment table
    n_table = counts / counts.sum(axis=0)
    print('counts loaded')

    # Performs the maximum liklihood reconstruction on a per-sample basis. 
    # Im not sure if this could be refined to optimize the alogirthm
    # to allow multiple samples ot be solved together, but... eh?
    sample_ids = n_table.columns.values
    n_samples = len(sample_ids)
    print('start normalization')

    rel_abund = _solve_iterative_noisy(align_mat=align_mat, 
                                       table=n_table,
                                       min_abund=min_abund,
                                       seq_summary=database_summary)
    database_summary = database_summary.loc[rel_abund.ids(axis='observation')]
    print('Relative abundance calculated')

    # Puts together the regional normalized counts
    count_table = _scale_relative_abundance(align_mat=align_mat,
                                            relative=rel_abund,
                                            counts=counts,
                                            region_normalize=region_normalize,
                                            num_regions=num_regions,
                                            seq_summary=database_summary)
    count_table = count_table.filter(lambda v, id_,  md: v.sum() > 0,  
                                     axis='observation')
    return count_table


def _construct_align_mat(matches, sequence_map, seq_summary, 
    nucleotide_error=0.005, kmer_name='kmer', asv_name='asv',
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

    # Prepares the delayed align maps by filtering, replacing the region order
    # and addding the database and clean names
    cleaned_matches = []
    for match in matches:
        expanded = dask.delayed(_expand_duplicate_sequences)(match, kmer_name)
        cleaned = \
            dask.delayed(_get_db_and_clean)(expanded, sequence_map, kmer_name, seq_name)
        cleaned_matches.append(cleaned)

    align_mat = dd.from_delayed(
        cleaned_matches,
        meta=[('db-seq', 'str'), (kmer_name, 'str'), (asv_name, 'str'), 
              ('length', int),  (miss_col, int), ('max-mismatch', int), 
              ('region', object), (seq_name, str)]
        )
    n_aligned = len(align_mat)
    align_mat = align_mat.repartition(int(np.ceil(n_aligned / blocksize)))
    align_mat.drop_duplicates([seq_name, asv_name, 'region'], inplace=True)

    # Collapses into ASV-ref seq combos where for degenerates, the closet
    # match is retained as the sequence.
    asv_ref = align_mat.groupby([seq_name, asv_name, 'region', 'max-mismatch'], 
                                sort=False)
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
        np.power(1 - nucleotide_error, asv_ref['length'] - asv_ref['max-mismatch']) * 
        np.power(nucleotide_error / 3, asv_ref['max-mismatch'])
    )

    # Calculates the maximum error rate to get the mismatch. 0.1 multiplier comes
    # directly from the SMURF code
    # https://github.com/NoamShental/SMURF/blob/master/mFiles/solve_iterative_noisy.m#L53
    asv_ref['error_thresh'] = \
        (asv_ref['perfect_prob'] - 0.1*asv_ref['max_error_prob'] - 
            np.spacing(1))
    asv_ref['perf_match'] = \
        (asv_ref['match_prob'] > asv_ref['error_thresh'])

    ### Computes Mij: Pr(kmer = j | read = j)
    # Normalizes the alignment probability based on the pre-region 
    # amplification
    region_amp = seq_summary['num-regions']

    asv_ref = asv_ref.set_index(seq_name)
    asv_ref['region-norm'] = region_amp.astype(int)

    asv_ref['norm'] = (asv_ref['match_prob'] / asv_ref['region-norm'])

    asv_ref  = asv_ref.reset_index().compute()
    asv_ref.drop(columns='region-norm', inplace=True)

    return asv_ref


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


def _get_db_and_clean(df, clean, kmer_name, seq_name):
    """
    A delayaed function to extract a database sequence and map it to a 
    cleaned database name
    """
    df['db-seq'] = df[kmer_name].apply(lambda x: x.split("@")[0])
    df.set_index('db-seq', inplace=True)
    df[seq_name] = clean.loc[df.index]
    return df.reset_index()


def _scale_relative_abundance(align_mat, relative, counts, seq_summary,
    num_regions, region_normalize='average', seq_name='clean_name', 
    asv_name='asv'):
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
    align_mat = align_mat.loc[
        align_mat[seq_name].isin(relative.ids(axis='observation')) & 
        align_mat[asv_name].isin(counts.index)
        ]
    align = align_mat.pivot_table(values='norm', index=asv_name, 
                                  columns=seq_name, fill_value=0
                                  )

    # And then we get indexing because I like to have the indexing
    # references
    align = align[relative.ids(axis='observation')]
    align_seqs = align.columns
    align_asvs = align.index
    align = align.values

    counts = counts.loc[align_asvs]
    samples = list(relative.ids(axis='sample'))
    spacer =  np.spacing(1)

    # A simple per-sample function to solve the counts because frequency is fun
    # and pandas rounding may be what's killing this
    @dask.delayed
    def _solve_sample(align, freq, count, align_seqs, sample):
        # Probability a count shows up in a particular ASV/reference pairing
        p_r_and_j = freq * align
        # Normalized for the number of asvs mapped to that reference
        p_r_given_j = \
            p_r_and_j / np.atleast_2d(p_r_and_j.sum(axis=1) + spacer).T
        # sums the counts
        count_j_given_r = count.values * p_r_given_j



        return biom.Table(
            np.atleast_2d(count_j_given_r.sum(axis=0)).T,
            sample_ids=[sample], observation_ids=align_seqs,
            )
    scaled_counts = []
    for sample in samples:
        non_zero_asvs = (counts[sample] > 0).values
        non_zero_seqs = (relative.data(sample) > 0)
        relative_seqs = relative.copy().filter(align_seqs[non_zero_seqs], 
                                               axis='observation')

        counted = _solve_sample(
            align[non_zero_asvs][:, non_zero_seqs],
            relative_seqs.data(sample),
            counts.loc[non_zero_asvs, [sample]],
            align_seqs[non_zero_seqs],
            sample
            )
        scaled_counts.append(counted)

    @dask.delayed
    def _combine_tables(*tables):
        tables = list(tables)
        if len(tables) == 1:
            return tables[0]
        else:
            return tables[0].concat(tables[1:], axis='sample')

    scaled_counts = _combine_tables(*scaled_counts).compute()
    scaled_counts.add_metadata(
        seq_summary[['num-regions']].to_dict(orient='index'),
        axis='observation')
    if region_normalize == 'average':
        scaled_counts.transform(
            lambda v, id_, md: np.round(v / md['num-regions'], 0),
            axis='observation',
            inplace=True
            )
    elif region_normalize == 'weighted':
        scaled_counts.transform(
            lambda v, id_, md: np.round(v * num_regions / md['num-regions'], 0),
            axis='observation',
            inplace=True
            )
    else:
        scaled_counts.transform(
            lambda v, id_, md: np.round(v, 0),
            axis='observation',
            inplace=True
            )

    return scaled_counts


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
    align_seqs = align.columns.values
    align_asvs = align.index.values

    recon = []
    for sample, col_ in table.iteritems():
        filt_align = align.loc[col_ > 0].values
        abund = col_[col_ > 0].values
        sub_seqs = copy.copy(align_seqs)[(filt_align > 0).any(axis=0)]
        filt_align = filt_align[:, (filt_align > 0).sum(axis=0) > 0]

        freq_ = dask.delayed(_solve_ml_em_iterative_1_sample)(
            align=filt_align,
            abund=abund,
            sample=sample,
            align_kmers=sub_seqs,
            num_iter=num_iter,
            tolerance=tolerance,
            min_abund=min_abund,
            )

        recon.append(freq_)
    recon = dask.compute(*recon)
    @dask.delayed
    def _combine_tables(*tables):
        tables = list(tables)
        if len(tables) == 1:
            return tables[0]
        else:
            return tables[0].concat(tables[1:], axis='sample')

    recon = _combine_tables(*recon)
    recon, = dask.compute(recon)

    recon.add_metadata(
        seq_summary.loc[recon.ids(axis='observation'), 
                       ['num-regions']].to_dict(orient='index'),
        axis='observation'
        )
    recon.transform(lambda v, id_, md: v / md['num-regions'],
                         axis='observation',
                         inplace=True)
    recon.norm(axis='sample', inplace=True)
    

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
    align : ndarray
        The alignment matrix describing  the probability a given ASV belongs 
        to a reference sequence
    abund : 1D-ndarray
        The relative  abundance of each ASV feature assigned to each region
    align_kmers: ndarray
        The names of the reference sequences
    sample : str
        The name of the sample
    n_iter: int
        the number of iterations to try before giving up
    tol : float
        The error tolerance for solving the data in optimization
    min_abund : float
        The minimum relative abundance  to retain  a  feature
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
        try:
            align_kmers = align_kmers[high_enough]
        except:
            print('Could not align!', sample)

    # And then we do hard threshholding
    bact_freq[bact_freq <= min_abund] = 0
    bact_freq = bact_freq / bact_freq.sum()

    return biom.Table(np.atleast_2d(bact_freq).T, 
                      observation_ids=align_kmers, 
                      sample_ids=[sample])


