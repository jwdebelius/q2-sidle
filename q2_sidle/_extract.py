import itertools as it
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)

import dask
import numpy as np
import pandas as pd
import regex

from qiime2 import Metadata
from q2_types.feature_data import (DNAFASTAFormat, DNAIterator)
from q2_sidle._utils import (_setup_dask_client, 
                             _convert_generator_to_seq_block,
                             _convert_generator_to_delayed_seq_block, 
                             _convert_seq_block_to_dna_fasta_format,
                             _count_degenerates,
                             )

complement = {'A': 'T', 'T': 'A', 'G': 'C',  'C': 'G',
              'R': 'Y', 'Y': 'R', 'S': 'S',  'W': 'W',
              'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
              'D': 'H', 'H': 'D', 'N': 'N'}
degen_sub = {'R': 'AG',
             'Y': 'CT',
             'S': 'CG',
             'W': 'AT',
             'K': 'GT',
             'M': 'AC',
             'B': 'CGT',
             'D': 'AGT',
             'H': 'ACT',
             'V': 'ACG',
             'N': 'ACGT',
             'A': np.nan,
             'C': np.nan,
             'G': np.nan,
             'T': np.nan,
             }


def prepare_extracted_region(sequences: DNAFASTAFormat, region:str, 
    trim_length:int, chunk_size:int=10000, 
    debug:bool=False, n_workers:int=0
    ) -> (DNAFASTAFormat, pd.DataFrame):
    """
    Prepares and extracted database for regional alignment

    This function takes an amplified region of the database, expands the
    degenerate sequences and collapses the duplciated sequences under a 
    single id that can be untangled later.

    Parameters
    ----------
    sequences: q2_type.DNAFASTAFormat
        The regional sequences to be collapsed
    region: str
        A unique name for the region being handled
    trim_length : int
        The length of final sequences to matched the trimmed kmers for 
        kmer-based alignment.
    chunk_size: int, optional
        The number of sequences to group for analysis
    debug: bool
        Whether the function should be run in debug mode (without a client)
        or not. `debug` superceeds all options
    n_workers: int, optional
        The number of jobs to initiate. When `n_workers` is 0, the cluster 
        will be able to access all avaliable resources.

    Returns
    -------
    q2_types.DNAFASTAFormat
        The reads with degenerate nucleotides expanded and duplicated 
        sequences collapsed.
    DataFrame
        Output location for the mapping between the kmer sequence name and the
        the full database sequence name.
    """

    # Sets up the client
    _setup_dask_client(debug=debug, cluster_config=None,  
                       n_workers=n_workers)
    # Reads in the sequences
    sequences = sequences.view(DNAIterator)
    seq_block = _convert_generator_to_delayed_seq_block(sequences, chunk_size)

    # Collapses the sequences
    collapse, map_ = _tidy_sequences(seq_block, region, trim_length)

    # Converts the sequences back to a nicely behaved qiime artifact
    ff = _convert_seq_block_to_dna_fasta_format([collapse])

    return ff, map_


def extract_regional_database(sequences: DNAFASTAFormat, 
    fwd_primer: str,
    rev_primer: str,
    trim_length:int=-1, 
    region:str=None, 
    primer_mismatch:int=2,
    trim_primers:bool=True,
    reverse_complement_rev:bool=True,
    chunk_size:int=10000, 
    debug:bool=False, 
    n_workers:int=1
    ) -> (DNAFASTAFormat,
         pd.DataFrame):
    """
    Performs in silico PCR to extract regions the database and collapse them
    to be used with sidle

    Parameters
    ----------
    sequences
        The full length database sequences to be extracted
    primer_fwd, primer_rev: str
        The forward and reverse primers to be amplified
    trim_length : int, optional
        The length sequences should be trimmed to match the asvs for '
        'kmer-based alignment.
    region : str
        The hypervariable region to profile.
    primer_mismatch: int, optional
        The number of mismatches between the primer and sequence allowed.
    trim_primers: bool, optional
        Whether the primers should be trimmed
    reverse_complement_rev:bool, optional
        When true, the reverse primer will be reverse complimented before
        extraction.

    Returns
    -------
    q2_types.DNAFASTAFormat
        The reads with degenerate nucleotides expanded and duplicated 
        sequences collapsed.
    DataFrame
        Output location for the mapping between the kmer sequence name and the
        the full database sequence name.
    """
    _setup_dask_client(debug=debug, cluster_config=None,  
                       n_workers=n_workers)
    # Reads in the sequences
    sequences = sequences.view(DNAIterator)
    seq_block = _convert_generator_to_seq_block(sequences, chunk_size)

    # Extracts the regional sequence
    regional_seqs = []
    for seq_ in seq_block:
        regional = _extract_region(
            full_seqs=seq_,
            primer_fwd=fwd_primer, 
            primer_rev=rev_primer, 
            mismatch=primer_mismatch, 
            trim_len=trim_length, 
            trim_fwd=trim_primers, 
            trim_rev=trim_primers, 
            reverse_complement_rev=reverse_complement_rev,
            )
        regional_seqs.append(regional)

    if region is None:
        region = '%s-%s' % (fwd_primer, rev_primer)
    collapse, map_ = _tidy_sequences(regional_seqs, region)

    ff = _convert_seq_block_to_dna_fasta_format([collapse])

    return ff, map_


@dask.delayed
def _extract_region(full_seqs, primer_fwd, primer_rev, mismatch=2, 
    trim_len=None, trim_fwd=True, trim_rev=True, 
    reverse_complement_rev=True):
    """
    Extracts a primer region from a full sequence

    Parameters
    ----------
    full_seqs: DataFrame
        An array of sequences where each row is a sequence (mapped to a
        sequence ID) and each column is a basepair. Missing values for 
        varying sequence lengths can be represented by a nan
    primer_fwd, primer_rev: str
        The forward and reverse primers to be amplified
    mismatch: int, optional
        The number of mismatches between the primer and sequence allowed.
    trim_len: None, 'min', int
        The length of the sequence to be retained. If the trim_len is `None`,
        then no trimming will be performed. `"min"` will result in sequences
        being trimmed to match the shortest sequence between the two primers.
        Otherwise, sequences will be trimmed to the length specified by the
        value.
        If the new trim is longer than the 
    trim_primers: bool, optional
        Whether the primers should be trimmed
    verbose: bool, optional
        Whether the progress of the data should be displayed or not.

    Returns
    -------
    DataFrame
        A sequence array as a dataframe containing the trimmed region

    """
    id_ = full_seqs.index.name
    if id_ is None:
        id_ = 'index'

    # This is so ugly, but it gets the reverse compliment of the primer
    # so we can pass it to the align_ref_seq function. There's probably
    # a far smarter way to do this. 
    if reverse_complement_rev:
        primer_rev = _reverse_complement(
            pd.DataFrame(list(primer_rev)).T
            ).apply(lambda x: ''.join(x), axis=1)[0]

    fwd_ = '((?<=(%s)))((%s){e<=%s})' % (primer_fwd[0], primer_fwd[1:], '%i')
    rev_ = '((?<=(%s)))((%s){e<=%s})' % (primer_rev[0], primer_rev[1:], '%i')
    for bp, rep in degen_sub.items():
        if bp in {'A', 'T', 'G', 'C'}:
            continue
        fwd_ = fwd_.replace(bp, '[%s]' % rep)
        rev_ = rev_.replace(bp, '[%s]' % rep)

    # Determines the read lengths
    read_length = (~full_seqs.isna()).sum(axis=1)
    read_length.name = 'read_length'

    # Determines the sequence positions
    seq_collapse = full_seqs.fillna('').apply(lambda x: ''.join(x), axis=1)
    
    ### Gets the foward
    if trim_fwd:
        fwd_func = _find_primer_end
    else:
        fwd_func = _find_primer_start


    ### Finds the reverse position
    if trim_rev:
        rev_func =_find_primer_start
    else:
        rev_func = _find_primer_end

    # Finds all possible match position and then finds the first position
    # with the smallest mismatch. We'll do this going backward: if we can't
    # find a first match with two mismatches, there's no reason to search the 
    # space again with more percison. 
    seqs_to_searchf = seq_collapse.copy()
    seqs_to_searchr = seq_collapse.copy()
    fwd_match = []
    rev_match = []
    for i in np.arange(mismatch, -1, -1):
        fwd_rough = seqs_to_searchf.apply(fwd_func, primer=(fwd_ % i), 
                                          prefix='fwd_')
        search_gate = ((fwd_rough['fwd_mis'] <= i) & (fwd_rough['fwd_mis'] > 0))
        search_ids_f = fwd_rough.loc[search_gate].index
        fwd_match.append(fwd_rough)
        seqs_to_searchf = seqs_to_searchf.loc[search_ids_f]

        rev_rough = seqs_to_searchr.apply(rev_func, primer=(rev_ % i), 
                                           prefix='rev_')
        search_gate = ((rev_rough['rev_mis'] <= i) & (rev_rough['rev_mis'] > 0))
        search_ids_r = rev_rough.loc[search_gate].index
        rev_match.append(rev_rough)
        seqs_to_searchr = seqs_to_searchr.loc[search_ids_r]
        

    fwd_match = pd.concat(axis=0, objs=fwd_match)
    fwd_match.reset_index(inplace=True)
    fwd_match.sort_values([id_, 'fwd_mis', 'fwd_pos'], ascending=True, 
                          inplace=True)
    fwd_match.drop_duplicates([id_], keep='first', inplace=True)
    fwd_match.rename(columns={id_: 'kmer'}, inplace=True)
    fwd_match.set_index('kmer', inplace=True)

    rev_match = pd.concat(objs=rev_match)
    rev_match.reset_index(inplace=True)
    rev_match.sort_values([id_, 'rev_mis', 'rev_pos'], ascending=True, 
                          inplace=True)
    rev_match.drop_duplicates([id_], keep='first', inplace=True)
    rev_match.rename(columns={id_: 'kmer'}, inplace=True)
    rev_match.set_index('kmer', inplace=True)
    rev_match['rev_pos'] = rev_match['rev_pos'] - 1

    ### Sets up positions for trimming

    ## Matches forward and reverse regions
    # Finds sequence mapped to forward and reverse primers
    fwd_and_rev = pd.concat(axis=1, sort=False, objs=[
        fwd_match, rev_match, read_length
        ])

    # determines the trim length and final position
    if trim_len == 'min':
        trim_len = (fwd_and_rev['rev_pos'] - fwd_and_rev['fwd_pos']).min()

    if trim_len is not None:
        fwd_and_rev['new_rev_pos'] = fwd_and_rev['fwd_pos'] + (trim_len - 1)
    else:
        fwd_and_rev['new_rev_pos'] = fwd_and_rev['rev_pos']

    # And then we filter out any read that can't amplify
    fwd_and_rev.loc[fwd_and_rev['read_length' ] < fwd_and_rev['new_rev_pos'],
                    'new_rev_pos'] = np.nan
    fwd_and_rev.dropna(inplace=True, how='any')

    ### Extracts reads
    trimmed = _trim_masked(full_seqs.loc[fwd_and_rev.index], 
                           fwd_and_rev['fwd_pos'], 
                           fwd_and_rev['new_rev_pos'])

    return trimmed.dropna(how='any')


def _find_primer_end(seq_, primer, prefix=''):
    """
    Finds the last position of a primer sequence
    """
    match = regex.search(primer, seq_)
    if match is not None:
        return pd.Series({'%spos' % prefix: match.end(),
                          '%smis' % prefix: sum(match.fuzzy_counts)})
    else:
        return pd.Series({'%spos' % prefix: np.nan,
                          '%smis' % prefix: np.nan})


def _find_primer_start(seq_, primer, adj=1, prefix=''):
    """
    Finds the first position of a primer sequence
    """
    match = regex.search(primer, seq_)
    if match is not None:
        return pd.Series({'%spos' % prefix: match.start() - adj,
                          '%smis' % prefix: sum(match.fuzzy_counts)})
    else:
        return pd.Series({'%spos' % prefix: np.nan,
                          '%smis' % prefix: np.nan})


def _tidy_sequences(seqs, region, trim_length=None):
    """
    Expands and collapses sequences to give a final summarized set

    Parameters
    ----------
    seqs : list of DataFrames or delayed DataFrames
        A dataframe of size n x m where each of n rows is a sequence and 
        each column (m columns) is a basepair in that sequence.
    region : str
        The hypervariable region to profile.
    trim_length : int, optional
        The length sequences should be trimmed to match the asvs for '
        'kmer-based alignment.

    Returns
    -------
    pd.DataFrame
        A dataframe of size n x m where each of n rows is a sequence and 
        each column (m columns) is a basepair in that sequence without 
        duplicated sequences
    DataFrame
        A mapping betweent the `db_seq`, its associated kmers, and the
        individual sequence name, and the region.
    """


    if trim_length is None:
        # @dask.delayed
        def _trim_seqs(x):
            return x
    else:
        @dask.delayed
        def _trim_seqs(x):
            return x[np.arange(0, trim_length)]
    shrunk = []
    for seq in seqs:
        trimmed = _trim_seqs(seq)
        expanded = _expand_degenerates(trimmed)
        contract = _collapse_duplicates(expanded)
        shrunk.append(contract)

    imploded = _collapse_duplicates(*shrunk).compute()

    map_ = _build_id_map(imploded.index, region)

    return imploded, map_


@dask.delayed
def _collapse_duplicates(*seq_array):
    """
    Collapses a sequexnce array including without degenerate nucleotides

    This will collapse a sequence array based on string matching. Its not
    set up to handle degenerates, but its a quicker function than the 
    version which handles degenerates since it avoids a hamming distance
    calculation which can be expensive. 

    Paramaters
    ----------
    seq_array: DataFrame
        A dataframe of size n x m where each of n rows is a sequence and 
        each column (m columns) is a basepair in that sequence.
    verbose: bool, optional
        If a progress message should be displayed

    Returns
    -------
    pd.DataFrame
        A dataframe of size n x m where each of n rows is a sequence and 
        each column (m columns) is a basepair in that sequence without 
        duplicated sequences

    """

    seq_array = pd.concat(seq_array)

    if seq_array.index.name is None:
        seq_name = 'index'
    else:
        seq_name = seq_array.index.name
    _, seq_len = seq_array.shape
    
    # seq_array.index.name = seq_name
    seq_array = seq_array.fillna('')

    # Converts the sequence array to a string and finds duplicates
    seqs = seq_array.apply(lambda x: ''.join(x), axis=1).reset_index()
    seq_array.tail()
    def rep_f(x):
        if x == seq_name: 
            return 'seq_name' 
        else:
            return 'sequence'
    
    groups = seqs.groupby(0)[seq_name].apply(
        lambda x: '|'.join(sorted(x.values))
        ).reset_index()
    groups = groups.rename(columns=rep_f)
    groups = groups.set_index('seq_name')
    # Prevents resource warning about open files with groupby? We hope
    groups.tail()
    group2 = groups.apply(lambda x: pd.Series(list(x['sequence'])), axis=1)

    return group2.sort_index()


@dask.delayed
def _expand_degenerates(seq_array, pad=4):
    """
    Expands a sequence array with degenerate nucleotides

    Parameters
    ----------
    seq_array: Dataframe
        A dataframe where each row is a sequence and each column is a basepair
        in that sequence. Can contain degenerate (non-canonical) basepairs

    pads: int
        The number of zeros to use to pad the identifier

    Returns
    -------
    DataFrame
        A dataframe where each row is a sequence and each column is a basepair
        in that sequence. Will not contain non-canonical basepairs
    """
    degen_count = _count_degenerates(seq_array)
    if (degen_count == 0).all():
        seq_array.index = seq_array.index.astype(str)
        return seq_array

    # Gets only the degenerate sequences and pulls off the sequence index
    yes_degen = seq_array.loc[degen_count > 0].copy()
    index_name = seq_array.index.name
    if index_name == None:
        index_name = 'index'
    yes_degen = yes_degen.reset_index()


    # Metls to the get the ordered degenate nucleotides
    degen_pos = yes_degen.melt(id_vars=index_name,
                               var_name='nt_pos',
                               value_name='degen_nt', 
                               )
    degen_pos = degen_pos.set_index(index_name)
    degen_pos['degen_rep'] = degen_pos['degen_nt'].replace(degen_sub)
    degen_pos.dropna(inplace=True)
    degen_pos['degen_count'] = 1
    degen_pos['degen_count'] = degen_pos.groupby(index_name).cumsum() - 1
    degen_pos['degen_rep'] = \
        degen_pos['degen_rep'].apply(lambda x: tuple(list(x)))
    degen_pos.reset_index(inplace=True)
    degen_lists = degen_pos.pivot(index=index_name, 
                                  columns='degen_count', 
                                  values='degen_rep')

    # Casts the sequence to a string and sets up the degenerates for 
    # replacement
    tmp_degen_rep = {'R': '%s', 'Y': '%s', 'S': '%s', 'W': '%s', 'K': '%s',
                     'M': '%s', 'B': '%s', 'H': '%s', 'V': '%s', 'N': "%s", 
                     'D': '%s',}
    yes_degen.set_index(index_name, inplace=True)
    yes_degen.replace(tmp_degen_rep, inplace=True)

    # We cast to string to make expansion easier
    yes_str = yes_degen.fillna('').apply(lambda x: ''.join(x), axis=1)
    # Performs the expansion. Because fo the way bags get cast, we 
    # need to add the ids later. 
    nums = np.arange(1, 1000).astype(str)
    def expand_degen_str(id_, seq_):
        reps = degen_lists.loc[id_].dropna().values
        if len(reps) == 1:
            seqs = pd.Series({'%s@%s' % (id_, i.zfill(pad)) : seq_ % a 
                              for i, a in zip(*(nums, reps[0]))})
        else:
            seqs = pd.Series({'%s@%s' % (id_, i.zfill(pad)) : seq_ % a
                              for i, a in zip(*(nums, it.product(*reps)))})
        return seqs

     # Expands the sequences
    degen_expand_str = pd.concat([expand_degen_str(id_, seq_) 
                                  for id_, seq_ in yes_str.iteritems()])

    # Stacks the expanded sequence array
    full_expand = pd.concat([
        seq_array.loc[degen_count == 0],
        degen_expand_str.apply(lambda x: pd.Series(list(x)))
    ])
    full_expand = full_expand.T.reset_index(drop=True).T

    # Returns the expanded array
    return full_expand


def _reverse_complement(seq_array):
    """
    Reverse complements a sequence array

    Parameters
    ----------
    seq_array: Dataframe
        A dataframe of size n x m where each of n rows is a sequence and 
        each column (m columns) is a basepair in that sequence. Can contain 
        degenerates.

    Returns
    -------
    DataFrame
        Reverse complemenet of `seq_array`
    """

    # Gets the complement
    comp_seq = seq_array.replace(complement)
    
    #  Filps the complement direction 
    rc = comp_seq.rename(
        columns={i: j for i, j in enumerate(comp_seq.columns[::-1])}
        )
    rc = rc[np.arange(len(comp_seq.columns))]

    return rc


def _trim_masked(seqs, start, end):
    """
    Extracts the sequence between `start` and `end`

    Paramters
    ---------
    seqs: DataFrame
        An n x m array of sequences where each row is a sequence (mapped to a
        sequence ID) and each column is a basepair. Missing values for 
        varying sequence lengths can be represented by a nan
    start, end: Series
        A series of length `n` (where `n` corresponds to the number of rows
        [sequences] in `seqs`) and the value represents a nucleotide in that
        sequence. EVery sequence must have a defined starting and ending 
        position.

    Returns
    -------
    DataFrame
        An n x k array of sequences where each row is a sequence (mapped to a
        sequence ID) and each column is a basepair. The represented sequence
        will be the sequence between `start_` and `end_`
    """
   # Looks only at places where the position has been defined.

    # Gets the positions
    pos_ref = (np.ones((seqs.shape[0], 1)) * 
        np.atleast_2d(np.arange(0, seqs.shape[1])))
    start_mask = \
        (np.ones((seqs.shape[1], 1)) * np.atleast_2d(start.values)).T
    end_mask = \
        (np.ones((seqs.shape[1], 1)) * np.atleast_2d(end.values)).T
    
    # Generates a mask of places where the sequence is not contained
    pos_mask = ((pos_ref < start_mask) | (pos_ref > end_mask) 
                & ~(np.isnan(start_mask) | np.isnan(end_mask)))
    
    # Masks the sequence
    wrk_seq = seqs.mask((pos_ref < start_mask) | (pos_ref > end_mask), '')
    wrk_seq.dropna(inplace=True)

    # String comprenhsion hack to get what we need... 
    sub_seq = wrk_seq.apply(lambda x: pd.Series(list(''.join(x))), axis=1)

    return sub_seq


def _build_id_map(kmers, region):
    """
    A silly function to map a collapsed kmer back to a summary

    Parameters
    ----------
    kmers : array-like
        A listing of the kmers to be seperated
    region: str
        The hypervariable region to profile.

    Returns
    -------
    DataFrame
        A mapping betweent the `db_seq`, its associated kmers, and the
        individual sequence name, and the region.
    """
    expansion = pd.DataFrame.from_dict(orient='index', data={
        kmer: pd.Series(kmer.split("|")) for kmer in np.unique(kmers)
        })
    expansion.index.set_names('kmer', inplace=True)
    expansion.reset_index(inplace=True)
    long_ = expansion.melt(id_vars=['kmer'], value_name='seq_name').dropna()
    long_.drop(columns=['variable'], inplace=True)
    long_['db_seq'] = long_['seq_name'].apply(lambda x: x.split("@")[0])
    long_['region'] = region
    long_.sort_values(['db_seq', 'seq_name'], inplace=True)
    long_.drop_duplicates(['db_seq', 'seq_name'], inplace=True)
    return long_.set_index('db_seq')[['seq_name', 'kmer', 'region']]

