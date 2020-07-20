import itertools as it
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)

import dask
import numpy as np
import pandas as pd
import sparse as sp

from q2_types.feature_data import (DNAFASTAFormat, DNAIterator)

from q2_sidle._utils import (_setup_dask_client, 
                             _convert_generator_to_seq_block,
                             _convert_generator_to_delayed_seq_block, 
                             _convert_seq_block_to_dna_fasta_format,
                             degen_reps,
                             )


def align_regional_kmers(kmers: DNAFASTAFormat, 
    rep_seq: DNAFASTAFormat, 
    region: str, 
    max_mismatch: int=2, 
    chunk_size:int=1000, 
    debug:bool=False, 
    n_workers:int=0,
    client_address:str=None) -> (pd.DataFrame, DNAFASTAFormat):
    """
    Performs regional alignment between database "kmers" and ASVs

    Parameters
    ----------
    kmers : DNAFastaFormat
        The set of reference sequences extracted from the database. These are
        assumes to be start in the same position of the 16s rRNA sequence as 
        the sequence being tested and assumed to be the same length as the
        ASVs being aligned.
    rep_seq: DNAFastaFormat
        The representative sequences for the regional ASV table being aligned.
        These are assumed to start at the same position as the kmers and 
        should be trimmed to the same length.
    region: str
        An identifier for the region. Ideally, this matches the identifier 
        used in the reference region map
    max_mismatch: int
        the maximum number of mismatched nucleotides allowed in mapping 
        between a sequence and kmer.
    debug: bool
        Whether the function should be run in debug mode (without a client)
        or not. `debug` superceeds all options
    n_workers: int, optional
        The number of jobs to initiate. When `n_workers` is 0, the cluster 
        will be able to access all avaliable resources.

    Returns
    -------
    DataFrame
        A mapping between the kmer (`kmer`) and the asv (`asv`), including 
        the region (`region`), number of mismatched basepairs (`mismatch`) and 
        the sequence length (`length`).
    DNAFASTAFormat
        The ASVs which could not be aligned to kmers

    """
     # Sets up the client
    _setup_dask_client(debug=debug, cluster_config=None,  
                       n_workers=n_workers, address=client_address)
    # Gets the sequences
    kmers = _convert_generator_to_delayed_seq_block(
        kmers.view(DNAIterator), chunk_size
        )
    rep_set = _convert_generator_to_seq_block(
        rep_seq.view(DNAIterator), chunk_size
        )

    # Performs the alignment
    aligned = np.hstack([
        dask.delayed(_align_kmers)(kmer, asv, max_mismatch)
        for kmer, asv in it.product(kmers, rep_set)
        ])
    aligned = pd.concat(axis=0, sort=False, objs=dask.compute(*aligned))

    aligned['region'] = region
    aligned.sort_values(['kmer', 'asv'], inplace=True)
    aligned.set_index('kmer', inplace=True)

    rep_set = pd.concat(rep_set, axis=0, sort=False)
    aligned_asvs = aligned['asv'].unique()
    
    discard = _convert_seq_block_to_dna_fasta_format(
        [rep_set.copy().drop(aligned_asvs)]
    )

    return aligned, discard


def _align_kmers(reads1, reads2, allowed_mismatch=2, read1_label='kmer', 
    read2_label='asv', allow_degen1=False, allow_degen2=False, 
    expand_match=False, block_size=500):
    """
    Performs a kmer-based alignment between two groups of n-mers

    Parameters
    ----------
    reads1, reads2 : Dataframe
        The sequences to be aligned where the sequence identifer is given in
        the index and each column contains a nucleotide.
    mismatch : int, optional
        The number of mismatches allowed between the two sets of sequences.
    read1_label, read2_label: str, optional
        A way to refer to the sequences in each alignment set
    allow_degen1, allow_degen2 : bool, optional
        Whether degeneracy should be allowed in the alignment. If 
        `allow_degen1` and `allow_degen2`, then sequences
    expand_match : bool, optional
        If `allow_degen1 = True` and `allow_degen2 = True` then whether we
        shoudl allow watches for degenerate positions that map to the same 
        basepair. For example, when `expand_match == True` the ambigious 
        nucleotide, "W" (A and T) could be matched to "N" (A, T, C, G) but 
        not "S" (G, C).

    Returns
    -------
    pd.DataFrame
        A long-form dataframe giving the two read identifiers and the number
        of nt that do not match.

    Raises
    ------
    ValueError
        If the two sets fo reads are of different lengths. (We can't align 
        kmers here that are different lengths).

    """
    # Gets the dimensions of the reads
    num_reads1, kmer_len = (reads1.shape)
    num_reads2, length2 = (reads2.shape)

    # Raises an error if the 
    if kmer_len != length2:
        raise ValueError('The reads are different lengths. Alignment cannot'
            ' be preformed.')

    # Determines the mismatch threshhold
    miss_thresh = (kmer_len - allowed_mismatch)

    # # gets the sequence identifiers
    read1_ids = reads1.index
    read2_ids = reads2.index


    # A couple of helper functions that just increase readability and throw things
    # to sparse because that will hopefully decreaes the memory requirements
    # because things hould be sparse 
    def _expand1(x):
        return sp.as_coo(np.dstack([x] * num_reads2).swapaxes(1, 2))
    def _expand2(x):
        return sp.as_coo(
            np.dstack([x] * num_reads1).swapaxes(0, 2).swapaxes(2, 1))

    if (allow_degen1 | allow_degen2):

        def check_nt(reads, nt, allow):
            return ((reads.isin([nt])) | 
                    (reads.isin(degen_reps[nt]) & allow)).values

        # Expands read1 to get positions and degeneracies
        reads1_d = (_expand1(~reads1.isin(['A', 'C', 'T', 'G']).values) & 
                             allow_degen1)
        reads1_A = _expand1(check_nt(reads1, 'A', allow_degen1))
        reads1_C = _expand1(check_nt(reads1, 'C', allow_degen1))
        reads1_G = _expand1(check_nt(reads1, 'G', allow_degen1))
        reads1_T = _expand1(check_nt(reads1, 'T', allow_degen1))

        # Expands read2 to get positions and degeneracies
        reads2_d = (_expand2(~reads2.isin(['A', 'C', 'T', 'G']).values) & 
                             allow_degen2)
        reads2_A = _expand2(check_nt(reads2, 'A', allow_degen2))
        reads2_C = _expand2(check_nt(reads2, 'C', allow_degen2))
        reads2_G = _expand2(check_nt(reads2, 'G', allow_degen2))
        reads2_T = _expand2(check_nt(reads2, 'T', allow_degen2))

        readsA = reads1_A & reads2_A
        readsC = reads1_C & reads2_C
        readsG = reads1_G & reads2_G
        readsT = reads1_T & reads2_T

        # Determines whether we allow mismatch between two degeneracies.
        # We want to keep any position where there are not two degeneracies OR
        # positions where there are two degeneries IF we flagged that with 
        # `expand_match`. Basically, its a stupid complex way to multiply by 1.
        # But, I think it will work here.
        extend_mat = ((reads1_d & reads2_d & expand_match) | 
                      ~(reads1_d & reads2_d))

        # Builds matching array and determines the total matches in each sequence
        # then filters down ot only positions that match above the threshhold
        pos_match = ((readsA | readsC | readsG | readsT) & 
                     extend_mat).sum(axis=-1)
    else:
        reads1_exp = np.dstack([reads1.values] * num_reads2).swapaxes(1, 2)
        reads2_exp = np.dstack([reads2.values] * 
                               num_reads1).swapaxes(0, 2).swapaxes(2, 1)

        pos_match = sp.as_coo((reads1_exp == reads2_exp).sum(axis=-1))


    pos_match = pos_match * (pos_match >= miss_thresh)

    match_values = pd.DataFrame(
        data=(kmer_len - pos_match.todense()), 
        index=pd.Index(read1_ids, name=read1_label),
        columns=pd.Index(read2_ids, name=read2_label)
        ).reset_index()

    match_long = match_values.melt(id_vars=read1_label, 
                                   var_name=read2_label,
                                   value_name='mismatch' 
                                   )
    match_long = match_long[match_long['mismatch'] <= allowed_mismatch].copy()
    match_long['length'] = kmer_len
    match_long.sort_values([read1_label, read2_label], inplace=True)

    return match_long[[read1_label, read2_label, 'mismatch', 'length']]


def _to_dask_array_3d(df, nstack):
    """
    A function to convert from a dask dataframe to an array.
    """

    return da.dstack([(df.to_dask_array()).compute_chunk_sizes()] * nstack)