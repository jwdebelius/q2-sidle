import itertools as it
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
import regex

from q2_types.feature_data import (DNAFASTAFormat, DNAIterator)

from q2_sidle._utils import (_setup_dask_client, 
                             _convert_generator_to_seq_block,
                             _convert_generator_to_delayed_seq_block, 
                             _convert_seq_block_to_dna_fasta_format,
                             degen_sub2
                             )


def align_regional_kmers(kmers: pd.Series, 
    rep_seq: pd.Series, 
    region: str, 
    max_mismatch: int=2, 
    chunk_size:int=100, 
    debug:bool=False, 
    n_workers:int=0,
    client_address:str=None) -> (dd.DataFrame, pd.Series):
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

    # Checks the kmer size
    num_kmers, kmer_length = _check_read_lengths(kmers, 'kmer')
    kmers = dd.from_pandas(kmers.astype(str), 
                           chunksize=chunk_size)

    # Converts the representative sequences to a delayed object
    num_asvs, asv_length = _check_read_lengths(rep_seq, 'rep_seq')
    rep_seq_ids  = rep_seq.index.values
    rep_seq = dd.from_pandas(rep_seq.astype(str),
                             chunksize=chunk_size)

    if kmer_length != asv_length:
        raise ValueError('The kmer and ASV sequences must be the same length')


    # Performs the alignment
    aligned = np.hstack([
        dask.delayed(_align_kmers)(kmer, asv, max_mismatch)
        for kmer, asv in it.product(kmers.to_delayed(), rep_seq.to_delayed())
        ])
    aligned = dd.from_delayed(aligned)
    aligned_asvs = rep_seq_ids[~np.isin(rep_seq_ids,
                                       aligned['asv'].unique().compute())]
    print(aligned_asvs)

    aligned['region'] = region


    discard = rep_seq.loc[aligned_asvs]

    return aligned, discard.compute()


def _align_kmers(reads1, reads2, allowed_mismatch=2, read1_label='kmer', 
    read2_label='asv', block_size=500):
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
    reads1 = reads1.apply(lambda x: np.array(list(x)))
    reads2 = reads2.apply(lambda x: np.array(list(x)))
    length = len(reads1.iloc[0])
    match = pd.Series({
        (id1, id2) : np.sum(seq1 == seq2)
        for (id1, seq1), (id2, seq2) 
        in it.product(reads1.items(), reads2.items())
        })
    match = length - match
    match = match.loc[match <= allowed_mismatch]
    match.index.set_names([read1_label, read2_label], inplace=True)
    match.name = 'mismatch'
    match = match.reset_index()
    match['length'] = length

    return match[[read1_label, read2_label, 'length', 'mismatch']]


def _check_read_lengths(reads, read_label):
    """
    Checks the length of the sequences
    """
    read_lengths = reads.apply(lambda x: len(x))
    length_dist = read_lengths.value_counts()
    if len(length_dist) > 1:
        raise ValueError('The %s reads are not a consistent length' 
                          % read_label)

    return (length_dist.values[0], length_dist.index[0])
    
