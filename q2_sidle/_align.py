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
    chunk_size:int=250, 
    debug:bool=False, 
    n_workers:int=0,
    client_address:str=None) -> (pd.DataFrame, pd.Series):
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

    # Checks the read lengths 
    _check_read_lengths(kmers, 'kmers')
    _check_read_lengths(rep_seq, 'rep_seq')

    if len(kmers.iloc[0]) != len(rep_seq.iloc[0]):
        raise ValueError('The  kmer sequence length must match the '
                         'rep-seq sequence length')
    
    # Gets the sequences
    kmers = kmers.astype(str)
    kmers = [kmers[i:(i + chunk_size)] 
             for i in np.arange(0, len(kmers), chunk_size)]
    rep_seq = rep_seq.astype(str)
    rep_seq = [rep_seq[i:(i + chunk_size)] 
               for i in np.arange(0, len(rep_seq), chunk_size)]

    # Performs the alignment
    aligned = np.hstack([
        dask.delayed(_align_kmers)(kmer, asv, max_mismatch)
        for kmer, asv in it.product(kmers, rep_seq)
        ])
    aligned = pd.concat(dask.compute(*aligned))
    aligned['region'] = region

    # Gets the discarded sequences
    discard = aligned.groupby('asv')['discard'].all()
    discard = discard.index[discard].values
    
    # Keeps only the aligned sequences
    aligned = aligned.loc[~aligned['discard']]
    aligned.drop(columns=['discard'], inplace=True)
    aligned = aligned[['kmer', 'asv', 'length', 'mismatch', 'region']]
    
    aligned.sort_values(['kmer', 'asv'], inplace=True)
    
    return aligned, pd.concat(rep_seq).loc[discard]


def _align_kmers(reads1, reads2, allowed_mismatch=2, read1_label='kmer', 
    read2_label='asv'):
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
    reads1 = reads1.apply(lambda x: pd.Series(list(x)))
    reads2 = reads2.apply(lambda x: pd.Series(list(x)))
    num_reads1, kmer_len = (reads1.shape)
    num_reads2, length2 = (reads2.shape)

    # gets the sequence identifiers
    read1_ids = reads1.index
    read2_ids = reads2.index

    reads1_exp = np.dstack([reads1.values] * num_reads2).swapaxes(1, 2)
    reads2_exp = np.dstack([reads2.values] * 
                            num_reads1).swapaxes(0, 2).swapaxes(2, 1)

    pos_match = ((reads1_exp == reads2_exp).sum(axis=-1))

    match_values = pd.DataFrame(
        data=(kmer_len - pos_match), 
        index=pd.Index(read1_ids, name=read1_label),
        columns=pd.Index(read2_ids, name=read2_label)
        ).reset_index()

    match_long = match_values.melt(id_vars=read1_label, 
                                   var_name=read2_label,
                                   value_name='mismatch' 
                                   )
    match_long['discard'] = match_long['mismatch'] > allowed_mismatch
    match_long['length'] = kmer_len
    match_long.sort_values([read1_label, read2_label], 
                           inplace=True)

    return match_long[[read1_label, read2_label, 
                      'length', 'mismatch', 'discard']]


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
    
