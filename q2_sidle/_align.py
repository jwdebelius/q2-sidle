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
    chunk_size:int=1000, 
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
    
    # Gets the sequences
    kmers = kmers.astype(str)

    # Performs the alignment
    aligned = _align_kmers(kmers, rep_seq.astype(str), 
                           allowed_mismatch=max_mismatch,
                           block_size=chunk_size)
    aligned['region'] = region
    aligned = aligned.compute()
    
    # Gets the discarded sequences
    discard = aligned.groupby('asv')['discard'].all()
    discard = discard.index[discard].values
    
    # Keeps only the aligned sequences
    aligned = aligned.loc[~aligned['discard']]
    aligned.drop(columns=['discard'], inplace=True)
    
    aligned.sort_values(['kmer', 'asv'], inplace=True)
    
    return aligned, rep_seq.loc[discard]


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
    
    # Checks the read lengths of both reads
    num_reads1, length1 = _check_read_lengths(reads1, read1_label)
    num_reads2, length2 = _check_read_lengths(reads2, read2_label)
    
    # The read lengths must also be the same length
    if length1 != length2:
        raise ValueError('The reads are different lengths. Alignment cannot'
                         ' be preformed.')

    # Gets the threshhoold
    miss_thresh = (length1 - allowed_mismatch)
    
    # Handles degenerate sequences if they should be replaced.
    if allow_degen1:
        reads1 = reads1.apply(lambda x: pd.Series(list(x)))
        reads1.replace(degen_sub2, inplace=True)
        reads1 = reads1.apply(lambda x: ''.join(x), axis=1)

    if allow_degen2:
        reads2 = reads2.apply(lambda x: pd.Series(list(x)))
        reads2.replace(degen_sub2, inplace=True)
        reads2 = reads2.apply(lambda x: ''.join(x), axis=1)
    
    # Puts together a long dask delayed array for matching the sequences by
    # first, converting to a delayed object and then introducing the secondary
    # sequences as columns, and then melting them
    read_joint = dd.from_pandas(reads1.reset_index(), chunksize=block_size)
    read_joint.columns = [read1_label, 'sequence1']
    read_joint = dd.concat(axis=1, dfs=[
        read_joint,
        read_joint[read1_label].apply(lambda x: reads2, 
                                      meta={i: 'object' for i in reads2.index})
    ])
    read_joint = read_joint.melt(id_vars=[read1_label, 'sequence1'],
                                 var_name=read2_label,
                                 value_name='sequence2'
                                 ) 
    read_joint['length'] = length1

    # Computes the number of sequences that are different between the two 
    # sequences
    def _compute_joint(x):
        seq1 = '(%s){e<=%i}' % (x['sequence1'], x['length'])
        seq2 = x['sequence2']
        m = regex.match(seq1, seq2)
        return m.fuzzy_counts[0]
    read_joint['mismatch'] = \
        read_joint.apply(_compute_joint, axis=1, meta=(None, 'int64'))
    
    # Drops out the sequences 
    read_joint = read_joint.drop(columns=['sequence1', 'sequence2'])
    
    read_joint['discard'] = read_joint['mismatch'] > allowed_mismatch

    return read_joint

    
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
    
