import itertools as it
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
import skbio
from skbio import DNA

from qiime2 import Metadata
from q2_types.feature_data import (DNAFASTAFormat, DNAIterator)
from q2_sidle._utils import (_setup_dask_client, 
                             )
from q2_feature_classifier._skl import _chunks


complement = {'A': 'T', 'T': 'A', 'G': 'C',  'C': 'G',
              'R': 'Y', 'Y': 'R', 'S': 'S',  'W': 'W',
              'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
              'D': 'H', 'H': 'D', 'N': 'N'}


def prepare_extracted_region(sequences: DNAFASTAFormat, 
    region:str, 
    trim_length:int, 
    fwd_primer:str, 
    rev_primer:str, 
    reverse_complement_rev:bool=True,
    reverse_complement_result:bool=False,
    chunk_size:int=10000, 
    debug:bool=False, 
    n_workers:int=0,
    client_address:str=None,
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
        A mapping between the kmer sequence name and the the full database 
        sequence name, along with regional information
    """

    # Sets up the client
    _setup_dask_client(debug=debug, cluster_config=None,  
                       n_workers=n_workers, address=client_address)

    # Reverse complements the reverse primer
    if reverse_complement_rev:
        rev_primer = str(DNA(rev_primer).reverse_complement())

    # Reads in the sequences
    sequences = sequences.view(DNAIterator)
    seq_blocks = [_block_seqs(seq)
                  for seq in _chunks(sequences, int((chunk_size)))]
    # Makes the fake extraction position based on the trim length
    fragment = [_artifical_trim(seq, trim_length) for seq in seq_blocks]
    # Prepares the amplicon for collapsing
    condensed = dd.from_delayed([
        _condense_seqs(seq) for seq in fragment],
        meta=[('amplicon', 'str'), ('seq-name', 'str')]
    )
    # Writes the 
    ff, group2 = _collapse_all_sequences(condensed, reverse_complement_result)
    ids = _expand_ids(group2, fwd_primer, rev_primer, region, trim_length,
                      chunk_size)

    return (ff, ids.compute().set_index('db-seq').sort_index())


@dask.delayed
def _artifical_trim(seqs, trim_length):
    seqs['extract-length'] = seqs['sequence'].apply(lambda x: len(x))
    keep = (np.absolute(seqs['extract-length']) >= np.absolute(trim_length))
    seqs = seqs.loc[keep].copy()

    if trim_length > 0:
        seqs['amplicon'] = seqs['sequence'].apply(lambda x: x[:trim_length])
    else:
        seqs['amplicon'] = seqs['sequence'].apply(lambda x: x[trim_length:])

    seqs['db-seq'] = seqs['seq-name'].apply(lambda x: x.split('@')[0])
    seqs.drop_duplicates(['db-seq', 'amplicon'], inplace=True)

    return seqs[['seq-name',  'amplicon']]


@dask.delayed
def _block_seqs(seqs, degen_thresh=3):
    """
    Converts the sequences into an expanded sequence block
    """
    s2 = pd.concat([_expand_degenerate_gen(seq, degen_thresh=degen_thresh) 
                    for seq in seqs]).astype(str)
    s2.index.set_names('seq-name', inplace=True)
    s2.name = 'sequence'
    s3 = s2.reset_index()
    s3['db-seq'] = s3['seq-name'].apply(lambda x: x.split("@")[0])
    return s3


def _collapse_all_sequences(condensed, reverse_complement_result):
    """
    Collapse and reverse complelent the results and tidies it
    """
    # Collapses the data into grouped data and prints the computed result
    grouped = condensed.groupby('amplicon')['seq-name'].apply(
        lambda x: '>%s' % "|".join(np.sort(x.values)), meta=('seq-name', str))  
    group2 = grouped.compute().reset_index()
    group2.sort_values('seq-name', inplace=True)
    # Reverse complemnents the sequences if desired
    if reverse_complement_result:
        def _rc(x):
            seq = skbio.DNA(x['amplicon'])
            return seq.reverse_complement()
        group2['seq'] = group2.apply(_rc, axis=1).astype(str)
    else:
        group2['seq'] = group2['amplicon']
    
    # Saves the sequences
    ff = DNAFASTAFormat()
    group2.to_csv(str(ff), header=None, sep='\n', index=False,
                  columns=['seq-name', 'seq'])
    return ff, group2


@dask.delayed
def _condense_seqs(seqs):
    """
    Collapses duplicate sequences before processing
    """
    seqs.sort_values(['amplicon', 'seq-name'], inplace=True)
    fragment = \
        seqs.groupby(['amplicon'])['seq-name'].apply(lambda x: "|".join(x))
    return fragment.sort_index().reset_index()


def _expand_degenerate_gen(seq_, degen_thresh=3):
    """
    Expands the degenerate sequences in the seq blocks
    """
    id_ = seq_.metadata['id']
    if seq_.has_degenerates():
        expand = pd.Series([s for s in seq_.expand_degenerates()]).astype(str)
        expand.sort_values(inplace=True)
        expand.reset_index(inplace=True, drop=True)
        expand.rename(index={i: '%s@%s' % (id_, str(i + 1).zfill(4)) 
                             for i in expand.index},
                      inplace=True)
    else:
        expand = pd.Series({id_: seq_}).astype(str)
    return expand


def _expand_ids(group2, fwd_primer, rev_primer, region, trim_length, 
    chunk_size):
    """
    Expands the seq-name column from grouped IDs into an ID map to be used
    in alignment

    Parameters
    ----------
    group2: DataFrame
        A dataframe which maps the trimmed sequence identifer (`seq-name`)
        to the sequence (`amplicon`)
    fwd_primer, rev_primer: str
        The forward and reverse primers to be amplified
    trim_length : int, optional
        The length sequences should be trimmed to match the asvs for '
        'kmer-based alignment.
    region : str
        The hypervariable region to profile.
    chunk_size: int, optional
        The number of sequences to group for analysis

    Returns
    -------
    DataFrame
        A mapping between the kmer sequence name and the the full database 
        sequence name, along with regional information
    """

    # Collapses the data into grouped data

    
    # Sorts out the ids
    ids = dd.from_pandas(group2['seq-name'].apply(lambda x: x.strip('>')),
                                                  chunksize=chunk_size)
    ids.name = 'kmer'
    ids = dd.from_delayed([_split_ids(seq_) for seq_ in ids.to_delayed()],
                          meta=[('kmer', str), ('seq-name', str)]
                          )
    ids['db-seq'] = ids['seq-name'].apply(lambda x: x.split("@")[0], 
                                          meta=('db-seq', str))
    ids['fwd-primer'] = fwd_primer
    ids['rev-primer'] = rev_primer
    ids['region'] = region
    ids['kmer-length'] = trim_length

    return ids


@dask.delayed
def _split_ids(ids):
    """
    Splits collapsed kmers into single sequences for a database map
    """
    kmers = pd.concat(axis=1, objs=[
        ids, ids.apply(lambda x: pd.Series(x.split("|")))
    ])
    kmers = kmers.melt(id_vars='kmer', value_name='seq-name')
    kmers.drop(columns=['variable'], inplace=True)
    return kmers.dropna()
