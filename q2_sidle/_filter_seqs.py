import dask
import dask.dataframe as dd
import numpy as np
import regex
import pandas as pd
from skbio import DNA

from q2_sidle._utils import (_setup_dask_client, 
                             _convert_generator_to_delayed_seq_block,
                             _convert_generator_to_seq_block, 
                             _convert_seq_block_to_dna_fasta_format,
                             _count_degenerates,
                             database_params,
                             )
from q2_types.feature_data import (DNAFASTAFormat,
                                   DNAIterator, 
                                   DNASequencesDirectoryFormat
                                   )


def filter_degenerate_sequences(sequences: DNASequencesDirectoryFormat, 
    max_degen:int=3,
    chunk_size:int=10000,
    debug:bool=False, 
    n_workers:int=0,
    client_address: str=None,
    ) -> pd.Series:
    """
    Prefilters the database to remove sequences with too many degenerates

    The original SMURF paper recommends filtering sequences with more than 
    3 degenerate nucleotides along their length to limit alignment issues and
    memory.

    Parameters
    ----------
    sequences: q2_types.DNASequencesDirectoryFormat
        The directory of reference sequences
    max_degen : int, optional
        The maximum number of degenerate nucleotides necessary to retain the
        sequence in the database
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
        The fitlered reads
    # """
    # # Sets up the client
    _setup_dask_client(debug=debug, cluster_config=None,  
                       n_workers=n_workers, address=client_address)
    
    # Reads in the sequences
    sequences = dd.from_pandas(
        sequences.file.view(pd.Series).reset_index(), 
        chunksize=chunk_size)
    sequences.columns = ['id', 'sequence']
    sequences['sequence'] = sequences['sequence'].astype(str)
    def _count_degen(x):
        return len(regex.findall('[RYSWKMBDHVN]', x)) <= max_degen

    sequences['count'] = sequences['sequence'].apply(_count_degen, 
                                                     meta=('sequence', bool))
    sequences = sequences.loc[sequences['count']]

    sequences['skbio'] = \
        sequences.apply(
            lambda x: DNA(x['sequence'], metadata={'id': x['id']}), 
            axis=1,
            meta=(None, 'object')
            )
    sequences = sequences.compute()
    sequences = sequences.set_index('id')

    return sequences['skbio']


