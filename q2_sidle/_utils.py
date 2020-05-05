
import dask
import numpy as np
import pandas as pd
import skbio

from dask.distributed import Client

from q2_feature_classifier._skl import _chunks
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

### Nucleotides
degenerate_map = {"R": ['A', 'G'],
                  'Y': ['C', 'T'],
                  'S': ['G', 'C'],
                  'W': ['A', 'T'],
                  'K': ['G', 'T'],
                  'M': ['A', 'C'],
                  'B': ['C', 'G', 'T'],
                  'D': ['A', 'G', 'T'],
                  'H': ['A', 'C', 'T'],
                  'V': ['A', 'C', 'G'],
                  'N': ['A', 'C', 'G', 'T'],
                  } 

defined = ['A', 'C', 'T', 'G']
degen = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
degen_reps = {'A': ['R', 'W', 'M', 'D', 'H', 'V', 'N'],
              'G': ['R', 'S', 'K', 'B', 'D', 'V', 'N'],
              'T': ['Y', 'W', 'K', 'B', 'D', 'H', 'N'],
              'C': ['Y', 'S', 'M', 'B', 'H', 'V', 'N'],
              }


def _setup_dask_client(debug=False, cluster_config=None, n_workers=1):
    """
    Sets up a Dask client and daskboard

    Parameters
    ----------
    debug: bool
        Whether the function should be run in debug mode (without a client)
        or not. `debug` superceeds all options
    client_config: dict, optional
        A dictionary describing configuration parameters for the dask client.
        More information about configuring the dask scheduler and dask client 
        can be found at
            https://docs.dask.org/en/latest/setup/single-distributed.html
        The client_config sueprceeds the n_workers value, so if you want 
        multi threading, that should be specified here.
    n_workers: int, optional
        The number of jobs to initiate. When `n_workers` is 0, the cluster 
        will be able to access all avalaibel resources.
    """

    if debug:
        pass
    elif cluster_config is not None:
        client = Client(**client_config.to_dict())
    elif n_workers == 0:
        client = Client()
    else:
        client = Client(n_workers=n_workers, processes=True)


def _convert_seq_block_to_dna_fasta_format(seqs):
    """
    Converts to a DNA fasta format
    """
    seqs = pd.concat(axis=0, sort=True, objs=seqs).fillna('')
    seqs = seqs.apply(lambda x: ''.join(x), axis=1)
    ff = DNAFASTAFormat()
    with ff.open() as f:
        for id_, seq_ in seqs.items():
            sequence = skbio.DNA(seq_, metadata={'id': id_})
            skbio.io.write(sequence, format='fasta', into=f)
    return ff


def _convert_generator_to_delayed_seq_block(generator, chunksize=5000):
    """
    Converts from a generator to a block of sequences
    """
    seq_block = [dask.delayed(_to_seq_array)(seqs) 
                 for seqs in _chunks(generator, chunksize)]

    return seq_block


def _convert_generator_to_seq_block(generator, chunksize=5000):
    """
    Converts from a generator to a block of sequences
    """
    seq_block = [_to_seq_array(seqs) for seqs in _chunks(generator, chunksize)]
    return seq_block


def _to_seq_array(x):
    """
    Converts a list of sequences from a generator to a DataFrame of sequences
    """
    seq_series = pd.Series({s.metadata['id']: str(s) for s in x})
    return seq_series.apply(lambda x: pd.Series(list(x)))



def _count_degenerates(seq_array):
    """
    Determines the number of degenerate nt in a sequence

    Parameters
    ----------
    seq_array: Dataframe
        A dataframe where each row is a sequence and each column is a basepair
        in that sequence

    Returns
    -------
    Series
        The number of degenerate basepairs in each array

    """
    any_degen = (~seq_array.isin(['A', 'G', 'T', 'C', np.nan]))
    num_degen = any_degen.sum(axis=1)
    
    return num_degen

    

database_params = {
    'greengenes': {
        'delim': '; ',
        # 'levels': ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'],
        'defined': lambda x: len(x) > 3,
        'inherient': lambda x: 'unsp. %s' % x.replace('__', '. '),
        'contested': lambda x: x.replace('[', 'cont. ').replace(']', '')
    },
    'silva': {
        'delim': ';',
        # 'levels': [],
        'defined': lambda x: ~(('uncul' in  x) | ('metagenome' in x)),
        'inherient': lambda x: x,
        'contested': lambda x: x.replace('[', 'cont. ').replace(']', ''),
    },
    'homd': {
        'delim': ';',
        # 'levels': [],
        'defined': lambda x: True,
        'inherient': lambda x: x,
        'contested': lambda x: x.replace('[', 'cont. ').replace(']', ''),
    },
    'none': {
        'delim': ';',
        # 'levels': [],
        'defined': lambda x: True,
        'inherient': lambda x: x,
        'contested': lambda x: x,
    },
    }