
import dask
import pandas as pd
import skbio

from dask.distributed import Client

from q2_feature_classifier._skl import _chunks
from q2_types.feature_data import DNAIterator, DNAFASTAFormat



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
        client = Client(n_workers=threads, processes=True)


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
    @dask.delayed
    def _to_seq_array(x):
        return pd.DataFrame({s.metadata['id']:  list(str(s)) for s in x}).T

    seq_block = [_to_seq_array(seqs) 
                 for seqs in _chunks(generator, chunksize)]

    return seq_block


def _convert_generator_to_seq_block(generator, chunksize=5000):
    """
    Converts from a generator to a block of sequences
    """
    def _to_seq_array(x):
        return  pd.DataFrame({s.metadata['id']:  list(str(s)) for s in x}).T

    seq_block = [_to_seq_array(seqs) for seqs in _chunks(generator, chunksize)]

    return seq_block

