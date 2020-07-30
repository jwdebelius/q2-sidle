import os

import dask
import numpy as np
import pandas as pd
import regex
import skbio

from dask.distributed import Client

from qiime2 import Artifact, Metadata
from qiime2.plugin import ValidationError
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
degen_sub2 = {k: '[%s]' % v for k, v in degen_sub.items() if (k in degen)}
degen_undo = {'AG': 'R',
              'CT': 'Y',
              'CG': 'S',
              'AT': 'W',
              'GT': 'K',
              'AC': 'M',
              'CGT': 'B',
              'AGT': 'D',
              'ACT': 'H',
              'ACG': 'V',
              'ACGT': 'N'}

def _setup_dask_client(debug=False, cluster_config=None, n_workers=1,
    address=None):
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
    address: str, optional
        The IP address for the client
    """

    if debug:
        pass
    elif cluster_config is not None:
        client = Client(**client_config.to_dict())
        print(client.dashboard_link)
    elif address is not None:
        client = Client(address)
        print(client.dashboard_link)
    elif n_workers == 0:
        client = Client()
        print(client.dashboard_link)
    else:
        client = Client(n_workers=n_workers, processes=True)
        print(client.dashboard_link)



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


def _check_manifest(manifest):
    """
    Makes sure that everything is in the manifest

    Parameters
    ---------
    Manifest : qiime2.Metadata
        A manifest file describing the relationship between regions and their
        alignment mapping. The manifest must have at least three columns
        (`kmer-map`, `alignment-map` and `frequency-table`) each of which
        contains a unique filepath. 
    """
    manifest = manifest.to_dataframe()
    cols = manifest.columns
    if not (('kmer-map' in cols) & ('alignment-map' in cols) & 
            ('frequency-table' in cols) & ('region-order' in cols)):
        raise ValidationError('The manifest must contain the columns '
                              'kmer-map, alignment-map and frequency-table.\n'
                              'Please check the manifest and make sure all'
                              ' column names are spelled correctly')
    manifest = manifest[['kmer-map', 'alignment-map', 'frequency-table']]
    if pd.isnull(manifest).any().any():
        raise ValidationError('All regions must have a kmer-map, '
                              'alignment-map and frequency-table. Please '
                              'check and make sure that you have provided '
                              'all the files you need')
    if (len(manifest.values.flatten()) != 
            len(np.unique(manifest.values.flatten()))):
        raise ValidationError('All paths in the manifest must be unique.'
                             ' Please check your filepaths')
    if not np.all([os.path.exists(fp_) for fp_ in manifest.values.flatten()]):
        raise ValidationError('All the paths in the manifest must exist.'
                             ' Please check your filepaths')
    if not np.all([os.path.isfile(fp_) for fp_ in manifest.values.flatten()]):
        raise ValidationError('All the paths in the manifest must be files.'
                             ' Please check your filepaths')


def _read_manifest_files(manifest, dataset, semantic_type=None, view=None):
    """
    Extracts files from the manifest and turns them into a list of objects
    for analysis
    """
    paths = manifest.get_column(dataset).to_series()
    artifacts = [Artifact.load(path) for path in paths]
    if semantic_type is not None:
        type_check = np.array([str(a.type) == semantic_type 
                               for a in artifacts])
        if not np.all(type_check):
            err_ = '\n'.join([
                'Not all %s Artifacts are of the %s semantic type.' 
                    % (dataset.replace('-', ' '),   semantic_type),
                'Please review semantic types for these regions:',
                '\n'.join(paths.index[type_check == False])
                ])
            raise TypeError(err_)
    if view is not None:
        return [a.view(view) for a in artifacts]
    else:
        return artifacts



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