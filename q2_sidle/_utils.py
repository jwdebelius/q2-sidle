import os

import dask
import numpy as np
import pandas as pd
import regex
import skbio

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



def _check_regions(region):
    """
    Converts the region to a numeric assignment
    """
    region, region_idx = np.unique(region, return_index=True)
    region_order = {region: i for (i, region) in zip(*(region_idx, region))}
    region_names = {i: r for r, i in region_order.items()}
    num_regions = len(region_order)
    
    return region_order, region_names, num_regions


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