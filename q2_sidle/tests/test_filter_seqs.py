from unittest import TestCase, main

import os

import dask
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact
from q2_sidle._filter_seqs import (filter_degenerate_sequences,
                                   _count_degenerates
                                   )

class FilterTest(TestCase):
    def setUp(self):
        self.seqs = pd.DataFrame(
            data=[list('AGNNNNNNN'), 
                  list('CATCATCAT'), 
                  list('CATCATCAT')],
            index=['0', '1', '2'],
            )
        self.dna_fasta = Artifact.import_data(
            'FeatureData[Sequence]', 
            pd.Series({'0': skbio.DNA('AGNNNNNNN', metadata={'id': '0'}),
                       '1': skbio.DNA('CATCATCAT', metadata={'id': '1'}),
                       '2': skbio.DNA('CATCATCAT', metadata={'id': '2'}),
                       }))
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')

    def test_count_degenerates(self):
        seq_array = pd.DataFrame.from_dict(orient='index', data={
            0: {'id_': '0', 0: 'A', 1: 'G', 2: 'T', 3: 'C'},
            1: {'id_': '1', 0: 'A', 1: 'R', 2: 'W', 3: 'S'},
            2: {'id_': '2', 0: 'C', 1: 'T', 2: 'W', 3: 'K'},
            3: {'id_': '3', 0: 'G', 1: 'T', 2: 'C', 3: 'M'},
            4: {'id_': '4', 0: 'A', 1: 'T', 2: 'G', 3: 'N'},
            }).set_index('id_')
        known = pd.Series(data=[0, 3, 2, 1, 1], 
                          index=pd.Index(['0', '1', '2', '3', '4'], name='id_'), 
                          dtype=int)
        test = _count_degenerates(seq_array)
        pdt.assert_series_equal(known, test)




if __name__ == '__main__':
    main()