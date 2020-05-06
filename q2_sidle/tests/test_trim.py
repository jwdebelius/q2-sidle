from unittest import TestCase, main

import warnings

import biom
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.util.testing as pdt
import skbio
import sparse as sp

from qiime2 import Artifact

from q2_sidle._trim import (trim_dada2_posthoc
                                )
import q2_sidle.tests.test_set as ts

class TrimTest(TestCase):
    def setUp(self):
        self.table = ts.region1_counts.view(biom.Table)
        self.rep_seq = ts.region1_rep_seqs

    def test_trim_dada2_posthoc_min_no_hash(self):
        test_table, test_seqs = trim_dada2_posthoc(self.table, self.rep_seq, 
                                                   hashed_feature_ids=False)
        npt.assert_array_equal(self.rep_seq.view(pd.Series).astype(str).values,
                               test_seqs.view(pd.Series).astype(str).values)
        npt.assert_array_equal(self.rep_seq.view(pd.Series).astype(str).values,
                               test_seqs.view(pd.Series).astype(str).index)
        npt.assert_array_equal(test_table.matrix_data.todense(),
                               self.table.matrix_data.todense())
        npt.assert_array_equal(test_table.ids(axis='sample'),
                               self.table.ids(axis='sample'))
        npt.assert_array_equal(test_table.ids(axis='observation'),
                               self.rep_seq.view(pd.Series).astype(str).values)

    def test_trim_dada2_posthoc_trim_hash(self):
        test_table, test_seqs = trim_dada2_posthoc(self.table, self.rep_seq,
                                                   trim_length=5)
        known_seqs = pd.Series({
            '21f1ded98563f87f5699c8aae7e1639e': 'GCGAA',
            '17e8ebf20d380363085629b22ea39722': 'ATCCG',
            '57a67d2f8243b5c4573798582153bd48': 'TTCCG',
            '74f4d3364948ee9313062e0f0b43a5e9': 'CGTTT',
            })
        known_counts = np.array([[150,  0,   0, 100],
                                 [125, 50,  50,  50],
                                 [100,  0, 100, 100],
                                 ]).T
        pdt.assert_series_equal(known_seqs, 
                                test_seqs.view(pd.Series).astype(str))
        npt.assert_array_equal(known_seqs.index, 
                               test_table.ids(axis='observation'))
        npt.assert_array_equal(known_counts, test_table.matrix_data.todense())
        npt.assert_array_equal(self.table.ids(axis='sample'), 
                               test_table.ids(axis='sample'))


if __name__ == '__main__':
    main()
