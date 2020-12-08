from unittest import TestCase, main

import os

import biom
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt

import q2_sidle.tests.test_set as ts

from q2_sidle._accounting import (check_alignment_discard,
                                  track_aligned_counts,
                                  _alignment_accounting
                                  )
from qiime2 import Artifact, Metadata


class ReconstructTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        self.manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(self.base_dir, 'region1_db_map.qza'),
                   os.path.join(self.base_dir, 'region1_align.qza'),
                   os.path.join(self.base_dir, 'region1_counts_discard.qza'), 0],
                  [os.path.join(self.base_dir, 'region2_db_map.qza'),
                   os.path.join(self.base_dir, 'region2_align.qza'),
                   os.path.join(self.base_dir, 'region2_counts.qza'), 1]],
            columns=['kmer-map', 'alignment-map', 'frequency-table', 
                     'region-order'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
            ))
        self.table1 = ts.region1_alt_counts.view(biom.Table)
        self.table2 = ts.region2_counts.view(biom.Table)
        self.align1 = ts.region1_align.view(pd.DataFrame)
        self.align2 = ts.region2_align.view(pd.DataFrame)

    def test_check_alignment_discard(self):
        known = np.array([[50,  25, 0]])

        test = check_alignment_discard(self.align1, self.table1)
        npt.assert_array_equal(test.ids(axis='sample'), 
                               np.array(['sample1', 'sample2', 'sample3']))
        npt.assert_array_equal(test.ids(axis='observation'),
                               ['asv20']
                               )
        npt.assert_array_equal(test.matrix_data.toarray(),
                               known)

    def test_check_alignment_max_discard(self):
        known = np.array([[50, 25, 50],
                          [50,  25, 0]])
        self.align1 = self.align1.loc[self.align1['kmer'] != 'seq5']
        test = check_alignment_discard(
            alignment=self.align1,
            table=self.table1,
            max_mismatch=0,
            )
        npt.assert_array_equal(test.ids(axis='sample'), 
                               np.array(['sample1', 'sample2', 'sample3']))
        npt.assert_array_equal(test.ids(axis='observation'),
                               ['asv04', 'asv20']
                               )
        npt.assert_array_equal(test.matrix_data.toarray(),
                               known)

    def test_alignment_accounting(self):
        known = pd.DataFrame(
                data=np.array([[300, 250, 0.83333333],
                               [300, 275, 0.91666667],
                               [300, 300, 1.00000000]]),
                index=['sample1', 'sample2', 'sample3'],
                columns=['Bludhaven starting counts', 
                         'Bludhaven aligned counts',
                         'Bludhaven aligned percentage']
            )
        test = _alignment_accounting('Bludhaven', self.align1, self.table1)
        pdt.assert_frame_equal(test, known)

    def test_track_sample_counts(self):
        known = pd.DataFrame(
            data=np.array([[600, 550, 0.9166667, 300, 250, 0.833333, 300, 300, 1],
                           [600, 575, 0.9583333, 300, 275, 0.916667, 300, 300, 1],
                           [600, 600, 1, 300, 300, 1, 300, 300, 1]]),
            columns=['total starting counts', 'total aligned counts', 
                     'total aligned percentage', 'Bludhaven starting counts', 
                     'Bludhaven aligned counts',
                     'Bludhaven aligned percentage', 'Gotham starting counts',
                     'Gotham aligned counts', 'Gotham aligned percentage'],
            index=pd.Index(['sample1', 'sample2', 'sample3'], name='sample-id')
            )
        test = track_aligned_counts(
            regions=['Bludhaven', 'Gotham'],
            regional_alignments=[self.align1, self.align2],
            regional_tables=[self.table1, self.table2]
            )
        self.assertTrue(isinstance(test, Metadata))
        pdt.assert_frame_equal(test.to_dataframe().round(5), known.round(5))

if __name__ == '__main__':
    main()