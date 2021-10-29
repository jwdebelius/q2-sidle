from unittest import TestCase, main

import os
import tempfile
import warnings

import matplotlib.pyplot as plt
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio

import qiime2
from qiime2 import Artifact, Metadata

from q2_sidle._first_position._visualization import (_build_cover_matrix,
                                                     _make_alignment_heatmap,
                                                     _regional_alignment_results,
                                                     )
from q2_sidle import summarize_alignment_positions

import q2_sidle.tests.test_set as ts

class FirstPositionTest(TestCase):
    def setUp(self):
        self.output_dir_obj = \
            tempfile.TemporaryDirectory(prefix='q2-sidle-summary-')
        self.output_dir = self.output_dir_obj.name

        self.expanded_alignment = ts.extra_alignment.view(pd.Series)
        self.coverage = ts.coverage
        self.summary = pd.DataFrame(
            data=np.array([
                [ 12,  12,  12,  12,  12,  52,  52,  52,   0,  28],
                [300, 100, 100, 100, 100,  20,  10, 100, 525,  13]
                ]).T,
            columns=['starting-position', 'sequence-counts'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                           name='feature-id'),
            )
        self.summary['direction'] = 'fwd'
        self.weighted_coverage = \
            (self.coverage.T * np.log10(self.summary['sequence-counts'])).T
        self.start_ordered = ['asv09', 'asv01', 'asv02', 'asv03', 'asv04',
                              'asv05', 'asv10', 'asv06', 'asv07', 'asv08']
    def tearDown(self):
        self.output_dir_obj.cleanup()

    def assertBasicVizValidity(self, viz_dir, normalize=True):
        index_fp = os.path.join(viz_dir, 'index.html')
        self.assertTrue(os.path.exists(index_fp))

        with open(index_fp) as fh:
            index_html = fh.read()

        # normalize_str = '(normalized)' if normalize else '(not normalized)'
        # self.assertTrue(normalize_str in index_html)

        for ext in ['png', 'svg']:
            fp = os.path.join(viz_dir, 'alignment-start-heatmap.%s' % ext)
            self.assertTrue(os.path.exists(fp))

    def test_summarize_alignment_positions_basic(self):
        summarize_alignment_positions(self.output_dir,
                                      self.expanded_alignment,
                                      Metadata(self.summary),
                                      )
        self.assertBasicVizValidity(self.output_dir)

    def test_build_cover_matrix_error(self):
        self.summary.rename({"asv01": 'grayson',
                             'asv02': 'todd', 
                             'asv03': 'drake', 
                             'asv04': 'brown', 
                             'asv05': 'wayne'},
                       inplace=True)
        with self.assertRaises(ValueError) as err:
            _build_cover_matrix(alignment=self.expanded_alignment,
                                summary=self.summary)
        self.assertEqual(str(err.exception),
                         ('The alignmnent does not contain all summarized '
                         'sequences. Please check that you have the correct '
                         'alignment and summary.')
                         )

    def test_build_cover_matrix_unweighted(self):
        test = _build_cover_matrix(self.expanded_alignment,
                                   summary=self.summary,
                                   weight=False)
        pdt.assert_frame_equal(self.coverage.loc[self.start_ordered],
                               test)

    def test_build_cover_matrix_weighted(self):
        test = _build_cover_matrix(self.expanded_alignment,
                                   summary=self.summary,
                                   weight=True)
        pdt.assert_frame_equal(self.weighted_coverage.loc[self.start_ordered], 
                              test)

    def test_regional_alignment_results(self):
        known = pd.DataFrame(
            data=np.array([['fwd', 1., 525., 525., 525.],
                           ['fwd', 5., 100., 100., 300.],
                           ['fwd', 1.,  13.,  13.,  13.],
                           ['fwd', 3.,  10.,  20., 100.],
                           ], dtype=object),
            index=pd.Index([0, 12, 28, 52], name='Starting Position'),
            columns=['Sequence Direction',
                     'number of mapped ASVs', 
                     'Frequency of ASV with fewest counts', 
                     'Median Frequency of mapped ASVs',
                     'Frequency of ASV with most counts', 
                     ],
            )
        num_cols = ['number of mapped ASVs', 
                     'Frequency of ASV with fewest counts', 
                     'Median Frequency of mapped ASVs',
                     'Frequency of ASV with most counts']
        known[num_cols] = known[num_cols].astype(float)
        test = _regional_alignment_results(self.summary)
        pdt.assert_frame_equal(test, known)

    def test_make_alignment_heatmap_vigilante(self):
        known_mask = self.coverage == 0
        fig, test_mask = _make_alignment_heatmap(self.weighted_coverage,
                                                 maskcolor='r',
                                                 grid=True,
                                                 test=True,
                                                 tick_interval=10,
                                                 )
        pdt.assert_frame_equal(test_mask, known_mask)
        self.assertTrue(isinstance(fig, plt.Figure))
        self.assertEqual(fig.axes[0].get_facecolor(), (1, 0, 0, 1))
        npt.assert_array_equal(np.arange(0, 86, 10), fig.axes[0].get_xticks())

    def test_make_alignment_heatmap_unmasked(self):
        known_mask = self.coverage.isna()
        fig, test_mask = _make_alignment_heatmap(self.weighted_coverage,
                                                 maskcolor=None,
                                                 grid=True,
                                                 test=True,
                                                 )
        pdt.assert_frame_equal(test_mask, known_mask)
        self.assertTrue(isinstance(fig, plt.Figure))
        self.assertEqual(fig.axes[0].get_facecolor(), (1, 1, 1, 1))
        npt.assert_array_equal(np.arange(0, 86, 100), fig.axes[0].get_xticks())


if __name__ == '__main__':
    main()

