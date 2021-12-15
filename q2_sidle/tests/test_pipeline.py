from unittest import TestCase, main

from copy import copy
import os
import shutil
import warnings


import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio
from skbio import DNA

from qiime2 import Artifact, Metadata
from qiime2.plugins.sidle import actions as sidle
from q2_sidle.plugin_setup import plugin

import q2_sidle.tests.test_set as ts


class PipelineTest(TestCase):
	setUp(self):
        self.align1 = ts.region1_align
        self.align2 = ts.region2_align
        self.kmer_map1 = ts.region1_db_map
        self.kmer_map2 = ts.region2_db_map
        self.table1 = ts.region1_counts
        self.table2 = ts.region2_counts
        self.taxonomy = ts.taxonomy
        self.seq_map = ts.seq_map
        self.database_summary = ts.db_summary

    def test_sidle_reconstruction(self):
        mapping, summary, counts, taxonomy  = sidle.sidle_reconstruction(
            region=['Bludhaven', 'Gotham'],
            kmer_map=[self.kmer_map1, self.kmer_map2],
            regional_alignment=[self.align1, self.align2],
            regional_table=[self.table1, self.table2],
            reference_taxonomy=self.taxonomy,
            database='greengenes',
            define_missing='ignore',
            min_counts=10,
            debug=True
            )
        pdt.assert_frame_equal(mapping.view(pd.DataFrame), 
                               self.seq_map.view(pd.DataFrame))
        pdt.assert_frame_equal(summary.view(pd.DataFrame),
                               self.database_summary.view(pd.DataFrame))
        pdt.assert_frame_equal(
            counts.view(pd.DataFrame),
            pd.DataFrame( 
                data=np.array([[100.,  50,   0,  50,  50, 50],
                               [100.,  25, 100,  25,  25, 25],
                               [  0., 100, 100,   0,  50, 50]]),
                index=pd.Index(['sample1', 'sample2', 'sample3']),
                columns=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6']
            )
        )
        pdt.assert_series_equal(self.taxonomy.view(pd.Series),
                                taxonomy.view(pd.Series))

    # def test 

if __name__ == '__main__':
    main()
