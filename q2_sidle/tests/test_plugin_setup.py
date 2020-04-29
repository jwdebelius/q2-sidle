from unittest import TestCase, main

import os

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from qiime2 import Artifact, Metadata
from q2_sidle.plugin_setup import plugin


class PluginSetupTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        self.ref_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'full_db.qza'))
        self.region1_extract = \
            Artifact.load(os.path.join(self.base_dir, 'region1_db_seqs.qza'))
        self.region2_extract = \
            Artifact.load(os.path.join(self.base_dir, 'region2_db_seqs.qza'))
    
    def test_plugin_setup(self):
        self.assertEqual(plugin.name, 'sidle')

    def test_filter_degenerate_sequences(self):
        known = self.ref_seqs.view(pd.Series).copy().astype(str).drop(['seq3'])

        test = plugin.plugin.methods.filter_degenerate_sequences(
            self.ref_seqs, 
            max_degen=0, 
            debug=True
            )
        test = test.filtered_sequences.view(pd.Series).astype(str)
        pdt.assert_series_equal(known, test)

    # def test_extract_regional_database(self):
    #     known = ref_1 = pd.Series(
    #         data=[('GCGAAGCGGCTCAGG'), 
    #               ('ATCCGCGTTGGAGTT'),
    #               ('TTCCGCGTTGGAGTT'),
    #               ('CGTTTATGTATGCCC'),
    #               ('CGTTTATGTATGCCT')],
    #     index=['seq1 | seq2', 'seq3@0001', 'seq3@0002', 'seq5', 'seq6'],
    #     )
    #     test_seqs, test_map = extract_regional_database(self.ref_seqs,
    #                                                     fwd_primer='WANTCAT',
    #                                                     rev_primer='CATCATCAT',
    #                                                     trim_length=15,
    #                                                     debug=True,
    #                                                     )


# class IntegrationTest(TestCase):


if __name__ == '__main__':
    main()