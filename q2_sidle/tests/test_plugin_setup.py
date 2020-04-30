from unittest import TestCase, main

import os

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from qiime2 import Artifact, Metadata
from qiime2.plugins.sidle import methods as sidle
from q2_sidle.plugin_setup import plugin

import q2_sidle.tests.test_set as ts


class PluginSetupTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        self.ref_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'full_db.qza'))
        self.region1_db_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'region1_db_seqs.qza'))
        self.region2_db_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'region2_db_seqs.qza'))
        self.region1_db_map = \
            Artifact.load(os.path.join(self.base_dir, 'region1_db_map.qza'))
        self.region2_db_map = \
            Artifact.load(os.path.join(self.base_dir, 'region2_db_map.qza'))
        self.rep_seqs1 = \
            Artifact.load(os.path.join(self.base_dir, 'region1_rep_set.qza'))
        self.align1 = \
            Artifact.load(os.path.join(self.base_dir, 'region1_align.qza'))
        self.manfiest = Metadata(pd.DataFrame(
            data=[[os.path.join(self.base_dir, 'region1_db_seqs.qza'),
                   os.path.join(self.base_dir, 'region1_align.qza'),
                   os.path.join(self.base_dir, 'region1_counts.qza')],
                  []],
            columns=['kmer-map', 'alignment-map', 'frequency-table'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
            ))
        # self.align2 = \
        #     Artifact.load(os.path.join(self.base_dir, 'region2_align.qza'))
        # self.counts1 = \
        #     Artifact.load(os.path.join(self.base_dir, 'region1_align.qza'))
        # self.counts2 = ts.region2_counts
        #     Artifact.load(os.path.join(self.base_dir, 'region1_rep_set.qza'))
        # self.align1 = \
        #     Artifact.load(os.path.join(self.base_dir, 'region1_align_map.qza'))
        # self.align2 = \
        #     Artifact.load('FeatureData[KmerAlignment]', pd.DataFrame(
        # #         ''
        # #     ))

    
    def test_plugin_setup(self):
        self.assertEqual(plugin.name, 'sidle')

    def test_filter_degenerate_sequences(self):
        known = self.ref_seqs.view(pd.Series).copy().astype(str).drop(['seq3'])

        test = sidle.filter_degenerate_sequences(
            self.ref_seqs, 
            max_degen=0, 
            debug=True
            )
        test = test.filtered_sequences.view(pd.Series).astype(str)
        pdt.assert_series_equal(known, test)

    def test_extract_regional_database(self):
        test_seqs, test_map = \
            sidle.extract_regional_database(self.ref_seqs,
                                            fwd_primer='WANTCAT',
                                            rev_primer='CATCATCAT',
                                            trim_length=15,
                                            primer_mismatch=1,
                                            debug=True,
                                            )
        pdt.assert_series_equal(
            self.region1_db_seqs.view(pd.Series).astype(str),
            test_seqs.view(pd.Series).astype(str)
            )
        pdt.assert_frame_equal(self.region1_db_map.view(pd.DataFrame),
                               test_map.view(pd.DataFrame))

    def test_prepare_extracted_region(self):
        test_seqs, test_map = \
            sidle.prepare_extracted_region(self.region2_db_seqs,
                                           region='Gotham',
                                           trim_length=15,
                                           debug=True,
                                           )
        pdt.assert_series_equal(
            test_seqs.view(pd.Series).astype(str),
            self.region2_db_seqs.view(pd.Series).astype(str)
            )
        pdt.assert_frame_equal(test_map.view(pd.DataFrame), 
                               self.region2_db_map.view(pd.DataFrame))

    def test_align_regional_kmers(self):
        test_align, test_discard = \
            sidle.align_regional_kmers(self.region1_db_seqs,
                                       self.rep_seqs1,
                                       region='WANTCAT-CATCATCAT',
                                       max_mismatch=2,
                                       debug=True,
                                       )
        self.assertEqual(len(test_discard.view(pd.Series)), 0)
        pdt.assert_frame_equal(self.align1.view(pd.DataFrame),
                               test_align.view(pd.DataFrame))

    def test_reconstruct_counts(self):
        
    #     with self.assertRaises(ValueError) as err:
    #         sidle.reconstruct_counts(kmer_map=self.region1_db_map,
    #                                  kmer_map=self.region2_db_map,
    #                                  regional_alignment=self.align1,
    #                                  count_table=self.count1,
    #                                  count_table=self.count2,
    #                                  )
    #     self.assertEqual(str(err.exception), 
    #                      ('Each database region must have a corresonding '
    #                       'alignment.'))


# class IntegrationTest(TestCase):


if __name__ == '__main__':
    main()