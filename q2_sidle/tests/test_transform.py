from unittest import TestCase, main

import os

from dask.delayed import Delayed
import dask.dataframe as dd
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio

from qiime2 import Metadata
from qiime2.plugin.testing import TestPluginBase
from q2_sidle import (KmerMapFormat,
                      KmerAlignFormat,
                      SidleReconFormat,
                      ReconSummaryFormat,
                      AlignmentLedgerFormat,
                      )
import q2_sidle._transformer as t

class TestTransform(TestCase):
    def setUp(self):
        ### The IO transformations get tested in the plugin set up and in the
        ### format testing. So, I just want to make sure that I can import/
        ### export formats correcting and I get those transformations correct.
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/types')

    def test_kmer_map_to_dataframe(self):
        known = pd.DataFrame(
            data=[['Batman', 'Batman', 'Gotham', 'WANTCAT', 'CATCATCAT', 50],
                  ['Superman', 'Superman', 'Metropolis', 'CATDAD', 'DADCAT', 
                   50]],
            columns=['seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer',
                     'kmer-length'],
            index=pd.Index(['Batman', 'Superman'], name='db-seq')
            )
        filepath = os.path.join(self.base_dir, 'kmer-map.tsv')
        format = KmerMapFormat(filepath, mode='r')
        test = t._1(format)
        self.assertTrue(isinstance(test, pd.DataFrame))
        pdt.assert_frame_equal(known, test)

    def test_kmer_map_to_metadata(self):
        known = pd.DataFrame(
            data=[['Batman', 'Batman', 'Batman', 'Gotham', 'WANTCAT', 
                   'CATCATCAT', 50.],
                  ['Superman', 'Superman', 'Superman', 'Metropolis', 'CATDAD',
                   'DADCAT', 50.]],
            columns=['db-seq', 'seq-name', 'kmer', 'region', 'fwd-primer', 
                     'rev-primer', 'kmer-length'],
            index=pd.Index(['0', '1'], name='id')
        )
        filepath = os.path.join(self.base_dir, 'kmer-map.tsv')
        format = KmerMapFormat(filepath, mode='r')
        test = t._2(format)
        self.assertTrue(isinstance(test, Metadata))
        columns = dict(test.columns)
        npt.assert_array_equal(list(columns.keys()),
                              ['db-seq', 'seq-name', 'kmer', 'region', 
                              'fwd-primer', 'rev-primer', 'kmer-length']
                               )
        for k, v in columns.items():
            if k == 'kmer-length':
                self.assertEqual(v.type, 'numeric')
            else:
                self.assertEqual(v.type, 'categorical')
        pdt.assert_frame_equal(known, test.to_dataframe())
        
    def test_kmer_map_delayed(self):
        known = pd.DataFrame(
            data=[['Batman', 'Batman', 'Gotham', 'WANTCAT', 'CATCATCAT', 50],
                  ['Superman', 'Superman', 'Metropolis', 'CATDAD', 'DADCAT', 
                   50]],
            columns=['seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer',
                     'kmer-length'],
            index=pd.Index(['Batman', 'Superman'], name='db-seq')
            )
        filepath = os.path.join(self.base_dir, 'kmer-map.tsv')
        format = KmerMapFormat(filepath, mode='r')
        test = t._3(format)
        self.assertTrue(isinstance(test, Delayed))
        pdt.assert_frame_equal(known, test.compute())
        
    def test_dataframe_to_kmer_map(self):
        input_ = pd.DataFrame(
            data=[['Batman', 'Batman', 'Gotham', 'WANTCAT', 'CATCATCAT', 
                    np.nan, np.nan, 50],
                  ['Superman', 'Superman', 'Metropolis', 'CATDAD', 'DADCAT', 
                   np.nan, np.nan, 50]],
            columns=['seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer',
                     'fwd-pos', 'rev-pos', 'kmer-length'],
            index=pd.Index(['Batman', 'Superman'], name='db-seq')
            )
        test = t._4(input_)
        self.assertTrue(isinstance(test, KmerMapFormat))

    def test_kmer_align_to_dataframe(self):
        known = pd.DataFrame(
            data=[['Batman', 'Bruce Wayne', 80, 2, 2, 'Gotham'],
                  ['Flash', 'Barry Allen', 50, 0, 2, 'Central City'],
                  ['GreenArrow', 'Oliver Queen', 30, 0, 2, 'Star City']],
            columns=['kmer', 'asv', 'length', 'mismatch',  
                     'max-mismatch', 'region']
            )
        filepath = os.path.join(self.base_dir, 'kmer-align.tsv')
        format = KmerAlignFormat(filepath, mode='r')
        test = t._5(format)
        self.assertTrue(isinstance(test, pd.DataFrame))
        pdt.assert_frame_equal(test, known)

    def test_kmer_align_to_metadata(self):
        known = pd.DataFrame(
            data=[['Batman', 'Bruce Wayne', 80., 2., 2., 'Gotham'],
                  ['Flash', 'Barry Allen', 50, 0, 2, 'Central City'],
                  ['GreenArrow', 'Oliver Queen', 30, 0, 2, 'Star City']],
            columns=['kmer', 'asv', 'length', 'mismatch', 
                     'max-mismatch', 'region']
            )
        known.index.set_names('feature-id', inplace=True)
        known.index = known.index.astype(str)
        filepath = os.path.join(self.base_dir, 'kmer-align.tsv')
        format = KmerAlignFormat(filepath, mode='r')
        test = t._6(format)
        self.assertTrue(isinstance(test, Metadata))
        columns = dict(test.columns)
        npt.assert_array_equal(
            list(columns.keys()),
            ['kmer', 'asv', 'length',  'mismatch', 
             'max-mismatch', 'region']
            )
        for k, v in columns.items():
            if k in {'mismatch', 'length', 'max-mismatch'}:
                self.assertEqual(v.type, 'numeric')
            else:
                self.assertEqual(v.type, 'categorical')
        pdt.assert_frame_equal(test.to_dataframe(), 
                               known)

    def test_kmer_align_to_delayed(self):
        known = pd.DataFrame(
            data=[['Batman', 'Bruce Wayne', 80, 2, 2, 'Gotham'],
                  ['Flash', 'Barry Allen', 50, 0, 2, 'Central City'],
                  ['GreenArrow', 'Oliver Queen', 30, 0, 2, 'Star City']],
            columns=['kmer', 'asv', 'length', 'mismatch', 
                     'max-mismatch', 'region']
            )
        filepath = os.path.join(self.base_dir, 'kmer-align.tsv')
        format = KmerAlignFormat(filepath, mode='r')
        test = t._7(format)
        self.assertTrue(isinstance(test, Delayed))
        pdt.assert_frame_equal(test.compute(), known)

    def test_dataframe_to_kmer_align(self):
        input_ = known = pd.DataFrame(
            data=[['Batman', 'Bruce Wayne', 80, 2, 2, 'Gotham'],
                  ['Flash', 'Barry Allen', 50, 0, 2, 'Central City'],
                  ['GreenArrow', 'Oliver Queen', 30, 0, 2, 'Star City']],
            columns=['kmer', 'asv', 'length', 'mismatch', 
                     'max-mismatch', 'region']
            )
        test = t._8(input_)
        self.assertTrue(isinstance(test, KmerAlignFormat))

    def test_recon_map_to_series(self):
        known = pd.Series({"Batman": 'Batman', 
                           'WonderWoman': 'WonderWoman', 
                           'Superman': 'Superman'},
                           name='clean_name')
        known.index.set_names('db-seq', inplace=True)

        filepath = os.path.join(self.base_dir, 
                                'sidle-reconstruction-mapping.tsv')
        format = SidleReconFormat(filepath, mode='r')
        test = t._9(format)
        pdt.assert_series_equal(test, known)

    def test_recon_map_to_dataframe(self):
        known = pd.DataFrame(
            data=[['Batman', 80., 50., 'BATDAD', 'DADBAT'],
                  ['WonderWoman', 60, 20, 'BATDAD', 'DADBAT'],
                  ['Superman', 90, 50, 'BATDAD', 'DADBAT']],
            columns=['clean_name', 'length-solo', 'length-justice-league', 
                     'fwd-primer', 'rev-primer'],
            index=pd.Index(['Batman', 'WonderWoman', 'Superman'],
                           name='db-seq')
        )
        filepath = os.path.join(self.base_dir, 
                                'sidle-reconstruction-mapping.tsv')
        format = SidleReconFormat(filepath, mode='r')
        test = t._10(format)
        pdt.assert_frame_equal(test, known)

    def test_dataframe_recon_map(self):
        input_ = pd.DataFrame(
            data=[['Batman', 80., 50., 'BATDAD', 'DADBAT'],
                  ['WonderWoman', 60, 20, 'BATDAD', 'DADBAT'],
                  ['Superman', 90, 50, 'BATDAD', 'DADBAT']],
            columns=['clean_name', 'length-solo', 'length-justice-league', 
                     'fwd-primer', 'rev-primer'],
            index=pd.Index(['Batman', 'WonderWoman', 'Superman'],
                           name='db-seq')
        )
        test = t._11(input_)
        self.assertTrue(isinstance(test, SidleReconFormat))

    def test_recon_summary_to_dataframe(self):
        known = pd.DataFrame(
            data=[[2, 4, 2., 0., 'Bruce|Wayne|Richard|Grayson'],
                  [1, 1, 1., 0., 'Kent']],
            index=pd.Index(['Batman', 'Superman'], 
                            name='feature-id'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
        )
        filepath = os.path.join(self.base_dir, 'sidle-summary.tsv')
        format = ReconSummaryFormat(filepath, mode='r')
        test = t._12(format)
        pdt.assert_frame_equal(test, known)

    def test_recon_summary_to_metadata(self):
        known = pd.DataFrame(
            data=[[2., 4., 2., 0., 'Bruce|Wayne|Richard|Grayson'],
                  [1., 1., 1., 0., 'Kent']],
            index=pd.Index(['Batman', 'Superman'], 
                            name='feature-id'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
        )
        filepath = os.path.join(self.base_dir, 'sidle-summary.tsv')
        format = ReconSummaryFormat(filepath, mode='r')
        test = t._13(format)
        self.assertTrue(isinstance(test, Metadata))
        columns = dict(test.columns)
        npt.assert_array_equal(list(columns.keys()),
                              ['num-regions', 'total-kmers-mapped', 
                               'mean-kmer-per-region', 'stdv-kmer-per-region', 
                               'mapped-asvs']
                               )
        for k, v in columns.items():
            if k == 'mapped-asvs':
                self.assertEqual(v.type, 'categorical')
            else:
                self.assertEqual(v.type, 'numeric')
        pdt.assert_frame_equal(test.to_dataframe(), known)

    def test_metadata_to_recon_summary(self):
        input_ = Metadata(pd.DataFrame(
            data=[[2., 4., 2., 0., 'Bruce|Wayne|Richard|Grayson'],
                  [1., 1., 1., 0., 'Kent']],
            index=pd.Index(['Batman', 'Superman'], 
                            name='feature-id'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
        ))
        test = t._14(input_)
        self.assertTrue(isinstance(test, ReconSummaryFormat))

    def test_delayed_df_recon_summary(self):
        input_ = dd.from_pandas(pd.DataFrame(
            data=[[2., 4., 2., 0., 'Bruce|Wayne|Richard|Grayson'],
                  [1., 1., 1., 0., 'Kent']],
            index=pd.Index(['Batman', 'Superman'], 
                            name='feature-id'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
        ), npartitions=1)
        test = t._15(input_)
        self.assertTrue(isinstance(test, KmerAlignFormat))

    def test_alignment_pos_to_metadata(self):
        known = pd.DataFrame(
            data=np.vstack([
                np.hstack([np.array([12.] * 5), np.array([52.] * 3), 0., 28.]),
                np.array([101.] * 10),
                ]).T,
            columns=['starting-position', 'sequence-counts'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                           name='feature-id'),
            )
        filepath = os.path.join(self.base_dir, 'position-summary.tsv')
        format = AlignmentLedgerFormat(filepath, mode='r')
        test = t._16(format)
        self.assertTrue(isinstance(test, Metadata))
        columns = dict(test.columns)
        npt.assert_array_equal(list(columns.keys()),
                              ['starting-position', 'sequence-counts'])
        for k, v in columns.items():
            self.assertEqual(v.type, 'numeric')
        pdt.assert_frame_equal(test.to_dataframe(), known)

    def test_metadata_to_alignment_pos(self):
        input_ = Metadata(pd.DataFrame(
            data=np.vstack([
                np.hstack([np.array([12.] * 5), np.array([52.] * 3), 0., 28.]),
                np.array([101.] * 10),
                ]).T,
            columns=['starting-position', 'sequence-counts'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                           name='feature-id'),
            ))
        test = t._17(input_)
        self.assertTrue(isinstance(test, AlignmentLedgerFormat))


if __name__ == '__main__':
    main()
