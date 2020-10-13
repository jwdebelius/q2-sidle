from unittest import TestCase, main

import os

import dask
import numpy as np
import pandas as pd
import pandas.util.testing as pdt
from skbio import DNA

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact, Metadata
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._tree import (reconstruct_fragment_rep_seqs,
                            _expand_primer,
                            _find_exact_forward,
                            _find_exact_reverse,
                            _find_approx_forward,
                            _find_approx_reverse,
                            _group_concensus,
                            )
from q2_sidle.tests import test_set as ts

class TreeTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        self.aligned_seqs = pd.Series({
            'seq01': DNA('-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC---------------------------------'),
            'seq02': DNA('ACTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC---------------------------------'),
            'seq03': DNA('CATAGTCATWTCCGCGTTGGAGTTATGATGATGAWACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCA-'),
            'seq04': DNA('------------------GGAGTTATGATGA--AGACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCAC'),
            'seq05': DNA('CATAGTCATCGTTTATGTATGCCCATGATGATGCGAGCACCTCGTATGGATGTAGAGCCACTGACGTGCGG'),
        })
        kmer1 = Artifact.load(os.path.join(self.base_dir, 'frag_r1_db_map.qza'))
        kmer2 = Artifact.load(os.path.join(self.base_dir, 'frag_r2_db_map.qza'))
        self.kmer_map1 = kmer1.view(pd.DataFrame)
        self.kmer_map2 = kmer2.view(pd.DataFrame)
        np.random.seed(5)

    def test_reconstruct_fragment_rep_seqs(self):
        recon_map = pd.Series(
            data=['seq01|seq02', 'seq01|seq02', 'seq03|seq04', 
                  'seq03|seq04', 'seq05'],
            index=pd.Index(['seq01', 'seq02', 'seq03', 'seq04', 'seq05'], 
                            name='db-seq'),
            name='clean_name'
            )

        recon_summary = pd.DataFrame(
            data=[[1, 2, 2, 0, 'asv01|asv02'],
                  [2, 3, 1.5, np.std([1, 2], ddof=1), 'asv03|asv04'],
                  [2, 2, 1, 0, 'asv07|asv08']],
            index=pd.Index(['seq01|seq02', 'seq03|seq04', 'seq05'], 
                            name='clean_name'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
            )
        known = pd.Series(
            data=['GCGAAGCGGCTCAGG',
                  'WTCCGCGTTGGAGTTATGATGATGAGACCACCTCGTCCCAGTTCCGCGCTTC'],
            index=pd.Index(['seq01|seq02', 'seq03|seq04'], name='clean_name'),
            )
        test = reconstruct_fragment_rep_seqs(
            region=['Bludhaven', 'Gotham'],
            kmer_map=[self.kmer_map1, self.kmer_map2],
            reconstruction_map=recon_map,
            reconstruction_summary=recon_summary,
            aligned_sequences=self.aligned_seqs,
            )
        pdt.assert_series_equal(test.view(pd.Series).astype(str), known)

    def test_expand_primer_miss(self):
        primer = 'WANTCAT'
        known = '((?<=([AT])))((A[ACGT]TCAT){e<=1})'
        test = _expand_primer(primer, 1)
        self.assertEqual(known, test)

    def test_expand_primer_none(self):
        primer = 'WANTCAT'
        known = '[AT]A[ACGT]TCAT'
        test = _expand_primer(primer, None)
        self.assertEqual(known, test)

    def test_find_exact_forward_match(self):
        args = pd.Series([
            '-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC--------------',
            '[AT]A[ACGT]TCAT'])
        test = _find_exact_forward(args)
        self.assertEqual(test, 9)

    def test_find_exact_forward_miss(self):
        args = pd.Series([
            '-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC--------------',
            'WANTCAT'])
        test = _find_exact_forward(args)
        self.assertTrue(np.isnan(test))

    def test_find_extract_reverse_match(self):
        args = pd.Series([
            '-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC--------------',
            'ATGATGATG'])
        test = _find_exact_reverse(args)
        self.assertTrue(test, 24)

    def test_find_extract_reverse_match(self):
        args = pd.Series([
            '-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC--------------',
            'CATCATCAT'])
        test = _find_exact_reverse(args)
        self.assertTrue(test, np.nan)    

    def test_find_approx_forward(self):
        args = pd.Series([
            DNA('-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC--------------'),
            'WANTCAT'])
        test = _find_approx_forward(args)
        self.assertEqual(test, 9)

    def test_find_approx_reverse(self):
        args = pd.Series([
            DNA('-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC--------------'),
            'ATGATGATG'])
        test = _find_approx_reverse(args)
        self.assertEqual(test, 24)

    def test_get_consus_seq(self):
        g = self.aligned_seqs.loc[['seq03', 'seq04']]

        known = ('CATAGTCATWTCCGCGTTGGAGTTATGATGATGAGACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCAC')

        test = _group_concensus(g)
        self.assertEqual(known, str(test))


if __name__ == '__main__':
    main()