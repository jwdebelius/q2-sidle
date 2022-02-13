from unittest import TestCase, main

import os

import dask
import numpy as np
import pandas as pd
import pandas.testing as pdt
from skbio import DNA

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact, Metadata
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._tree import (reconstruct_fragment_rep_seqs,
                            _find_approx_forward,
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

        self.recon_map = pd.DataFrame(
            data=np.array([['seq01|seq02', 0, 'WANTCAT', 9, 0, 'WANTCAT', np.nan, 15], 
                           ['seq01|seq02', 0, 'WANTCAT', 9, 0, 'WANTCAT', np.nan, 15], 
                           ['seq03|seq04', 0, 'WANTCAT', 9, 1, 'CACCTCGTN', np.nan, 15], 
                           ['seq03|seq04', 0, 'CACCTCGTN', np.nan, 1, 'CACCTCGTN', np.nan, 15], 
                           ['seq05', 0, 'WANTCAT', 9, 1, 'CACCTCGTN', np.nan, 15],
                           ],  dtype=object),
            index=pd.Index(['seq01', 'seq02', 'seq03', 'seq04', 'seq05'], name='db-seq'),
            columns=['clean_name', 'first-region', 'first-fwd-primer',  'first-fwd-pos',
                     'last-region', 'last-fwd-primer', 'last-fwd-pos', 'last-kmer-length'],
            )
        self.recon_summary = pd.DataFrame(
            data=[[1, 2, 2, 0, 'asv01|asv02'],
                  [2, 3, 1.5, np.std([1, 2], ddof=1), 'asv03|asv04'],
                  [2, 2, 1, 0, 'asv07|asv08']],
            index=pd.Index(['seq01|seq02', 'seq03|seq04', 'seq05'], 
                            name='clean_name'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
            )
        np.random.seed(5)

    def test_reconstruct_fragment_rep_seqs_all_unique(self):
        recon_summary = pd.DataFrame(
            data=[[1, 2, 2, 0, 'asv01|asv02'],
                  [2, 3, 1.5, np.std([1, 2], ddof=1), 'asv03|asv04'],
                  [2, 2, 1, 0, 'asv07|asv08']],
            index=pd.Index(['seq01', 'seq03', 'seq05'], 
                            name='clean_name'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
            )
        test = reconstruct_fragment_rep_seqs(self.recon_map, recon_summary, self.aligned_seqs)

    def test_reconstruct_fragment_rep_seqs(self):
        known = pd.Series(
            data=['GCGAAGCGGCTCAGG',
                  'WTCCGCGTTGGAGTTATGATGATGAGACCACCTCGTCCCAGTTCCGCGCTT'],
            index=pd.Index(['seq01|seq02', 'seq03|seq04'], name='clean_name'),
            )
        test = reconstruct_fragment_rep_seqs(
            reconstruction_map=self.recon_map,
            reconstruction_summary=self.recon_summary,
            aligned_sequences=self.aligned_seqs,
            )
        pdt.assert_series_equal(test.view(pd.Series).astype(str), known)

    def test_find_approx_forward(self):
        args = pd.Series([
            DNA('-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC--------------'),
            'WANTCAT'])
        test = _find_approx_forward(args)
        self.assertEqual(test, 9)

    def test_get_consus_seq(self):
        g = self.aligned_seqs.loc[['seq03', 'seq04']]

        known = ('CATAGTCATWTCCGCGTTGGAGTTATGATGATGAGACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCAC')

        test = _group_concensus(g)
        self.assertEqual(known, str(test))


if __name__ == '__main__':
    main()