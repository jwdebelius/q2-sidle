from unittest import TestCase, main

import os
import shutil

import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from qiime2.plugin import ValidationError
from q2_sidle._formats import (KmerMapFormat,
                               KmerMapDirFmt,
                               KmerAlignFormat,
                               KmerAlignDirFmt,
                               SidleReconFormat,
                               SidleReconDirFormat,
                               ReconSummaryFormat,
                               ReconSummaryDirFormat,
                               )

class PluginSetupTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/types')
        self.tmp = os.path.join(self.base_dir, 'tmp')
        os.makedirs(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_kmer_map_format_validate_pass(self):
        filepath = os.path.join(self.base_dir, 'kmer-map.tsv')
        format = KmerMapFormat(filepath, mode='r')
        format.validate()

    def test_kmer_map_format_validate_column_fail(self):
        filepath = os.path.join(self.base_dir, 'kmer-map-col-fail.tsv')
        format = KmerMapFormat(filepath, mode='r')
        with self.assertRaises(ValidationError) as err:
            format.validate()
        self.assertEqual(str(err.exception), 
                         'The KmerMap does not contain the correct columns')

    def test_kmer_map_format_validate_column_kmer_fail(self):
        filepath = os.path.join(self.base_dir, 'kmer-map-kmer-fail.tsv')
        format = KmerMapFormat(filepath, mode='r')
        with self.assertRaises(ValidationError) as err:
            format.validate()
        self.assertEqual(str(err.exception), 'The kmer-length column must be numeric')

    def test_kmer_map_dir_format_validate_pass(self):
        shutil.copy(os.path.join(self.base_dir, 'kmer-map.tsv'), self.tmp)
        format = KmerMapDirFmt(self.tmp, mode='r')
        format.validate()

    def test_kmer_align_format(self):
        filepath = os.path.join(self.base_dir, 'kmer-align.tsv')
        format = KmerAlignFormat(filepath, mode='r')
        format.validate()

    def test_kmer_map_dir_format_validate(self):
        shutil.copy(os.path.join(self.base_dir, 'kmer-align.tsv'), self.tmp)
        format = KmerAlignDirFmt(self.tmp, 'r')
        format.validate()

    def test_sidle_recon_format(self):
        filepath = os.path.join(self.base_dir, 
                                'sidle-reconstruction-mapping.tsv')
        format = SidleReconFormat(filepath, mode='r')
        format.validate()

    def test_sidle_recon_dir_format_validate(self):
        shutil.copy(os.path.join(self.base_dir, 
                                 'sidle-reconstruction-mapping.tsv'), 
                    self.tmp)
        format = SidleReconDirFormat(self.tmp, 'r')
        format.validate()

    def test_sidle_summary_format(self):
        filepath = os.path.join(self.base_dir, 'sidle-summary.tsv')
        format = ReconSummaryFormat(filepath, mode='r')
        format.validate()

    def test_sidle_summary_dir_format_validate(self):
        shutil.copy(os.path.join(self.base_dir, 'sidle-summary.tsv'), 
                    self.tmp)
        format = ReconSummaryDirFormat(self.tmp, 'r')
        format.validate()



if __name__ == '__main__':
    main()