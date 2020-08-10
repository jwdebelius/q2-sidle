from unittest import TestCase, main

import os

import dask
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact
from q2_sidle._filter_seqs import (filter_degenerate_sequences,
                                   )

class FilterTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        self.ref_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'full_db.qza'))

    def test_filter_degenerate_sequences(self):
      # Tested in plugin_setup.py
      pass


if __name__ == '__main__':
    main()