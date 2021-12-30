from unittest import TestCase, main

import os

import dask
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio

from qiime2 import Artifact, Metadata
from qiime2.plugin import ValidationError

from q2_sidle._utils import (
                             _check_regions,
                             )
import q2_sidle.tests.test_set as ts
from q2_types.feature_data import DNAIterator, DNAFASTAFormat



class UtilTest(TestCase):
    def setUp(self):
        self.seq_block = [pd.DataFrame(
            data=[list('CATS'), list('WANT'), list("CANS")],
            index=['0', '1', '2']
        )]
        self.skbio_series = pd.Series(data={
            "0": skbio.DNA('CATS', metadata={'id': '0'}),
            "1": skbio.DNA('WANT', metadata={'id': '1'}),
            "2": skbio.DNA('CANS', metadata={'id': '2'}),
             })
        self.seq_artifact = Artifact.import_data('FeatureData[Sequence]', 
                                                 self.skbio_series, 
                                                 pd.Series)

    def test_check_regions(self):
        regions = ['Bludhaven', 'Gotham']
        k_order = {'Bludhaven': 0, 'Gotham': 1}
        k_names = {0: 'Bludhaven', 1: "Gotham"}

        t_order, t_names, t_num = _check_regions(regions)
        self.assertEqual(t_num, 2)
        npt.assert_array_equal(list(k_order.keys()), list(t_order.keys()))
        for k, v in k_order.items():
            self.assertEqual(v, t_order[k])
        npt.assert_array_equal(list(k_names.keys()), list(t_names.keys()))
        for k, v in k_names.items():
            self.assertEqual(v, t_names[k])



if __name__ == '__main__':
    main()