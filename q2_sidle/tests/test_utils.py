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

from q2_sidle._utils import (_count_degenerates,
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

    def test_count_degenerates(self):
        seq_array = pd.DataFrame.from_dict(orient='index', data={
            0: {'id_': '0', 0: 'A', 1: 'G', 2: 'T', 3: 'C'},
            1: {'id_': '1', 0: 'A', 1: 'R', 2: 'W', 3: 'S'},
            3: {'id_': '3', 0: 'G', 1: 'T', 2: 'C', 3: 'M'},
            4: {'id_': '4', 0: 'A', 1: 'T', 2: 'G', 3: 'N'},
            }).set_index('id_')
        known = pd.Series(data=[0, 3, 1, 1], 
                          index=pd.Index(['0', '1', '3', '4'], name='id_'), 
                          dtype=int)
        test = _count_degenerates(seq_array)
        pdt.assert_series_equal(known, test)

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