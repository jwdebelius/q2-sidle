from unittest import TestCase, main

import dask
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._utils import (_convert_seq_block_to_dna_fasta_format,
                             _convert_generator_to_delayed_seq_block,
                             _convert_generator_to_seq_block,
                             _count_degenerates,
                             )


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

    def test_convert_seq_block_to_dna_fasta_format(self):
        ff = _convert_seq_block_to_dna_fasta_format(self.seq_block)
        test = ff.view(pd.Series)
        pdt.assert_series_equal(self.skbio_series.astype(str), 
                                self.skbio_series.astype(str))

    def test_convert_generator_to_delayed_seq_block(self):
        gen = self.seq_artifact.view(DNAIterator)
        test = _convert_generator_to_delayed_seq_block(gen, chunksize=5)
        self.assertTrue(test, list)
        self.assertEqual(len(test), 1)
        pdt.assert_frame_equal(self.seq_block[0], dask.compute(test[0])[0])

    def test_convert_generator_to_seq_block(self):
        gen = self.seq_artifact.view(DNAIterator)
        test = _convert_generator_to_seq_block(gen, chunksize=5)
        self.assertTrue(test, list)
        self.assertEqual(len(test), 1)
        pdt.assert_frame_equal(self.seq_block[0], test[0])

    def test_count_degenerates(self):
        seq_array = pd.DataFrame.from_dict(orient='index', data={
            0: {'id_': '0', 0: 'A', 1: 'G', 2: 'T', 3: 'C'},
            1: {'id_': '1', 0: 'A', 1: 'R', 2: 'W', 3: 'S'},
            2: {'id_': '2', 0: 'C', 1: 'T', 2: 'W', 3: 'K'},
            3: {'id_': '3', 0: 'G', 1: 'T', 2: 'C', 3: 'M'},
            4: {'id_': '4', 0: 'A', 1: 'T', 2: 'G', 3: 'N'},
            }).set_index('id_')
        known = pd.Series(data=[0, 3, 2, 1, 1], 
                          index=pd.Index(['0', '1', '2', '3', '4'], name='id_'), 
                          dtype=int)
        test = _count_degenerates(seq_array)
        pdt.assert_series_equal(known, test)

        
if __name__ == '__main__':
    main()