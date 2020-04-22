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

        
if __name__ == '__main__':
    main()