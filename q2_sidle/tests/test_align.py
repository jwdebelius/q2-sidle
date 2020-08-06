from unittest import TestCase, main

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._align import (align_regional_kmers,
                             _align_kmers,
                             _check_read_lengths,                        
                             )


class AlignTest(TestCase):
    def setUp(self):
        self.seq_array = pd.Series(['AGTC', 'ARWS', 'CTWK', 'GTCM', 'ATGN'],
                                   index=['A', 'B', 'C', 'D', 'E'])
        self.seq_array.index = self.seq_array.index.astype(str)
        self.in_mer = pd.Series(
            data=np.array(['AGTCCATGC', 'TACGAGTGA', 
                           'ACTCCATGC', 'AAAAAAAGT'])
        )
        self.reads2 = pd.Series(
            data=np.array(['AGTC', 'WGWN', 'AGTT']), 
            index=['r2.0', 'r2.1', 'r2.2'],
            )

    def test_align_kmers_length_error(self):
        with self.assertRaises(ValueError):
            align_regional_kmers(self.seq_array, self.in_mer, region='Gotham')


    def test_align_kmers_with_no_degen(self):
        known0 = pd.DataFrame(
          data=[['A', 'r2.0', 4, 0],
                ],
          columns=['kmer', 'asv', 'length', 'mismatch']
        )
        known0[['mismatch', 'length']] = \
            known0[['mismatch', 'length']].astype(int)
        test0 = _align_kmers(self.seq_array,
                             self.reads2, 
                             allowed_mismatch=0)
        known0.sort_values(['kmer', 'asv'], inplace=True)
        test0.sort_values(['kmer', 'asv'], inplace=True)
        known0.reset_index(drop=True, inplace=True)
        pdt.assert_frame_equal(known0, test0.reset_index(drop=True))

    def test_align_regional_kmers(self):
        kmers = Artifact.import_data('FeatureData[Sequence]', pd.Series({
            'seq1|seq2': skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq1 | seq2'}),
            'seq3@0001': skbio.DNA('ATCCGCGTTGGAGTT', 
                                   metadata={'id': 'seq3@0001'}),
            'seq3@0002': skbio.DNA('TTCCGCGTTGGAGTT', 
                                   metadata={'id': 'seq3@0002'}),
            'seq5': skbio.DNA('CGTTTATGTATGCCC', 
                              metadata={'id': 'seq5'}),
            'seq6':  skbio.DNA('CGTTTATGTATGCCT', 
                              metadata={'id': 'seq6'})
            }))
        rep_set = Artifact.import_data('FeatureData[Sequence]', pd.Series({
            'asv01': skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'asv01'}),
            'asv02': skbio.DNA('ATCCGCGTTGGAGTT', metadata={'id': 'asv02'}),
            'asv03': skbio.DNA('TTCCGCGTTGGAGTT', metadata={'id': 'asv03'}),
            'asv04': skbio.DNA('CGTTTATGTATGCCC', metadata={'id': 'asv04'}),
            'asv05': skbio.DNA('CGTTTATGTATGCCT', metadata={'id': 'asv05'}),
            'asv06': skbio.DNA('AGAGTTTCTGAATCC', metadata={'id': 'asv07'})
            }))

        known = pd.DataFrame(
            data=np.array([['seq1|seq2', 'asv01', 15, 0, 'Bludhaven'],
                           ['seq3@0001', 'asv02', 15, 0, 'Bludhaven'],
                           ['seq3@0001', 'asv03', 15, 1, 'Bludhaven'],
                           ['seq3@0002', 'asv02', 15, 1, 'Bludhaven'],
                           ['seq3@0002', 'asv03', 15, 0, 'Bludhaven'],
                           ['seq5', 'asv04', 15, 0, 'Bludhaven'],
                           ['seq5', 'asv05', 15, 1, 'Bludhaven'],
                           ['seq6', 'asv04', 15, 1, 'Bludhaven'],
                           ['seq6', 'asv05', 15, 0, 'Bludhaven']], 
                           dtype=object),
            columns=['kmer', 'asv', 'length', 'mismatch', 'region'],
            )
        known['length'] = known['length'].astype(int)
        known['mismatch'] = known['mismatch'].astype(int)
        
        match, discard = align_regional_kmers(kmers.view(pd.Series),
                                              rep_set.view(pd.Series),
                                              region='Bludhaven',
                                              debug=True)
        pdt.assert_frame_equal(known, match.reset_index(drop=True))
        known_discard = pd.Series(['AGAGTTTCTGAATCC'], 
                                   pd.Index(['asv06']))
        pdt.assert_series_equal(
            discard.view(pd.Series).astype(str), 
            known_discard, 
        )

    def test_check_read_length_pass(self):
        number_, length_ = _check_read_lengths(self.in_mer, 'inmer')
        self.assertEqual(length_, 9)
        self.assertEqual(number_, 4)

    def test_check_read_length_fail(self):
        reads = pd.concat([self.in_mer, self.seq_array])
        with self.assertRaises(ValueError):
            _check_read_lengths(reads, 'inmer')

if __name__ == '__main__':
    main()
