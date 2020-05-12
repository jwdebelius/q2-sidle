from unittest import TestCase, main

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._align import (align_regional_kmers,
                             _align_kmers,                             
                             )


class AlignTest(TestCase):
    def setUp(self):
        self.seq_array = pd.DataFrame.from_dict(orient='index', data={
            0: {'id_': '0', 0: 'A', 1: 'G', 2: 'T', 3: 'C'},
            1: {'id_': '1', 0: 'A', 1: 'R', 2: 'W', 3: 'S'},
            2: {'id_': '2', 0: 'C', 1: 'T', 2: 'W', 3: 'K'},
            3: {'id_': '3', 0: 'G', 1: 'T', 2: 'C', 3: 'M'},
            4: {'id_': '4', 0: 'A', 1: 'T', 2: 'G', 3: 'N'},
            }).set_index('id_')
        self.in_mer = pd.DataFrame(
            data=np.array([list('AGTCCATGC'), 
                          list('TACGAGTGA'),
                          list('ACTCCATGC'),
                          list('AAAAAAAGT')]),
            )
        self.reads2 = pd.DataFrame(
            data=np.array([list('AGTC'), list('WGWN'), list('AGTT')]), 
            index=['r2.0', 'r2.1', 'r2.2'],
            )

    def test_align_kmers_length_error(self):
        with self.assertRaises(ValueError):
            _align_kmers(self.seq_array, self.in_mer)

    def test_align_kmers_with_degen_expand(self):
        # We'll use seq_array as reads1
        known0 = pd.DataFrame(
            data=[['0', '0', '1', '1'],
                  ['r2.0', 'r2.1', 'r2.0', 'r2.1'],
                  [0, 0, 0, 0],
                  [4, 4, 4, 4],
                ],
            index=['kmer', 'asv', 'mismatch', 'length'],
            ).T
        known0[['mismatch', 'length']] = \
            known0[['mismatch', 'length']].astype(int)
        test0 = _align_kmers(self.seq_array,
                             self.reads2, 
                             allow_degen1=True, 
                             allowed_mismatch=0, 
                             allow_degen2=True, 
                             expand_match=True)
        pdt.assert_frame_equal(known0, test0.reset_index(drop=True))

    def test_align_kmers_with_degen_no_exp(self):
        known0 = pd.DataFrame(
            data=[['0', '0', '1'],
                  ['r2.0', 'r2.1', 'r2.0'],
                  [0, 0, 0],
                  [4, 4, 4],
                  ],
            index=['kmer', 'asv', 'mismatch', 'length']
            ).T
        known0[['mismatch', 'length']] = \
            known0[['mismatch', 'length']].astype(int)
        test0 = _align_kmers(self.seq_array,
                             self.reads2, 
                             allow_degen1=True,
                             allow_degen2=True,
                             allowed_mismatch=0)
        
    def test_align_kmers_with_one_degen(self):
        known0 = pd.DataFrame(
            data=[['0', '1'],
                  ['r2.0', 'r2.0'],
                  [0, 0],
                  [4, 4]
                ],
            index=['kmer', 'asv', 'mismatch', 'length']
            ).T
        known0[['mismatch', 'length']] = \
            known0[['mismatch', 'length']].astype(int)
        test0 = _align_kmers(self.seq_array,
                             self.reads2, 
                             allow_degen1=True,
                             allow_degen2=False, 
                             allowed_mismatch=0)
        pdt.assert_frame_equal(known0, test0.reset_index(drop=True))

    def test_align_kmers_with_no_degen(self):
        known0 = pd.DataFrame(
            data=[['0'],
                  ['r2.0'],
                  [0],
                  [4]
                ],
            index=['kmer', 'asv', 'mismatch', 'length']
            ).T
        known0[['mismatch', 'length']] = \
            known0[['mismatch', 'length']].astype(int)
        test0 = _align_kmers(self.seq_array,
                             self.reads2, 
                             allow_degen1=False, 
                             allowed_mismatch=0)
        pdt.assert_frame_equal(known0, test0.reset_index(drop=True))

    def test_align_regional_kmers(self):
        kmers = Artifact.import_data('FeatureData[Sequence]', pd.Series({
            'seq1|seq2': skbio.DNA('GCGAAGCGGCTCAGG', 
                                     metadata={'id': 'seq1 | seq2'}),
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
            data=np.array([['seq1|seq2', 'asv01', 0, 15, 'Bludhaven'],
                           ['seq3@0001', 'asv02', 0, 15, 'Bludhaven'],
                           ['seq3@0001', 'asv03', 1, 15, 'Bludhaven'],
                           ['seq3@0002', 'asv02', 1, 15, 'Bludhaven'],
                           ['seq3@0002', 'asv03', 0, 15, 'Bludhaven'],
                           ['seq5', 'asv04', 0, 15, 'Bludhaven'],
                           ['seq5', 'asv05', 1, 15, 'Bludhaven'],
                           ['seq6', 'asv04', 1, 15, 'Bludhaven'],
                           ['seq6', 'asv05', 0, 15, 'Bludhaven']], 
                           dtype=object),
            columns=['kmer', 'asv', 'mismatch', 'length', 'region'],
            )
        known['length'] = known['length'].astype(int)
        known['mismatch'] = known['mismatch'].astype(int)
        known.set_index('kmer', inplace=True)
        
        match, discard = align_regional_kmers(kmers, rep_set, 
                                              region='Bludhaven', 
                                              debug=True)
        pdt.assert_frame_equal(known, match)
        pdt.assert_series_equal(
            discard.view(pd.Series).astype(str), 
            pd.Series({'asv06': 'AGAGTTTCTGAATCC'}), 
        )


if __name__ == '__main__':
    main()
