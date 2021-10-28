from unittest import TestCase, main

import copy
import os
import warnings

import biom
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio

from qiime2 import Artifact, Metadata

from q2_sidle._primerless import (reverse_complement_sequence,
                                  find_first_alignment_position,
                                  _generate_align_mask,
                                  _get_first_align_pos,
                                  )
import q2_sidle.tests.test_set as ts

class PrimerlessTest(TestCase):
    def setUp(self):
        self.expanded_alignment = pd.Series(
            data=[
                skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTAAGAGAGGCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq1'}),
                skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTCAGAGTTTCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq2'}),
                skbio.DNA('TATGGTACTCATWTCCGCGTTGGAGTTATGATGATGGGGTGA-CACCTCGTTCCAGTTCCGCGCTTCATGACGTGCTGACC---', metadata={'id': 'seq3'}),
                skbio.DNA('------------------------------AAGGCGGGTGAG-CACCTCGTCCCGGAGACGAGAGGCATGACGTG---ATCCGT', metadata={'id': 'seq4'}),
                skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGTGGATGTAGAGCCACCTGACGTGC--ACCTG-', metadata={'id': 'seq5'}),
                skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGAAAATGTAGAGCCACCTGACGTGC--ACC---', metadata={'id': 'seq6'}),
                skbio.DNA('------------GCGAAGCGGCTCAGG---------------------------------------------------------', metadata={'id': 'asv01'}),
                skbio.DNA('------------ATCCGCGTTGGAGTT---------------------------------------------------------', metadata={'id': 'asv02'}),
                skbio.DNA('------------TTCCGCGTTGGAGTT---------------------------------------------------------', metadata={'id': 'asv03'}),
                skbio.DNA('------------CGTTTATGTATGCCC---------------------------------------------------------', metadata={'id': 'asv04'}),
                skbio.DNA('------------CGTTTATGTATGCCT---------------------------------------------------------', metadata={'id': 'asv05'}),
                skbio.DNA('----------------------------------------------------AGAGAGGCTGAATCC-----------------', metadata={'id': 'asv06'}),
                skbio.DNA('----------------------------------------------------AGAGTTTCTGAATCC-----------------', metadata={'id': 'asv07'}),
                skbio.DNA('----------------------------------------------------CCGGAGACGAGAGGC-----------------', metadata={'id': 'asv08'}),
                skbio.DNA('AAGGCTAGTCATCGT---------------------------------------------------------------------', metadata={'id': 'asv09'}),
                skbio.DNA('----------------------------TGAGTWWGGGAGGGA-----------------------------------------', metadata={'id': 'asv10'}),
                ],
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6', 'asv01', 'asv02',
                            'asv03', 'asv04', 'asv05', 'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                            dtype='object')
            )
        self.rep_seqs = pd.Series({
            'asv01': skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'asv01'}),
            'asv02': skbio.DNA('ATCCGCGTTGGAGTT', metadata={'id': 'asv02'}),
            'asv03': skbio.DNA('TTCCGCGTTGGAGTT', metadata={'id': 'asv03'}),
            'asv04': skbio.DNA('CGTTTATGTATGCCC', metadata={'id': 'asv04'}),
            'asv05': skbio.DNA('CGTTTATGTATGCCT', metadata={'id': 'asv05'}),
            'asv06': skbio.DNA('AGAGAGGCTGAATCC', metadata={'id': 'asv06'}),
            'asv07': skbio.DNA('AGAGTTTCTGAATCC', metadata={'id': 'asv07'}),
            'asv08': skbio.DNA('CCGGAGACGAGAGGC', metadata={'id': 'asv08'}),
            'asv09': skbio.DNA('AAGGCTAGTCATCGT', metadata={'id': 'asv09'}),
            'asv10': skbio.DNA('TGAGTWWGGGAGGGA', metadata={'id': 'asv10'}),
            })
        self.table = biom.Table(
            data=np.array([[100,  10,  90,  50,  50,  10,   0,  50,  25,   5],
                           [200,  90,  10,  50,  50,  10,  10,  50, 500,   8],
                           ], dtype=float).T,
            observation_ids=('asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                             'asv06', 'asv07', 'asv08', 'asv09', 'asv10'),
            sample_ids=('grayson', 'todd')
            )
        self.coverage = pd.DataFrame(
            data=np.array([
                list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'), 
                list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'),
                list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'),
                list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'),
                list('000000000000111111111111111000000000000000000000000000000000000000000000000000000000'), 
                list('000000000000000000000000000000000000000000000000000011111111111111100000000000000000'),
                list('000000000000000000000000000000000000000000000000000011111111111111100000000000000000'),
                list('000000000000000000000000000000000000000000000000000011111111111111100000000000000000'),
                list('111111111111111000000000000000000000000000000000000000000000000000000000000000000000'), 
                list('000000000000000000000000000011111111111111100000000000000000000000000000000000000000'),
            ]).astype(int),
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'])
            )

    def test_reverse_complement_sequence(self):
        input_ = pd.Series({
                'grayson': skbio.DNA('AAAA', metadata={'id': 'grayson'}),
                'todd': skbio.DNA('ATAT', metadata={'id': 'todd'}),
                'drake': skbio.DNA('GCGC', metadata={'id': 'drake'}),
                'brown': skbio.DNA('GGGG', metadata={'id': 'brown'}),
            })
        test = reverse_complement_sequence(input_)

        known = pd.Series({
            'grayson': skbio.DNA('TTTT', metadata={'id': 'grayson'}),
            'todd': skbio.DNA('ATAT', metadata={'id': 'todd'}),
            'drake': skbio.DNA('GCGC', metadata={'id': 'drake'}),
            'brown': skbio.DNA('CCCC', metadata={'id': 'brown'}),
            })
        pdt.assert_series_equal(test.astype(str), known.astype(str))

    def test_find_find_first_alignment_position(self):
        known = pd.DataFrame(
            data=np.vstack([
                np.hstack([np.array([12.] * 5), np.array([52.] * 3), 0., 28.]),
                np.array([101.] * 10),
                ]).T,
            columns=['starting-position', 'sequence-counts'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                           name='feature-id'),
            )
        known['direction'] = 'fwd'
        known['sequence-counts'] = known['sequence-counts'].astype(float)
        known['starting-position'] = \
            known['starting-position'].astype(int).astype(str)

        test = find_first_alignment_position(
            alignment=self.expanded_alignment,
            representative_sequences=self.rep_seqs,
            ).to_dataframe()

        pdt.assert_frame_equal(known, test)

    def test_find_first_alignment_position_error(self):
        test_seqs = pd.Series({
            'nightwing': skbio.DNA("AAAAA", metadata={'id': 'nightwing'}),
            'arsenel':  skbio.DNA("CCCCC", metadata={'id': 'arsenel'}),
            'troya':  skbio.DNA("GGGGG", metadata={'id': 'troya'}),
            })
        with self.assertRaises(ValueError) as err:
            find_first_alignment_position(
                alignment=self.expanded_alignment,
                representative_sequences=test_seqs,
            )
        self.assertEqual(str(err.exception),
                         ('All the representative sequences must be present '
                         'in the alignment. Please make sure that you have '
                         'the correct representative sequences and '
                         'alignment.')    
                        )

    def test_generate_align_mask(self):
        test = _generate_align_mask(self.expanded_alignment, 
                                    self.rep_seqs.index)
        pdt.assert_frame_equal(self.coverage, test)

    def test_get_first_align_pos_no_table(self):
        known = pd.DataFrame(
            data=np.vstack([
                np.hstack([np.array([12.] * 5), np.array([52.] * 3), 0., 28.]),
                np.array([101.] * 10),
                ]).T,
            columns=['starting-position', 'sequence-counts'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                           name='feature-id'),
            )
        test = _get_first_align_pos(self.coverage)
        pdt.assert_frame_equal(known.astype(float), test.astype(float))

    def test_get_first_align_pos_table(self):
        known = pd.DataFrame(
            data=np.array([
                [ 12,  12,  12,  12,  12,  52,  52,  52,   0,  28],
                [300, 100, 100, 100, 100,  20,  10, 100, 525,  13]
                ]).T,
            columns=['starting-position', 'sequence-counts'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                           name='feature-id'),
            )
        test = _get_first_align_pos(self.coverage, self.table)
        pdt.assert_frame_equal(known.astype(float), test.astype(float))


if __name__ == '__main__':
    main()