from unittest import TestCase, main

from copy import copy
import os
import shutil
import warnings


import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio
from skbio import DNA

from qiime2 import Artifact, Metadata, Visualization
from qiime2.plugins.sidle import actions as sidle
from q2_sidle.plugin_setup import plugin

import q2_sidle.tests.test_set as ts


class PipelineTest(TestCase):
    def setUp(self):
        self.align1 = ts.region1_align
        self.align2 = ts.region2_align
        self.kmer_map1 = ts.region1_db_map
        self.kmer_map2 = ts.region2_db_map
        self.table1 = ts.region1_counts
        self.table2 = ts.region2_counts
        self.taxonomy = ts.taxonomy
        self.seq_map = ts.seq_map
        self.database_summary = ts.db_summary
        self.extra_alignment = ts.extra_alignment

    def test_sidle_reconstruction(self):
        mapping, summary, counts, taxonomy  = sidle.sidle_reconstruction(
            region=['Bludhaven', 'Gotham'],
            kmer_map=[self.kmer_map1, self.kmer_map2],
            regional_alignment=[self.align1, self.align2],
            regional_table=[self.table1, self.table2],
            reference_taxonomy=self.taxonomy,
            database='greengenes',
            define_missing='ignore',
            min_counts=10,
            debug=True
            )
        pdt.assert_frame_equal(mapping.view(pd.DataFrame), 
                               self.seq_map.view(pd.DataFrame))
        pdt.assert_frame_equal(summary.view(pd.DataFrame),
                               self.database_summary.view(pd.DataFrame))
        pdt.assert_frame_equal(
            counts.view(pd.DataFrame),
            pd.DataFrame( 
                data=np.array([[100.,  50,   0,  50,  50, 50],
                               [100.,  25, 100,  25,  25, 25],
                               [  0., 100, 100,   0,  50, 50]]),
                index=pd.Index(['sample1', 'sample2', 'sample3']),
                columns=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6']
            )
        )
        pdt.assert_series_equal(self.taxonomy.view(pd.Series),
                                taxonomy.view(pd.Series))

    def test_map_alignment_positions(self):
        ori_alignment = Artifact.import_data('FeatureData[AlignedSequence]', pd.Series(
            data=[
                skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTAAGAGAGGCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq1'}),
                skbio.DNA('-----AAATCATGCGAAGCGGCTCAGGATGATGATGGGTGAGTCACCTCGTCAGAGTTTCTGAATCCATGACGTG---ACCAGC', metadata={'id': 'seq2'}),
                skbio.DNA('TATGGTACTCATWTCCGCGTTGGAGTTATGATGATGGGGTGA-CACCTCGTTCCAGTTCCGCGCTTCATGACGTGCTGACC---', metadata={'id': 'seq3'}),
                skbio.DNA('------------------------------AAGGCGGGTGAG-CACCTCGTCCCGGAGACGAGAGGCATGACGTG---ATCCGT', metadata={'id': 'seq4'}),
                skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGTGGATGTAGAGCCACCTGACGTGC--ACCTG-', metadata={'id': 'seq5'}),
                skbio.DNA('-AGGCTAGTCATCGTTTATGTATGCCCATGATGATGGGGTGAGCACCTCGTGAAAATGTAGAGCCACCTGACGTGC--ACC---', metadata={'id': 'seq6'}),
                ],
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
                            dtype='object')
            )
        )
        rep_seqs = Artifact.import_data('FeatureData[Sequence]', pd.Series({
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
        )
        known = pd.DataFrame(
            data=np.vstack([
                np.hstack([np.array([12] * 5), np.array([52] * 3), 0, 28]),
                np.array([101.] * 10),
                ]).T,
            columns=['starting-position', 'sequence-counts'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05', 
                            'asv06', 'asv07', 'asv08', 'asv09', 'asv10'],
                           name='feature-id'),
            )
        known['starting-position'] = \
            known['starting-position'].astype(int).astype(str)
        known['direction'] = 'fwd'
        known = Metadata(known)
        expanded, starts, viz_ = \
        sidle.map_alignment_positions(
            alignment=ori_alignment,
            sequences=rep_seqs,
            direction='fwd',
            add_fragments=True, # Improves testing behavior
            )
        self.assertEqual(str(expanded.type), 'FeatureData[AlignedSequence]')
        self.assertEqual(str(starts.type), 'FeatureData[AlignmentPosSummary]')
        self.assertTrue(isinstance(viz_, Visualization))

        pdt.assert_series_equal(self.extra_alignment.view(pd.Series).astype(str),
                               expanded.view(pd.Series).astype(str))
        pdt.assert_frame_equal(known.to_dataframe(), 
                               starts.view(Metadata).to_dataframe())




if __name__ == '__main__':
    main()
