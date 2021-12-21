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

from qiime2 import Artifact, Metadata
from qiime2.plugins.sidle import actions as sidle
from q2_sidle.plugin_setup import plugin

import q2_sidle.tests.test_set as ts


class PluginSetupTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        self.ref_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'full_db.qza'))
        self.region1_db_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'region1_db_seqs.qza'))
        self.region2_db_seqs = \
            Artifact.load(os.path.join(self.base_dir, 'region2_db_seqs.qza'))
        self.region1_db_map = ts.region1_db_map
        self.region2_db_map = ts.region2_db_map
        self.rep_seqs1 = \
            Artifact.load(os.path.join(self.base_dir, 'region1_rep_set.qza'))

        self.align1 = ts.region1_align
        self.align2 = ts.region2_align
        self.kmer_map1 = ts.region1_db_map
        self.kmer_map2 = ts.region2_db_map
        self.table1 = ts.region1_counts
        self.table2 = ts.region2_counts
        self.seq_map = ts.seq_map
        self.taxonomy = ts.taxonomy
        self.count1 = ts.region1_counts
        self.database_summary = ts.db_summary
    
    def test_plugin_setup(self):
        self.assertEqual(plugin.name, 'sidle')

    def test_reverse_complement_sequence(self):
        known = pd.Series(
            data=[DNA('ATGATGATG', metadata={'id': 'seq1'})],
            index=pd.Index(['seq1'], name=None),
            )
        input_ = Artifact.import_data(
            'FeatureData[Sequence]', 
            pd.Series(data=[DNA('CATCATCAT', metadata={'id': 'seq1'})],
                      index=pd.Index(['seq1'], name=None),
                      )
            )
        test = sidle.reverse_complement_sequence(input_).reverse_complement
        self.assertEqual(str(test.type), 'FeatureData[Sequence]')
        pdt.assert_series_equal(known.astype(str),
                                test.view(pd.Series).astype(str))

    def test_reverse_complement_aligned_sequence(self):
        known = pd.Series(
            data=[DNA('--ATGATGATG', metadata={'id': 'seq1'}),
                  DNA('TTATGATGA--', metadata={'id': 'seq2'})],
            index=pd.Index(['seq1', 'seq2'], name=None),
            )
        input_ = Artifact.import_data(
            'FeatureData[AlignedSequence]', 
            pd.Series(data=[DNA('CATCATCAT--', metadata={'id': 'seq1'}),
                            DNA('--TCATCATAA', metadata={'id': 'seq2'})],
                      index=pd.Index(['seq1', 'seq2'], name='feature-id'),
                      ),
            pd.Series,
            )
        test = sidle.reverse_complement_sequence(input_).reverse_complement
        self.assertEqual(str(test.type), 'FeatureData[AlignedSequence]')
        pdt.assert_series_equal(known.astype(str), 
                                test.view(pd.Series).astype(str))

    def test_find_first_alignment_position(self):
        known = pd.DataFrame(
            data=np.vstack([np.array([12.] * 5, dtype=float),
                            np.array([101.] * 5, dtype=float),
                            np.array(['fwd'] * 5, dtype=object),
                ]).T,
            columns=['starting-position', 'sequence-counts', 'direction'],
            index=pd.Index(['asv01', 'asv02', 'asv03', 'asv04', 'asv05'],
                           name='feature-id'),
            )
        known['starting-position'] = known['starting-position'].astype(int).astype(str)
        known['sequence-counts'] = known['sequence-counts'].astype(float)
        test = sidle.find_first_alignment_position(
            alignment=ts.extra_alignment,
            representative_sequences=self.rep_seqs1,
            ).position_summary
        pdt.assert_frame_equal(known, test.view(Metadata).to_dataframe())

    def test_prepare_extracted_region(self):
        test_seqs, test_map = \
            sidle.prepare_extracted_region(self.region2_db_seqs,
                                           region='Gotham',
                                           trim_length=15,
                                           fwd_primer='CACCTCGTN',
                                           rev_primer='CACGTCAK',
                                           debug=True,
                                           )
        pdt.assert_series_equal(
            test_seqs.view(pd.Series).astype(str),
            self.region2_db_seqs.view(pd.Series).astype(str)
            )
        pdt.assert_frame_equal(test_map.view(pd.DataFrame), 
                               self.region2_db_map.view(pd.DataFrame))

    def test_trim_dada2_posthoc(self):
        test_table, test_seqs = \
            sidle.trim_dada2_posthoc(self.count1,
                                     self.rep_seqs1,
                                     hashed_feature_ids=False,
                                     )

        known_seqs = self.rep_seqs1.view(pd.Series).astype(str).copy()
        known_seqs.index = known_seqs
        known_table = self.count1.view(pd.DataFrame).copy()
        known_table.columns = known_seqs.index

        pdt.assert_frame_equal(test_table.view(pd.DataFrame), known_table)
        pdt.assert_series_equal(test_seqs.view(pd.Series).astype(str), 
                                known_seqs)

    def test_align_regional_kmers(self):
        warnings.filterwarnings('ignore', 
                                category=skbio.io.FormatIdentificationWarning)
        test_align = \
            sidle.align_regional_kmers(self.region1_db_seqs,
                                       self.rep_seqs1,
                                       region='Bludhaven',
                                       max_mismatch=2,
                                       debug=True,
                                       ).regional_alignment
        pdt.assert_frame_equal(
            self.align1.view(pd.DataFrame),
            test_align.view(pd.DataFrame).sort_values(['kmer', 'asv'])
            )

    def test_reconstruct_database(self):
        mapping, summary = sidle.reconstruct_database(
            region=['Bludhaven', 'Gotham'],
            kmer_map=[self.kmer_map1, self.kmer_map2],
            regional_alignment=[self.align1, self.align2],
            debug=True,
            )
        pdt.assert_frame_equal(mapping.view(pd.DataFrame), 
                               self.seq_map.view(pd.DataFrame))
        pdt.assert_frame_equal(summary.view(pd.DataFrame),
                               self.database_summary.view(pd.DataFrame))

    def test_reconstruct_counts(self):
        count_table = \
            sidle.reconstruct_counts(
                region=['Bludhaven', 'Gotham'],
                regional_alignment=[self.align1, self.align2],
                regional_table=[self.table1, self.table2],
                database_map=self.seq_map,
                database_summary=self.database_summary,
                debug=True,
                min_abund=1e-2,
                min_counts=10
                ).reconstructed_table
        pdt.assert_frame_equal(
            count_table.view(pd.DataFrame),
            pd.DataFrame( 
                data=np.array([[100.,  50,   0,  50,  50, 50],
                               [100.,  25, 100,  25,  25, 25],
                               [  0., 100, 100,   0,  50, 50]]),
                index=pd.Index(['sample1', 'sample2', 'sample3']),
                columns=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6']
            )
        )

    def test_reconstruct_taxonomy(self):
        test = sidle.reconstruct_taxonomy(self.seq_map, 
                                          self.taxonomy,
                                          database='greengenes',
                                          define_missing='ignore'
                                          ).reconstructed_taxonomy
        pdt.assert_series_equal(self.taxonomy.view(pd.Series),
                                test.view(pd.Series))

    def test_reconstruct_fragment_rep_seqs(self):
        recon_map = Artifact.import_data(
            'FeatureData[SidleReconstruction]', 
            pd.DataFrame(
                data=np.array([['seq01|seq02', 0,  'WANTCAT', 0, 'WANTCAT', 15], 
                               ['seq01|seq02', 0, 'WANTCAT', 0, 'WANTCAT', 15], 
                               ['seq03|seq04', 0, 'WANTCAT', 1, 'CACCTCGTN', 15], 
                               ['seq03|seq04', 0, 'CACCTCGTN', 1, 'CACCTCGTN', 15], 
                               ['seq05', 0, 'WANTCAT', 1, 'CACCTCGTN', 15],
                               ],  dtype=object),
                index=pd.Index(['seq01', 'seq02', 'seq03', 'seq04', 'seq05'], 
                                name='db-seq'),
                columns=['clean_name', 'first-region', 'first-fwd-primer',  
                         'last-region', 'last-fwd-primer', 'last-kmer-length'],
            )
            )
        recon_summary = Artifact.import_data(
            'FeatureData[ReconstructionSummary]',
            Metadata(pd.DataFrame(data=[[1, 2, 2, 0, 'asv01|asv02'],
                                        [2, 3, 1.5, np.std([1, 2], ddof=1), 
                                         'asv03|asv04'],
                                        [2, 2, 1, 0, 'asv07|asv08']],
                                 index=pd.Index(['seq01|seq02', 'seq03|seq04', 
                                                 'seq05'], name='feature-id'),
                                columns=['num-regions', 'total-kmers-mapped', 
                                         'mean-kmer-per-region', 
                                         'stdv-kmer-per-region', 
                                         'mapped-asvs']))
        )
        aligned_seqs = Artifact.import_data(
            'FeatureData[AlignedSequence]', 
            skbio.TabularMSA([
                DNA('CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC-------------------'
                    '--------------', metadata={'id': 'seq01'}),
                DNA('CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC-------------------'
                    '--------------', metadata={'id': 'seq02'}),
                DNA('CATAGTCATWTCCGCGTTGGAGTTATGATGATGAWACCACCTCGTCCCAGTTCCGC'
                    'GCTTCTGACGTGC-', metadata={'id': 'seq03'}),
                DNA('------------------GGAGTTATGATGA--AGACCACCTCGTCCCAGTTCCGC'
                    'GCTTCTGACGTGCC', metadata={'id': 'seq04'}),
                DNA('CATAGTCATCGTTTATGTATGCCCATGATGATGCGAGCACCTCGTATGGATGTAGA'
                    'GCCACTGACGTGCG', metadata={'id': 'seq05'}),
            ])
        )
        known = pd.Series(
            data=['GCGAAGCGGCTCAGG',
                  'WTCCGCGTTGGAGTTATGATGATGAGACCACCTCGTCCCAGTTCCGCGCTTC'],
            index=pd.Index(['seq01|seq02', 'seq03|seq04']),
            )
        test = sidle.reconstruct_fragment_rep_seqs(
            reconstruction_map=recon_map, 
            reconstruction_summary=recon_summary, 
            aligned_sequences=aligned_seqs,
            ).representative_fragments
        pdt.assert_series_equal(known, test.view(pd.Series).astype(str))


if __name__ == '__main__':
    main()