from unittest import TestCase, main

import os
import shutil
import warnings


import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact, Metadata
from qiime2.plugins.sidle import methods as sidle
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
        self.manfiest = Metadata(pd.DataFrame(
            data=[[os.path.join(self.base_dir, 'region1_db_map.qza'),
                   os.path.join(self.base_dir, 'region1_align.qza'),
                   os.path.join(self.base_dir, 'region1_counts.qza')],
                  [os.path.join(self.base_dir, 'region1_db_map.qza'),
                   os.path.join(self.base_dir, 'region1_align.qza'),
                   os.path.join(self.base_dir, 'region1_counts.qza')]],
            columns=['kmer-map', 'alignment-map', 'frequency-table'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
            ))
        self.seq_map = ts.seq_map
        self.taxonomy = ts.taxonomy
        self.count1 = ts.region1_counts
    
    def test_plugin_setup(self):
        self.assertEqual(plugin.name, 'sidle')

    def test_filter_degenerate_sequences(self):
        known = self.ref_seqs.view(pd.Series).copy().astype(str).drop(['seq3'])

        test = sidle.filter_degenerate_sequences(
            self.ref_seqs, 
            max_degen=0, 
            debug=True
            )
        test = test.filtered_sequences.view(pd.Series).astype(str)
        pdt.assert_series_equal(known, test)

    def test_extract_regional_database(self):
        test_seqs, test_map = \
            sidle.extract_regional_database(self.ref_seqs,
                                            fwd_primer='WANTCAT',
                                            rev_primer='CATCATCAT',
                                            region='Bludhaven',
                                            trim_length=15,
                                            primer_mismatch=1,
                                            debug=True,
                                            )
        pdt.assert_series_equal(
            self.region1_db_seqs.view(pd.Series).astype(str),
            test_seqs.view(pd.Series).astype(str)
            )
        pdt.assert_frame_equal(self.region1_db_map.view(pd.DataFrame),
                               test_map.view(pd.DataFrame))

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
        test_align, test_discard = \
            sidle.align_regional_kmers(self.region1_db_seqs,
                                       self.rep_seqs1,
                                       region='Bludhaven',
                                       max_mismatch=2,
                                       debug=True,
                                       )
        self.assertEqual(len(test_discard.view(pd.Series)), 0)
        pdt.assert_frame_equal(self.align1.view(pd.DataFrame),
                               test_align.view(pd.DataFrame))

    def test_reconstruct_counts(self):

        known_summary = pd.DataFrame.from_dict(orient='index', data={
            'seq1': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1.,
                     'stdv-kmer-per-region': 0.,
                     'mapped-asvs': 'asv01|asv06'
                    },
            'seq2': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'mapped-asvs': 'asv01|asv07',
                    },
            'seq3': {'num-regions': 2, 
                     'total-kmers-mapped': 3, 
                     'mean-kmer-per-region': 1.5,
                     'stdv-kmer-per-region': np.std([1, 2], ddof=1),
                     'mapped-asvs': 'asv02|asv03|asv08'
                    },
            'seq4': {'num-regions': 1, 
                     'total-kmers-mapped': 1, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'mapped-asvs': 'asv09'
                    },
            'seq5': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'mapped-asvs': 'asv04|asv05|asv10',
                    },
            'seq6': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'mapped-asvs': 'asv04|asv05|asv11',
                    },
            })
        known_summary.index.set_names('feature-id', inplace=True)        
        manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(self.base_dir, 'region1_db_map.qza'),
                   os.path.join(self.base_dir, 'region1_align.qza'),
                   os.path.join(self.base_dir, 'region1_counts.qza'), '0'],
                  [os.path.join(self.base_dir, 'region2_db_map.qza'),
                   os.path.join(self.base_dir, 'region2_align.qza'),
                   os.path.join(self.base_dir, 'region2_counts.qza'), '1']],
            columns=['kmer-map', 'alignment-map', 'frequency-table', 'region-order'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
            ))
        count_table, summary, mapping = \
            sidle.reconstruct_counts(manifest, debug=True)
        pdt.assert_frame_equal(
            count_table.view(pd.DataFrame),
            pd.DataFrame( 
                data=np.array([[100.,  50,   0,  50,  50, 50],
                               [100.,  25, 100,  25,  25, 25],
                               [  0., 100, 100,   0,  50, 50]]),
                index=pd.Index(['sample1', 'sample2', 'sample3'], 'sample-id'),
                columns=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6']
            )
        )
        pdt.assert_series_equal(self.seq_map.view(pd.Series), mapping.view(pd.Series))
        pdt.assert_frame_equal(known_summary, summary.view(pd.DataFrame))

    def test_reconstruct_taxonomy(self):
        test = sidle.reconstruct_taxonomy(self.seq_map, 
                                          self.taxonomy,
                                          database='greengenes',
                                          define_missing='ignore'
                                          ).reconstructed_taxonomy
        pdt.assert_series_equal(self.taxonomy.view(pd.Series),
                                test.view(pd.Series))

    def test_integration(self):
        # This will run through a slightly more complex dataset...
        base_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/integration')
        test_dir = os.path.join(base_dir, 'test')
        known_dir = os.path.join(base_dir, 'known')
        data_dir = os.path.join(base_dir, 'data')
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.makedirs(test_dir)

        # Sets up the temporary manifest
        manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(test_dir, 'region1-kmer-db.qza'),
                   os.path.join(test_dir, 'region1-align-map.qza'),
                   os.path.join(data_dir, 'region1-counts.qza'), '0'],
                  [os.path.join(test_dir, 'region2-kmer-db.qza'),
                   os.path.join(test_dir, 'region2-align-map.qza'),
                   os.path.join(data_dir, 'region2-counts.qza'), '1'],
                  [os.path.join(test_dir, 'region3-kmer-db.qza'),
                   os.path.join(test_dir, 'region3-align-map.qza'),
                   os.path.join(data_dir, 'region3-counts.qza'), '2']],
            columns=['kmer-map', 'alignment-map', 'frequency-table', 'region-order'],
            index=pd.Index(['1', '2', '3'], name='id'),
        ))    

        ### Sequence extraction
        region1_seqs, region1_map = sidle.extract_regional_database(
            Artifact.load(os.path.join(data_dir, 'database_sequences.qza')),
            fwd_primer='TGGCGGACGGGTGAGTAA',
            rev_primer='CTGCTGCCTCCCGTAGGA',
            trim_length=50,
            primer_mismatch=1,
            region='1',
            debug=True,
            )
        known = \
            Artifact.load(os.path.join(known_dir, 'region1-kmer-seqs.qza'))
        pdt.assert_series_equal(region1_seqs.view(pd.Series).astype(str),
                                known.view(pd.Series).astype(str))
        known = \
            Artifact.load(os.path.join(known_dir, 'region1-kmer-map.qza'))
        pdt.assert_frame_equal(known.view(pd.DataFrame), 
                              region1_map.view(pd.DataFrame))

        region2_seqs, region2_map = sidle.extract_regional_database(
            Artifact.load(os.path.join(data_dir, 'database_sequences.qza')),
            fwd_primer='CAGCAGCCGCGGTAATAC',
            rev_primer='CGCATTTCACCGCTACAC',
            trim_length=50,
            primer_mismatch=1,
            region='2',
            debug=True,
            )
        known = \
            Artifact.load(os.path.join(known_dir, 'region2-kmer-seqs.qza'))
        pdt.assert_series_equal(region2_seqs.view(pd.Series).astype(str),
                                known.view(pd.Series).astype(str))
        known = \
            Artifact.load(os.path.join(known_dir, 'region2-kmer-map.qza'))
        pdt.assert_frame_equal(known.view(pd.DataFrame), 
                              region2_map.view(pd.DataFrame))
        region3_seqs, region3_map = sidle.extract_regional_database(
            Artifact.load(os.path.join(data_dir, 'database_sequences.qza')),
            fwd_primer='GCACAAGCGGTGGAGCAT',
            rev_primer='CGCTCGTTGCGGGACTTA',
            trim_length=50,
            primer_mismatch=1,
            region='3',
            debug=True,
            )
        known = \
            Artifact.load(os.path.join(known_dir, 'region3-kmer-seqs.qza'))
        pdt.assert_series_equal(region3_seqs.view(pd.Series).astype(str),
                                known.view(pd.Series).astype(str))
        known = \
            Artifact.load(os.path.join(known_dir, 'region3-kmer-map.qza'))
        pdt.assert_frame_equal(known.view(pd.DataFrame), 
                              region3_map.view(pd.DataFrame))

        region1_map.save(os.path.join(test_dir, 'region1-kmer-db.qza'))
        region2_map.save(os.path.join(test_dir, 'region2-kmer-db.qza'))
        region3_map.save(os.path.join(test_dir, 'region3-kmer-db.qza'))

        
        ### Regiomal Alignment
        align1, discard1 = sidle.align_regional_kmers(
            region1_seqs, 
            Artifact.load(os.path.join(data_dir, 'region1-rep-seq.qza')),
            region='1',
            max_mismatch=2,
            debug=True,
            )
        known = \
            Artifact.load(os.path.join(known_dir, 'region1-align-map.qza'))
        pdt.assert_frame_equal(align1.view(pd.DataFrame), 
                               known.view(pd.DataFrame))

        align2, discard2 = sidle.align_regional_kmers(
            region2_seqs, 
            Artifact.load(os.path.join(data_dir, 'region2-rep-seq.qza')),
            region='2',
            max_mismatch=2,
            debug=True,
            )
        known = \
            Artifact.load(os.path.join(known_dir, 'region2-align-map.qza'))
        pdt.assert_frame_equal(align2.view(pd.DataFrame), 
                               known.view(pd.DataFrame))
        
        align3, discard3 = sidle.align_regional_kmers(
            region3_seqs, 
            Artifact.load(os.path.join(data_dir, 'region3-rep-seq.qza')),
            region='3',
            max_mismatch=2,
             debug=True,
            )
        known = \
            Artifact.load(os.path.join(known_dir, 'region3-align-map.qza'))
        pdt.assert_frame_equal(align3.view(pd.DataFrame), 
                               known.view(pd.DataFrame))

        align1.save(os.path.join(test_dir, 'region1-align-map.qza'))
        align2.save(os.path.join(test_dir, 'region2-align-map.qza'))
        align3.save(os.path.join(test_dir, 'region3-align-map.qza'))


        ### Reconstruction
        table, summary, map_ = sidle.reconstruct_counts(
            manifest,
            debug=True,
            count_degenerates=False,
        )
        known = \
            Artifact.load(os.path.join(known_dir, 'reconstructed-table.qza'))
        pdt.assert_frame_equal(known.view(pd.DataFrame), 
                               table.view(pd.DataFrame))
        known = \
            Artifact.load(os.path.join(known_dir, 'reconstructed-summary.qza'))
        # ASV mapping was optional in the  original sidle. This is  tested
        # elsewhere  and dealing w ith it is going to suck. 
        pdt.assert_frame_equal(
            known.view(pd.DataFrame),
            summary.view(pd.DataFrame).drop(columns=['mapped-asvs'])
            )
        known = \
            Artifact.load(os.path.join(known_dir, 'sidle-reconstruction.qza'))
        pdt.assert_frame_equal(
            known.view(pd.DataFrame),
            map_.view(pd.DataFrame)
            )

        shutil.rmtree(test_dir)

if __name__ == '__main__':
    main()