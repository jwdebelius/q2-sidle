from unittest import TestCase, main

import os

import dask
import dask.dataframe as dd
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact, Metadata
from qiime2.plugin import ValidationError

from q2_sidle._reconstruct import (reconstruct_counts,
                                   _construct_align_mat,
                                   _count_mapping,
                                   _detangle_names,
                                   _expand_duplicate_sequences,
                                   _get_clean,
                                   _get_shared_seqs,
                                   _get_unique_kmers,
                                   _scale_relative_abundance,
                                   _solve_ml_em_iterative_1_sample,
                                   _solve_iterative_noisy,
                                   _sort_untidy,
                                   _tidy_sequence_set,

                                   _untangle_database_ids,
                                   
                                   )
import q2_sidle.tests.test_set as ts


class ReconstructTest(TestCase):
    def setUp(self):
        self.tangle = pd.DataFrame(
            data=[['seq01', 0, 'seq01', 'seq02'],
                  ['seq02', 0, 'seq01', 'seq02'],
                  ['seq01', 1, 'seq01', np.nan],
                  ['seq02', 1, 'seq02', np.nan],
                  ['seq05', 0, 'seq05', np.nan],
                  ['seq06', 0, 'seq06', np.nan],
                  ['seq05', 1, 'seq05', np.nan],
                  ['seq05', 1, 'seq05', 'seq06'],
                  ['seq06', 1, 'seq05', 'seq06']
                  ],
            columns=['db-seq', 'region', '0', '1'],
            )
        self.id_set = [np.array(['seq05', 'seq06'])]
        self.pair = pd.DataFrame(
            data=[['seq05', 0, 'seq05', 'seq05'],
                  ['seq05@001', 1, 'seq05', 'seq05'],
                  ['seq05@002', 1, 'seq05', 'seq05'],
                  ['seq06', 0, 'seq06', 'seq06'],
                  ['seq06', 1, 'seq06', 'seq06']
                  ],
            columns=['kmer', 'region', 'db-seq', 'clean_name'],
            )
        self.match1 = ts.region1_align.view(pd.DataFrame).copy()
        self.match2 = ts.region2_align.view(pd.DataFrame).copy()
        self.kmer1 = ts.region1_db_map.view(pd.DataFrame).copy()

        self.align1 = pd.DataFrame(
            data=np.array([['seq1', 'asv01', 'Bludhaven', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq2', 'asv01', 'Bludhaven', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq3', 'asv02', 'Bludhaven', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq5', 'asv04', 'Bludhaven', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq5', 'asv05', 'Bludhaven', 1, 15, 0.00155, 0.92757, 0.0, 0.92757, False, 0.0008],
                           ['seq6', 'asv04', 'Bludhaven', 1, 15, 0.00155, 0.92757, 0.0, 0.92757, False, 0.0008],
                           ['seq6', 'asv05', 'Bludhaven', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638]], dtype=object),
            columns=['clean_name', 'asv', 'region', 'mismatch', 'length', 
                     'match_prob', 'perfect_prob',  
                     'max_error_prob', 'error_thresh', 'perf_match', 'norm']
            )
        self.align2 = pd.DataFrame(
            data=np.array([['seq1', 'asv06', 'Gotham', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq2', 'asv07', 'Gotham', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq3', 'asv08', 'Gotham', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq4', 'asv09', 'Gotham', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.9276],
                           ['seq5', 'asv10', 'Gotham', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq6', 'asv11', 'Gotham', 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638]], dtype=object),
            columns=['clean_name', 'asv', 'region', 'mismatch', 'length', 
                     'match_prob', 'perfect_prob',  
                     'max_error_prob', 'error_thresh', 'perf_match', 'norm']
            )
        self.int_cols = ['length', 'mismatch']
        self.float_cols = ['match_prob', 'perfect_prob', 'max_error_prob', 
                           'error_thresh', 'norm']
        self.align1[self.int_cols] = self.align1[self.int_cols].astype(int)
        self.align2[self.int_cols] = self.align2[self.int_cols].astype(int)
        self.align1[self.float_cols] = self.align1[self.float_cols].astype(float)
        self.align2[self.float_cols] = self.align2[self.float_cols].astype(float)
        self.align1['perf_match'] = self.align1['perf_match'] == True
        self.align2['perf_match'] = self.align2['perf_match'] == True
        self.align1.sort_values(['clean_name', 'asv'], inplace=True)
        self.align2.sort_values(['clean_name', 'asv'], inplace=True)
        self.align1.reset_index(drop=True, inplace=True)
        self.align2.reset_index(drop=True, inplace=True)
        self.align1 = self.align1.round(4)
        self.align2 = self.align2.round(4)

        self.seq_summary = pd.DataFrame.from_dict(orient='index', data={
            'seq1': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batgirl;Gordon;Barbara',
                     'mapped-asvs': 'asv01|asv06'
                    },
            'seq2': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                    'stdv-kmer-per-region': 0,
                    'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batgirl;Gordon;Barbara',
                    'mapped-asvs': 'asv01|asv07',
                    },
            'seq3': {'num-regions': 2, 
                     'total-kmers-mapped': 3, 
                     'mean-kmer-per-region': 1.5,
                     'stdv-kmer-per-region': np.std([1, 2], ddof=1),
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batman;Wayne;',
                     'mapped-asvs': 'asv02|asv03|asv08'
                    },
            'seq4': {'num-regions': 1, 
                     'total-kmers-mapped': 1, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Civillian;Doctor;Thompson;Leslie',
                     'mapped-asvs': 'asv09'
                    },
            'seq5': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Robin;Grayson;Dick',
                     'mapped-asvs': 'asv04|asv05|asv10',
                    },
            'seq6': {'num-regions': 2, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Robin;Todd;Jason',
                     'mapped-asvs': 'asv04|asv05|asv11',
                    },
            })
        self.seq_summary.index.set_names('clean_name', inplace=True)
        self.table = pd.DataFrame(
            data=np.array([[20, 10, 10, 10, 10, 10, 10, 10, 10, 10]]),
            index=['s.1'],
            columns=['asv01', 'asv02', 'asv04', 'asv05', 'asv06', 
                     'asv07', 'asv08', 'asv09', 'asv10', 'asv11']
            ).T
        self.table = self.table / self.table.sum()
        self.freq = pd.Series(
            np.array([0.1818, 0.1818, 0.1818, 0.0909, 0.1818, 0.1818]),
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], 
                            name='clean_name'),
            name='sample.1'
            )
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')

        self.shared_long = pd.DataFrame(
            data=[['seq1', 'Bludhaven', 'seq1'],
                  ['seq1', 'Bludhaven', 'seq2'],
                  ['seq2', 'Bludhaven', 'seq1'],
                  ['seq2', 'Bludhaven', 'seq2'],
                  ['seq3', 'Bludhaven', 'seq3'],
                  ['seq5', 'Bludhaven', 'seq5'],
                  ['seq6', 'Bludhaven', 'seq6']],
            columns=['db-seq', 'region', 'value']
            )

    def test_reconstruct_counts(self):
        known_map = pd.DataFrame(
            data=[['seq1', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq2', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq3', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq4', 'CACCTCGTN', 'CACCTCGTN', 15],
                  ['seq5', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq6', 'WANTCAT', 'CACCTCGTN', 15]],
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], 
                            name='db-seq'),
            columns=['clean_name', 'first-fwd-primer', 'last-fwd-primer', 
                     'last-kmer-length'],
            )

        known_summary = pd.DataFrame.from_dict(orient='index', data={
            'seq1': {'num-regions': 2., 
                     'total-kmers-mapped': 2., 
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
                   os.path.join(self.base_dir, 'region1_counts.qza'), 0],
                  [os.path.join(self.base_dir, 'region2_db_map.qza'),
                   os.path.join(self.base_dir, 'region2_align.qza'),
                   os.path.join(self.base_dir, 'region2_counts.qza'), 1]],
            columns=['kmer-map', 'alignment-map', 'frequency-table', 
                     'region-order'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
            ))
        count_table, summary, mapping = \
            reconstruct_counts(manifest, debug=True, min_abund=1e-2)

        npt.assert_array_equal(
            count_table.matrix_data.todense(),
            np.array([[100,  50,   0,  50,  50, 50],
                      [100,  25, 100,  25,  25, 25],
                      [  0, 100, 100,   0,  50, 50]]).T
           )
        npt.assert_array_equal(
            np.array(list(count_table.ids(axis='sample'))),
            np.array(['sample1', 'sample2', 'sample3'])
        )
        npt.assert_array_equal(
            np.array(list(count_table.ids(axis='observation'))),
            np.array(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6']),
        )
        pdt.assert_frame_equal(known_map, mapping)
        pdt.assert_frame_equal(known_summary, summary.to_dataframe())


    def test_reconstruct_counts_unweighted(self):
        known_map = pd.DataFrame(
            data=[['seq1', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq2', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq3', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq4', 'CACCTCGTN', 'CACCTCGTN', 15],
                  ['seq5', 'WANTCAT', 'CACCTCGTN', 15],
                  ['seq6', 'WANTCAT', 'CACCTCGTN', 15]],
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], 
                            name='db-seq'),
            columns=['clean_name', 'first-fwd-primer', 'last-fwd-primer', 
                     'last-kmer-length'],
            )

        known_summary = pd.DataFrame.from_dict(orient='index', data={
            'seq1': {'num-regions': 1., 
                     'total-kmers-mapped': 2., 
                     'mean-kmer-per-region': 1.,
                     'stdv-kmer-per-region': 0.,
                     'mapped-asvs': 'asv01|asv06'
                    },
            'seq2': {'num-regions': 1, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                    'stdv-kmer-per-region': 0,
                    'mapped-asvs': 'asv01|asv07',
                    },
            'seq3': {'num-regions': 1, 
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
            'seq5': {'num-regions': 1, 
                     'total-kmers-mapped': 2, 
                     'mean-kmer-per-region': 1,
                     'stdv-kmer-per-region': 0,
                     'mapped-asvs': 'asv04|asv05|asv10',
                    },
            'seq6': {'num-regions': 1, 
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
                   os.path.join(self.base_dir, 'region1_counts.qza'), 0],
                  [os.path.join(self.base_dir, 'region2_db_map.qza'),
                   os.path.join(self.base_dir, 'region2_align.qza'),
                   os.path.join(self.base_dir, 'region2_counts.qza'), 1]],
            columns=['kmer-map', 'alignment-map', 'frequency-table', 
                     'region-order'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
            ))
        count_table, summary, mapping = \
            reconstruct_counts(manifest, debug=True, min_abund=1e-2, 
                               region_normalize=False)
        npt.assert_array_equal(
            count_table.matrix_data.todense(),
            np.array([[200, 100,   0,  50, 100, 100],
                      [200,  50, 200,  25,  50,  50],
                      [  0, 200, 200,   0, 100, 100]]).T * 1.
           )
        npt.assert_array_equal(
            np.array(list(count_table.ids(axis='sample'))),
            np.array(['sample1', 'sample2', 'sample3'])
        )
        npt.assert_array_equal(
            np.array(list(count_table.ids(axis='observation'))),
            np.array(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6']),
        )
        pdt.assert_frame_equal(known_map, mapping)
        pdt.assert_frame_equal(known_summary, summary.to_dataframe())



    def test_construct_align_mat(self):
        sequence_map = pd.Series({'seq1': 'seq1',
                                  'seq2': 'seq2',
                                  'seq3': 'seq3',
                                  'seq4': 'seq4',
                                  'seq5': 'seq5',
                                  'seq6': 'seq6'}, name='clean_name')
        sequence_map.index.set_names('db-seq', inplace=True)
        test_mat = _construct_align_mat(self.match2, 
                                        sequence_map,
                                        self.seq_summary).compute()
        # Rounds for the sake of sanity
        test_mat[self.float_cols] = test_mat[self.float_cols].round(4)
        test_mat[self.int_cols] = test_mat[self.int_cols].astype(int)

        pdt.assert_frame_equal(self.align2, test_mat)

    def test_count_mapping_degen(self):
        known = pd.DataFrame.from_dict(orient='index', data={
            'seq05': {'num-regions': 2., 
                      'total-kmers-mapped': 3., 
                      'mean-kmer-per-region': 1.5,
                      'stdv-kmer-per-region': np.std([2, 3], ddof=1),
                      },
            'seq06': {'num-regions': 2, 
                      'total-kmers-mapped': 2,
                      'mean-kmer-per-region': 1,
                      'stdv-kmer-per-region': 0},
            })
        known.index.set_names('clean_name', inplace=True)
        test = _count_mapping(self.pair, count_degen=True)
        pdt.assert_frame_equal(test, known)

    def test_count_mapping_no_degen(self):
        known = pd.DataFrame.from_dict(orient='index', data={
            'seq05': {'num-regions': 2., 
                      'total-kmers-mapped': 2., 
                      'mean-kmer-per-region': 1.,
                      'stdv-kmer-per-region': 0.,
                      },
            'seq06': {'num-regions': 2, 
                      'total-kmers-mapped': 2,
                      'mean-kmer-per-region': 1,
                      'stdv-kmer-per-region': 0},
            })
        known.index.set_names('clean_name', inplace=True)
        test = _count_mapping(self.pair, count_degen=False)
        pdt.assert_frame_equal(test, known)

    def test_detangle_names(self):
        long_ = pd.DataFrame(
            data=[['seq00', 'seq00', 0],
                  ['seq00', 'seq02', 1],
                  ['seq01', 'seq01', 0], 
                  ['seq01', 'seq02', 1], 
                  ['seq02', 'seq00', 0],
                  ['seq02', 'seq01', 1],
                  ['seq02', 'seq02', 2],
                  ['seq03', 'seq03', 0],
                  ['seq03', 'seq04', 1],
                  ['seq04', 'seq03', 0],
                  ['seq04', 'seq04', 1],
                  ['seq04', 'seq05', 2],
                  ['seq05', 'seq04', 0],
                  ['seq05', 'seq05', 1],
                  ['seq06', 'seq06', 1],
                  ['seq07', 'seq07', 0],
                  ['seq07', 'seq08', 1],
                  ['seq08', 'seq07', 0],
                  ['seq08', 'seq08', 1],
                  ['seq10', 'seq10', 0],
                  ['seq10', 'seq12', 1],
                  ['seq10', 'seq13', 2],
                  ['seq10', 'seq14', 3],
                  ['seq11', 'seq11', 0],
                  ['seq11', 'seq12', 1],
                  ['seq11', 'seq13', 2],
                  ['seq11', 'seq14', 3],
                  ['seq11', 'seq15', 4],
                  ['seq12', 'seq10', 0],
                  ['seq12', 'seq11', 1],
                  ['seq12', 'seq12', 2],
                  ['seq12', 'seq13', 3],
                  ['seq12', 'seq14', 4],
                  ['seq13', 'seq10', 0],
                  ['seq13', 'seq11', 1],
                  ['seq13', 'seq12', 2],
                  ['seq13', 'seq13', 3],
                  ['seq13', 'seq14', 4],
                  ['seq14', 'seq10', 0],
                  ['seq14', 'seq11', 1],
                  ['seq14', 'seq12', 2],
                  ['seq14', 'seq13', 3],
                  ['seq14', 'seq14', 4],
                  ['seq15', 'seq11', 0],
                  ['seq15', 'seq15', 1]
                  ],
            columns=['db-seq', 'clean_name', 'counter'],
            )
        known = pd.Series(
            data=['seq00', 'seq01', 'seq02', 'seq03', 
                  'seq04',  'seq05', 'seq06', 'seq07|seq08', 
                  'seq07|seq08', 'seq10', 'seq11', 'seq12|seq13|seq14', 
                  'seq12|seq13|seq14', 'seq12|seq13|seq14', 'seq15'],
            index=pd.Index(['seq00', 'seq01', 'seq02', 'seq03', 'seq04', 
                            'seq05', 'seq06', 'seq07', 'seq08', 'seq10',
                            'seq11', 'seq12', 'seq13', 'seq14', 'seq15'],
                            name='db-seq'),
            name='clean_name'
            )
        known.index.set_names(['db-seq'], inplace=True)
        test = _detangle_names(long_) 
        pdt.assert_series_equal(known, test)

    def test_expand_duplicate_sequences(self):
        original = pd.DataFrame(data=[['1', '2|3', '3|4|5', '6'],
                                      [10, 20, 30, 40]],
                                index=['ids', 'counts']).T
        known = pd.DataFrame([[ '1', '2', '3', '3', '4', '5', '6'],
                              [10, 20, 20, 30, 30, 30, 40]],
                             index=['ids', 'counts']).T
        test = _expand_duplicate_sequences(original, id_col='ids')

        pdt.assert_frame_equal(test, known)

    def test_get_clean(self):
        known = pd.DataFrame(
            data=[['seq1', 'Bludhaven', 'seq1|seq2'],
                  ['seq2', 'Bludhaven', 'seq1|seq2'],
                  ['seq3', 'Bludhaven', 'seq3'],
                  ['seq5', 'Bludhaven', 'seq5'],
                  ['seq6', 'Bludhaven', 'seq6']],
            columns=['db-seq', 'region', 'value']
            )
        test = _get_clean(self.shared_long).compute()
        pdt.assert_frame_equal(test, known)

    def test_get_shared(self):
        test = _get_shared_seqs(self.kmer1.reset_index()).compute()
        test.sort_values(['db-seq', 'value'], inplace=True)
        pdt.assert_frame_equal(test.reset_index(drop=True), 
                               self.shared_long)

    def test_get_unique_kmers(self):
        test = pd.Series(['seq00', 'seq00|seq01@001', 'seq01@002', 
                          'seq01|seq02'])
        known = np.array(['seq00', 'seq01', 'seq02'])
        test = _get_unique_kmers(test)
        npt.assert_array_equal(test, known)

    def test_scale_relative_abundance(self):
        known = pd.DataFrame(
            data=np.array([[10., 10., 10., 10., 10., 10.]]),
            index=['sample.1'],
            columns=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
                             name='clean_name'),
            ).T
        counts = pd.DataFrame(
            data=np.array([[20, 10, 10, 10, 10, 10, 10, 10, 10, 10]]),
            index=['sample.1'],
            columns=['asv01', 'asv02', 'asv04', 'asv05', 'asv06', 
                     'asv07', 'asv08', 'asv09', 'asv10', 'asv11']
            ).T
        test = _scale_relative_abundance(
            pd.concat([self.align1, self.align2]),
            relative=pd.DataFrame(self.freq),
            counts=counts,
            seq_summary=self.seq_summary,
            )
        pdt.assert_frame_equal(known, test)

    def test_solve_iterative_noisy(self):
        known = pd.DataFrame(
            data=np.array([[1] * 6]) / 6,
            index=['s.1'],
            columns=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], 
                             name='clean_name'),
        ).T
        test = _solve_iterative_noisy(
            align_mat=pd.concat([self.align1, self.align2]),
            table=self.table / self.table.sum(),
            seq_summary=self.seq_summary,
            )
        pdt.assert_frame_equal(known, test)

    def test_solve_ml_em_iterative_1_sample(self):
        abund = np.array([0.18181818, 0.09090909, 0.09090909,
                          0.09090909,  0.09090909, 0.09090909, 
                          0.09090909, 0.09090909,  0.09090909, 
                          0.09090909])
        align = np.array([
            [0.4638, 0.4638, 0,      0,      0,      0     ],
            [0,      0,      0.4638, 0,      0,      0     ],
            [0,      0,      0,      0,      0.4638, 0.0008],
            [0,      0,      0,      0,      0.0008, 0.4638],
            [0.4638, 0,      0 ,     0,      0,      0     ],
            [0,      0.4638, 0,      0,      0,      0     ],
            [0,      0,      0.4638, 0,      0,      0     ],
            [0,      0,      0,      0.9276, 0,      0     ],
            [0,      0,      0,      0,      0.4638, 0     ],
            [0,      0,      0,      0,      0,      0.4638]])

        t_freq = _solve_ml_em_iterative_1_sample(
            align=align, 
            abund=abund,
            align_kmers=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 
                                  'seq6'], 
                                 name='clean_name'),
            sample='sample.1'
            )
        pdt.assert_series_equal(self.freq, t_freq.compute(), 
                                check_less_precise=True)

    def test_sort_untidy(self):
        pass

    def test_tidy_sequence_set(self):
        clean_kmer = dd.from_pandas(npartitions=1, data=pd.DataFrame(
            columns=['db-seq', 'region', 'shared-set'],
            data=[['seq00', 0, {'seq00'}], ['seq00', 1, {'seq00'}],
                  ['seq01', 0, {'seq01', 'seq02'}], ['seq01', 1, {'seq01'}],
                  ['seq02', 0, {'seq01', 'seq02'}], ['seq02', 1, {'seq02'}],
                  ['seq05', 0, {'seq05'}], ['seq05', 1, {'seq05', 'seq06'}],
                  ['seq06', 0, {'seq06'}], ['seq06', 1, {'seq05', 'seq06'}],
                  ['seq07', 0, {'seq07'}], 
                  ['seq08', 0, {'seq08', 'seq09'}],
                  ['seq09', 0, {'seq08', 'seq09'}], ['seq09', 1, {'seq09'}],
                  ['seq10', 0, {'seq10', 'seq11'}], ['seq10', 1, {'seq10', 'seq11'}],
                  ['seq11', 0, {'seq10', 'seq11'}], ['seq11', 1, {'seq10', 'seq11'}],
                  ['seq12', 0, {'seq12', 'seq13', 'seq14'}], ['seq12', 1, {'seq12', 'seq14', 'seq15'}],
                  ['seq13', 0, {'seq12', 'seq13', 'seq14'}], ['seq13', 1, {'seq13', 'seq14', 'seq15'}],
                  ['seq14', 0, {'seq12', 'seq13', 'seq14'}], 
                  ['seq14', 1, {'seq12', 'seq13', 'seq14', 'seq15'}],
                  ['seq15', 1, {'seq12', 'seq13', 'seq14', 'seq15'}],
                  ['seq16', 0, {'seq16'}], ['seq16', 1, {'seq16', 'seq17'}],
                  ['seq17', 0, {'seq17', 'seq18'}], ['seq17', 1, {'seq16', 'seq17'}],
                  ['seq18', 0, {'seq17', 'seq18'}], 
                  ['seq19', 0, {'seq19', 'seq20'}], ['seq20', 0, {'seq19', 'seq20'}],
                  ['seq21', 0, {'seq21'}], ['seq21', 1, {'seq21', 'seq22', 'seq23'}],
                  ['seq22', 1, {'seq21', 'seq22', 'seq23'}],
                  ['seq23', 1, {'seq21', 'seq22', 'seq23'}], ['seq23', 2, {'seq23'}]]
            ))
        clean_kmer['tidy'] = False
        known_seqs = set(['seq00', 'seq01', 'seq02', 'seq05', 'seq06',
                          'seq07', 'seq09', 'seq10', 'seq11', 'seq16', 'seq17', 
                          'seq19', 'seq20', 'seq21', 'seq23'])
        known_kmers = pd.DataFrame(
            columns=['db-seq', 'region', 'shared-set', 'tidy'],
            data=[['seq00', 0, 'seq00', True], 
                  ['seq00', 1, 'seq00', True],
                  ['seq01', 0, 'seq01|seq02', True], 
                  ['seq01', 1, 'seq01', True],
                  ['seq02', 0, 'seq01|seq02', True],
                  ['seq02', 1, 'seq02', True],
                  ['seq05', 0, 'seq05', True],
                  ['seq05', 1, 'seq05|seq06', True],
                  ['seq06', 0, 'seq06', True],
                  ['seq06', 1, 'seq05|seq06', True],
                  ['seq07', 0, 'seq07', True], 
                  ['seq08', 0, 'seq08', False],
                  ['seq09', 0, 'seq08|seq09', True],
                  ['seq09', 1, 'seq09', True],
                  ['seq10', 0, 'seq10|seq11', True],
                  ['seq10', 1, 'seq10|seq11', True],
                  ['seq11', 0, 'seq10|seq11', True], 
                  ['seq11', 1, 'seq10|seq11', True],
                  ['seq12', 0, 'seq12|seq13|seq14', False], 
                  ['seq12', 1, 'seq12|seq14|seq15', False],
                  ['seq13', 0, 'seq12|seq13|seq14', False], 
                  ['seq13', 1, 'seq13|seq14|seq15', False],
                  ['seq14', 0, 'seq12|seq13|seq14', False], 
                  ['seq14', 1, 'seq12|seq13|seq14|seq15', False],
                  ['seq15', 1, 'seq12|seq13|seq14|seq15', False],
                  ['seq16', 0, 'seq16', True], 
                  ['seq16', 1, 'seq16|seq17', True],
                  ['seq17', 0, 'seq17|seq18', True], 
                  ['seq17', 1, 'seq16|seq17', True],
                  ['seq18', 0, 'seq18', False], 
                  ['seq19', 0, 'seq19|seq20', True], 
                  ['seq20', 0, 'seq19|seq20', True],
                  ['seq21', 0, 'seq21', True], 
                  ['seq21', 1, 'seq21|seq22|seq23', True],
                  ['seq22', 1, 'seq22', False],
                  ['seq23', 1, 'seq21|seq22|seq23', True], 
                  ['seq23', 2, 'seq23', True]]
            )
        test_kmers, test_seqs = _tidy_sequence_set(clean_kmer, set([]))
        self.assertEqual(test_seqs, known_seqs)
        # Makes testing easier
        test_kmers = test_kmers.compute()
        test_kmers['shared-set'] =  \
            test_kmers['shared-set'].apply(lambda x: "|".join(sorted(x)))
        pdt.assert_frame_equal(known_kmers, test_kmers)

    def test_untangle_database_ids(self):
        matches = dd.from_pandas(npartitions=1, data=pd.DataFrame(
            data=np.array([['seq00', 'seq00', 'seq00', '0'],
                           ['seq01', 'seq01', 'seq01|seq02', '0'],
                           ['seq02', 'seq02', 'seq01|seq02', '0'],
                           ['seq03', 'seq03@001', 'seq03@001', '0'],
                           ['seq03', 'seq03@002', 'seq03@002', '0'],
                           ['seq04', 'seq04@001', 'seq04@001', '0'],
                           ['seq04', 'seq04@002', 'seq04@002', '0'],
                           ['seq05', 'seq05', 'seq05', '0'],
                           ['seq06', 'seq06', 'seq06', '0'],
                           ['seq07', 'seq07', 'seq07', '0'],
                           ['seq08', 'seq08', 'seq08|seq09', '0'],
                           ["seq09", 'seq09', 'seq08|seq09', '0'],
                           ['seq10', 'seq10', 'seq10|seq11', '0'],
                           ['seq11', 'seq11', 'seq10|seq11', '0'],
                           ['seq12', 'seq12', 'seq12|seq13|seq14', '0'],
                           ['seq13', 'seq13', 'seq12|seq13|seq14', '0'],
                           ['seq14', 'seq14', 'seq12|seq13|seq14', '0'],
                           ['seq16', 'seq16', 'seq16', '0'],
                           ['seq17', 'seq17', 'seq17|seq18', '0'],
                           ['seq18', 'seq18', 'seq17|seq18', '0'],
                           ['seq19', 'seq19', 'seq19|seq20', '0'],
                           ['seq20', 'seq20', 'seq19|seq20', '0'],
                           ['seq21', 'seq21', 'seq21', '0'],
                           ['seq00', 'seq00', 'seq00', '1'],
                           ['seq01', 'seq01', 'seq01', '1'],
                           ['seq02', 'seq02', 'seq02', '1'],
                           ['seq03', 'seq03', 'seq03', '1'],
                           ['seq04', 'seq04@001', 'seq04@001', '1'],
                           ['seq04', 'seq04@002', 'seq04@002', '1'],
                           ['seq04', 'seq04@003', 'seq04@003', '1'],
                           ['seq05', 'seq05@001', 'seq05@001', '1'],
                           ['seq05', 'seq05@002', 'seq05@002|seq06', '1'],
                           ['seq06', 'seq06', 'seq05@002|seq06', '1'],
                           ['seq09', 'seq09', 'seq09', '1'],
                           ['seq10', 'seq10', 'seq10|seq11', '1'],
                           ['seq11', 'seq11', 'seq10|seq11', '1'],
                           ['seq12', 'seq12@0001', 'seq12@0001', '1'],
                           ['seq12', 'seq12@0002', 'seq12@0002|seq14|seq15', '1'],
                           ['seq13', 'seq13@0001', 'seq13@0001|seq14|seq15', '1'],
                           ['seq13', 'seq13@0002', 'seq13@002', '1'],
                           ['seq14', 'seq14', 'seq12@0002|seq13@0001|seq14|seq15', '1'],
                           ['seq15', 'seq15', 'seq12@0002|seq13@0001|seq14|seq15', '1'],
                           ['seq16', 'seq16', 'seq16|seq17',  '1'],
                           ['seq17', 'seq17', 'seq16|seq17', '1'],
                           ['seq21', 'seq21', 'seq21|seq22|seq23', '1'],
                           ['seq22', 'seq22', 'seq21|seq22|seq23', '1'],
                           ['seq23', 'seq23', 'seq21|seq22|seq23', '1'],
                           ['seq22', 'seq22', 'seq22', '2'],
                           ], dtype=object),
            columns=['db-seq', 'seq-name', 'kmer', 'region'],
            ))
        matches['region'] =  matches['region'].astype(int)
        known_seq = pd.Series({'seq00': 'seq00', 
                                'seq01': 'seq01',
                                'seq02': 'seq02',
                                'seq03': 'seq03',
                                'seq04': 'seq04',
                                'seq05': 'seq05',
                                'seq06': 'seq06',
                                'seq07': 'seq07',
                                'seq08': 'seq08',
                                'seq09': 'seq09',
                                'seq10': 'seq10|seq11',
                                'seq11': 'seq10|seq11',
                                'seq12': 'seq12',
                                'seq13': 'seq13',
                                'seq14': 'seq14',
                                'seq15': 'seq15',
                                'seq16': 'seq16',
                                'seq17': 'seq17',
                                'seq18': 'seq18',
                                'seq19': 'seq19|seq20',
                                'seq20': 'seq19|seq20',
                                'seq21': 'seq21',
                                'seq22': 'seq22',
                                'seq23': 'seq23'
                                }, name='clean_name')
        known_seq.index.set_names('db-seq', inplace=True)
        # Generates the renaming
        seq_ = _untangle_database_ids(matches, num_regions=3)

        pdt.assert_series_equal(known_seq, seq_.sort_index())


if __name__ == '__main__':
    main()
