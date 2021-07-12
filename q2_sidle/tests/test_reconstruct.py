from unittest import TestCase, main

import copy
import os
import warnings

import biom
import dask
from dask.delayed import Delayed
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio

from qiime2 import Artifact, Metadata
from qiime2.plugin import ValidationError

from q2_sidle._reconstruct import (reconstruct_counts,
                                   _construct_align_mat,
                                   _expand_duplicate_sequences,
                                   _get_db_and_clean,
                                   _scale_relative_abundance,
                                   _solve_ml_em_iterative_1_sample,
                                   _solve_iterative_noisy,
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
            data=np.array([['seq1', 'asv01', 'Bludhaven', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq2', 'asv01', 'Bludhaven', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq3', 'asv02', 'Bludhaven', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq5', 'asv04', 'Bludhaven', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq5', 'asv05', 'Bludhaven', 2, 1, 15, 0.00155, 0.92757, 0.0, 0.92757, False, 0.0008],
                           ['seq6', 'asv04', 'Bludhaven', 2, 1, 15, 0.00155, 0.92757, 0.0, 0.92757, False, 0.0008],
                           ['seq6', 'asv05', 'Bludhaven', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638]], dtype=object),
            columns=['clean_name', 'asv', 'region', 'max-mismatch', 'mismatch', 'length', 
                     'match_prob', 'perfect_prob',  
                     'max_error_prob', 'error_thresh', 'perf_match', 'norm']
            )
        self.align2 = pd.DataFrame(
            data=np.array([['seq1', 'asv06', 'Gotham', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq2', 'asv07', 'Gotham', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq3', 'asv08', 'Gotham', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq4', 'asv09', 'Gotham', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.9276],
                           ['seq5', 'asv10', 'Gotham', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq6', 'asv11', 'Gotham', 2, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638]], dtype=object),
            columns=['clean_name', 'asv', 'region', 'max-mismatch', 'mismatch', 'length',
                     'match_prob', 'perfect_prob',  
                     'max_error_prob', 'error_thresh', 'perf_match', 'norm']
            )
        self.int_cols = ['length', 'mismatch',  'max-mismatch']
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
        self.freq = biom.Table(
            np.array([[0.1818, 0.1818, 0.1818, 0.0909, 0.1818, 0.1818]]).T,
            observation_ids=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
            sample_ids=['sample.1']
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

        self.manifest = Metadata(pd.DataFrame(
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
        self.datbase_map = pd.DataFrame(
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
        db_summary = pd.DataFrame.from_dict(orient='index', data={
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
        db_summary.index.set_names('feature-id', inplace=True)
        self.database_summary = Metadata(db_summary)

    def test_reconstruct_counts(self):
        count_table = reconstruct_counts(
              region=['Bludhaven', 'Gotham'],
              regional_alignment=[ts.region1_align.view(Delayed), 
                                  ts.region2_align.view(Delayed)],
              regional_table=[ts.region1_counts.view(biom.Table),
                              ts.region2_counts.view(biom.Table)],
              database_map=self.datbase_map['clean_name'],
              database_summary = self.database_summary,
              debug=True, 
              min_counts=10,
              min_abund=1e-2)
        npt.assert_array_equal(
            np.array(count_table.matrix_data.todense()),
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

    def test_reconstruct_counts_unweighted(self):
        count_table = reconstruct_counts(
            region=['Bludhaven', 'Gotham'],
            regional_alignment=[ts.region1_align.view(Delayed).copy(), 
                                ts.region2_align.view(Delayed).copy()],
            regional_table=[ts.region1_counts.view(biom.Table),
                            ts.region2_counts.view(biom.Table)],
            database_map=self.datbase_map['clean_name'],
            database_summary = self.database_summary,
            debug=True, 
            min_counts=10,
            min_abund=1e-2, 
            region_normalize='unweighted'
            )
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

    def test_reconstruct_counts_align_drop_samples_error(self):
        with self.assertRaises(ValueError) as err:
            count_table = reconstruct_counts(
                region=['Bludhaven', 'Gotham'],
                regional_alignment=[ts.region1_align.view(pd.DataFrame).copy(), 
                                    ts.region2_align.view(pd.DataFrame).copy()],
                regional_table=[ts.region1_counts.view(biom.Table),
                              ts.region2_counts.view(biom.Table)],
                database_map=self.datbase_map['clean_name'],
                database_summary = self.database_summary,
                debug=True, 
                min_abund=1e-2
                )
        self.assertEqual(
            str(err.exception),
            'None of the samples have more than the 1000 total '
            'sequences required for reconstruction.' 
            )

    def test_reconstruct_counts_align_drop_samples_warning(self):
        with warnings.catch_warnings(record=True) as w:
            count_table = reconstruct_counts(
                region=['Bludhaven', 'Gotham'],
                regional_alignment=[ts.region1_align.view(pd.DataFrame).copy(), 
                                    ts.region2_align.view(pd.DataFrame).copy()],
                regional_table=[ts.region1_counts.view(biom.Table),
                              ts.region2_counts.view(biom.Table)],
                database_map=self.datbase_map['clean_name'],
                database_summary = self.database_summary,
                debug=True, 
                min_counts=590,
                min_abund=1e-2
                )
        self.assertTrue(issubclass(w[-1].category, UserWarning))
        self.assertEqual(str(w[-1].message), 
                         'There are 2 samples with fewer than 590 '
                         'total reads. These samples will be discarded.')

    def test_construct_align_mat(self):
        sequence_map = pd.Series({'seq1': 'seq1',
                                  'seq2': 'seq2',
                                  'seq3': 'seq3',
                                  'seq4': 'seq4',
                                  'seq5': 'seq5',
                                  'seq6': 'seq6'}, name='clean_name')
        sequence_map.index.set_names('db-seq', inplace=True)
        test_mat = _construct_align_mat([ts.region2_align.view(Delayed)],
                                        sequence_map=sequence_map,
                                        seq_summary=self.seq_summary)
        # Rounds for the sake of sanity
        test_mat[self.float_cols] = test_mat[self.float_cols].round(4)
        test_mat[self.int_cols] = test_mat[self.int_cols].astype(int)

        pdt.assert_frame_equal(self.align2, test_mat)

    def test_expand_duplicate_sequences(self):
        original = pd.DataFrame(data=[['1', '2|3', '3|4|5', '6'],
                                      [10, 20, 30, 40]],
                                index=['ids', 'counts']).T
        known = pd.DataFrame([[ '1', '2', '3', '3', '4', '5', '6'],
                              [10, 20, 20, 30, 30, 30, 40]],
                             index=['ids', 'counts']).T
        test = _expand_duplicate_sequences(original, id_col='ids')

        pdt.assert_frame_equal(test, known)

    def test_get_db_and_clean(self):
        in_ = pd.DataFrame(np.array([['Robin@001', 'RedHood', 'Spoiler']]).T,
                           columns=['hero'])
        match = pd.Series({"Robin": 'DickGrayson', 'RedHood': "JasonTodd", 
                           'Spoiler': 'StephanieBrown'})
        known = pd.DataFrame(
            data=np.array([['Robin', 'Robin@001', 'DickGrayson'],
                           ['RedHood', 'RedHood', 'JasonTodd'],
                           ['Spoiler', 'Spoiler', 'StephanieBrown'],
                           ], dtype=str),
            columns=['db-seq', 'hero', 'civillian-id'],
            )
        test = _get_db_and_clean(in_, match, 'hero', 'civillian-id')
        pdt.assert_frame_equal(known, test)


    # def test_scale_relative_abundance_average(self):
    #     known = biom.Table(
    #         data=np.array([[10., 10., 10., 10., 10., 10.]]).T,
    #         sample_ids=['sample.1'],
    #         observation_ids=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
    #         )
    #     counts = pd.DataFrame(
    #         data=np.array([[20, 10, 10, 10, 10, 10, 10, 10, 10, 10]]),
    #         index=['sample.1'],
    #         columns=['asv01', 'asv02', 'asv04', 'asv05', 'asv06', 
    #                  'asv07', 'asv08', 'asv09', 'asv10', 'asv11']
    #         ).T
    #     test = _scale_relative_abundance(
    #         pd.concat([self.align1, self.align2]),
    #         relative=self.freq,
    #         counts=counts,
    #         num_regions=2,
    #         region_normalize='average',
    #         seq_summary=self.seq_summary,
    #         )
    #     npt.assert_array_equal(test.ids(axis='observation'),
    #                            known.ids(axis='observation'))
    #     npt.assert_array_equal(test.ids(axis='sample'),
    #                            known.ids(axis='sample'))
    #     npt.assert_array_equal(known.matrix_data.todense(),
    #                            test.matrix_data.todense())

    # def test_scale_relative_abundance_weighted(self):
    #     known = biom.Table(
    #         data=np.array([[20., 20., 20., 20., 20., 20.]]).T,
    #         sample_ids=['sample.1'],
    #         observation_ids=['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
    #         )
    #     counts = pd.DataFrame(
    #         data=np.array([[20, 10, 10, 10, 10, 10, 10, 10, 10, 10]]),
    #         index=['sample.1'],
    #         columns=['asv01', 'asv02', 'asv04', 'asv05', 'asv06', 
    #                  'asv07', 'asv08', 'asv09', 'asv10', 'asv11']
    #         ).T
    #     test = _scale_relative_abundance(
    #         pd.concat([self.align1, self.align2]),
    #         relative=self.freq,
    #         counts=counts,
    #         num_regions=2,
    #         region_normalize='weighted',
    #         seq_summary=self.seq_summary,
    #         )
    #     npt.assert_array_equal(test.ids(axis='observation'),
    #                            known.ids(axis='observation'))
    #     npt.assert_array_equal(test.ids(axis='sample'),
    #                            known.ids(axis='sample'))
    #     npt.assert_array_equal(known.matrix_data.todense(),
    #                            test.matrix_data.todense())

    # def test_solve_iterative_noisy(self):
    #     known = pd.DataFrame(
    #         data=np.array([[1] * 6]) / 6,
    #         index=['s.1'],
    #         columns=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'], 
    #                          name='clean_name'),
    #     ).T
    #     test = _solve_iterative_noisy(
    #         align_mat=pd.concat([self.align1, self.align2]),
    #         table=self.table / self.table.sum(),
    #         seq_summary=self.seq_summary,
    #         )
    #     npt.assert_array_equal(['s.1'], list(test.ids(axis='sample')))
    #     npt.assert_array_equal(np.array(['seq1', 'seq2', 'seq3', 
    #                                      'seq4', 'seq5', 'seq6']),
    #                            list(test.ids(axis='observation')))
    #     npt.assert_almost_equal(
    #         np.array([[1] * 6]) / 6,
    #         test.matrix_data.todense().T
    #     )

    # def test_solve_ml_em_iterative_1_sample(self):
    #     abund = np.array([0.18181818, 0.09090909, 0.09090909,
    #                       0.09090909,  0.09090909, 0.09090909, 
    #                       0.09090909, 0.09090909,  0.09090909, 
    #                       0.09090909])
    #     align = np.array([
    #         [0.4638, 0.4638, 0,      0,      0,      0     ],
    #         [0,      0,      0.4638, 0,      0,      0     ],
    #         [0,      0,      0,      0,      0.4638, 0.0008],
    #         [0,      0,      0,      0,      0.0008, 0.4638],
    #         [0.4638, 0,      0 ,     0,      0,      0     ],
    #         [0,      0.4638, 0,      0,      0,      0     ],
    #         [0,      0,      0.4638, 0,      0,      0     ],
    #         [0,      0,      0,      0.9276, 0,      0     ],
    #         [0,      0,      0,      0,      0.4638, 0     ],
    #         [0,      0,      0,      0,      0,      0.4638]])

    #     t_freq = _solve_ml_em_iterative_1_sample(
    #         align=align, 
    #         abund=abund,
    #         align_kmers=pd.Index(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 
    #                               'seq6'], 
    #                              name='clean_name').values,
    #         sample='sample.1'
    #         )
    #     npt.assert_array_equal(['sample.1'], list(t_freq.ids(axis='sample')))
    #     npt.assert_array_equal(np.array(['seq1', 'seq2', 'seq3', 
    #                                      'seq4', 'seq5', 'seq6']),
    #                            list(t_freq.ids(axis='observation')))
    #     npt.assert_array_equal(
    #         np.array([[0.1818, 0.1818, 0.1818, 0.0909, 0.1818, 0.1818]]).T,
    #         t_freq.matrix_data.todense().round(4)
    #     )

if __name__ == '__main__':
    main()
