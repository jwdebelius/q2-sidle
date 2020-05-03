from unittest import TestCase, main

import os

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.util.testing as pdt
import skbio
import sparse as sp

from qiime2 import Artifact, Metadata
from qiime2.plugin import ValidationError

from q2_sidle._reconstruct import (reconstruct_counts,
                                   _build_id_set,
                                   _check_manifest,
                                   _construct_align_mat,
                                   _count_mapping,
                                   _expand_duplicate_sequences,
                                   _get_unique_kmers,
                                   _map_id_set,
                                   _read_manifest_files,
                                   _solve_iterative_noisy,
                                   _solve_ml_em_iterative_1_sample,
                                   _untangle_database_ids,
                                   _scale_relative_abundance
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
            columns=['db_seq', 'region', '0', '1'],
            )
        self.id_set = [np.array(['seq05', 'seq06'])]
        self.pair = pd.DataFrame(
            data=[['seq05', 0, 'seq05', 'seq05'],
                  ['seq05@001', 1, 'seq05', 'seq05'],
                  ['seq05@002', 1, 'seq05', 'seq05'],
                  ['seq06', 0, 'seq06', 'seq06'],
                  ['seq06', 1, 'seq06', 'seq06']
                  ],
            columns=['kmer', 'region', 'db_seq', 'clean_name'],
            )
        self.match1 = pd.DataFrame(
            data=np.array([['seq1|seq2', 'asv01', 0, 15, 'Bludhaven'],
                           ['seq3@0001', 'asv02', 0, 15, 'Bludhaven'],
                           ['seq3@0002', 'asv02', 1, 15, 'Bludhaven'],
                           ['seq5', 'asv04', 0, 15, 'Bludhaven'],
                           ['seq5', 'asv05', 1, 15, 'Bludhaven'],
                           ['seq6', 'asv04', 1, 15, 'Bludhaven'],
                           ['seq6', 'asv05', 0, 15, 'Bludhaven']], 
                           dtype=object),
            columns=['kmer', 'asv', 'mismatch', 'length', 'region'],
            )
        self.match1[['length', 'mismatch']] = \
            self.match1[['length', 'mismatch']].astype(int)
        self.match2 = match2 = pd.DataFrame(
            data=np.array([['seq1', 'asv06', 0, 15, 'Gotham'],
                           ['seq2', 'asv07', 0, 15, 'Gotham'],
                           ['seq3', 'asv08', 0, 15, 'Gotham'],
                           ['seq4', 'asv09', 0, 15, 'Gotham'],
                           ['seq5', 'asv10', 0, 15, 'Gotham'],
                           ['seq6', 'asv11', 0, 15, 'Gotham']],
                           dtype=object),
            columns=['kmer', 'asv', 'mismatch', 'length', 'region']
            )
        self.match2[['length', 'mismatch']] = \
            self.match2[['length', 'mismatch']].astype(int)
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
            data=np.array([['seq1', 'asv06', 1, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq2', 'asv07', 1, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq3', 'asv08', 1, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq4', 'asv09', 1, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.9276],
                           ['seq5', 'asv10', 1, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638],
                           ['seq6', 'asv11', 1, 0, 15, 0.92757, 0.92757, 0.0, 0.92757,  True, 0.4638]], dtype=object),
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
            'seq1': {'num_regions': 2, 
                     'total_kmers_mapped': 2, 
                     'mean_kmer_per_region': 1,
                     'stdv_kmer_per_region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batgirl;Gordon;Barbara',
                     'mapped_asvs': 'asv01|asv06'
                    },
            'seq2': {'num_regions': 2, 
                     'total_kmers_mapped': 2, 
                     'mean_kmer_per_region': 1,
                    'stdv_kmer_per_region': 0,
                    'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batgirl;Gordon;Barbara',
                    'mapped_asvs': 'asv01|asv07',
                    },
            'seq3': {'num_regions': 2, 
                     'total_kmers_mapped': 3, 
                     'mean_kmer_per_region': 1.5,
                     'stdv_kmer_per_region': np.std([1, 2], ddof=1),
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Batman;Wayne;',
                     'mapped_asvs': 'asv02|asv03|asv08'
                    },
            'seq4': {'num_regions': 1, 
                     'total_kmers_mapped': 1, 
                     'mean_kmer_per_region': 1,
                     'stdv_kmer_per_region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Civillian;Doctor;Thompson;Leslie',
                     'mapped_asvs': 'asv09'
                    },
            'seq5': {'num_regions': 2, 
                     'total_kmers_mapped': 2, 
                     'mean_kmer_per_region': 1,
                     'stdv_kmer_per_region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Robin;Grayson;Dick',
                     'mapped_asvs': 'asv04|asv05|asv10',
                    },
            'seq6': {'num_regions': 2, 
                     'total_kmers_mapped': 2, 
                     'mean_kmer_per_region': 1,
                     'stdv_kmer_per_region': 0,
                     'taxonomy': 'DCU;Superhero;Gotham;Batfamily;Robin;Todd;Jason',
                     'mapped_asvs': 'asv04|asv05|asv11',
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

    def test_reconstruct_counts(self):
        known_map = pd.Series({'seq1': 'seq1', 'seq2': 'seq2', 'seq3': 'seq3', 
                              'seq4': 'seq4', 'seq5': 'seq5', 'seq6': 'seq6'}, 
                              name='clean_name')
        known_map.index.set_names('db_seq', inplace=True)

        known_summary = pd.DataFrame.from_dict(orient='index', data={
            'seq1': {'num_regions': 2., 
                     'total_kmers_mapped': 2., 
                     'mean_kmer_per_region': 1.,
                     'stdv_kmer_per_region': 0.,
                     'mapped_asvs': 'asv01|asv06'
                    },
            'seq2': {'num_regions': 2, 
                     'total_kmers_mapped': 2, 
                     'mean_kmer_per_region': 1,
                    'stdv_kmer_per_region': 0,
                    'mapped_asvs': 'asv01|asv07',
                    },
            'seq3': {'num_regions': 2, 
                     'total_kmers_mapped': 3, 
                     'mean_kmer_per_region': 1.5,
                     'stdv_kmer_per_region': np.std([1, 2], ddof=1),
                     'mapped_asvs': 'asv02|asv03|asv08'
                    },
            'seq4': {'num_regions': 1, 
                     'total_kmers_mapped': 1, 
                     'mean_kmer_per_region': 1,
                     'stdv_kmer_per_region': 0,
                     'mapped_asvs': 'asv09'
                    },
            'seq5': {'num_regions': 2, 
                     'total_kmers_mapped': 2, 
                     'mean_kmer_per_region': 1,
                     'stdv_kmer_per_region': 0,
                     'mapped_asvs': 'asv04|asv05|asv10',
                    },
            'seq6': {'num_regions': 2, 
                     'total_kmers_mapped': 2, 
                     'mean_kmer_per_region': 1,
                     'stdv_kmer_per_region': 0,
                     'mapped_asvs': 'asv04|asv05|asv11',
                    },
            })
        known_summary.index.set_names('feature-id', inplace=True)
        
        manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(self.base_dir, 'region1_db_map.qza'),
                   os.path.join(self.base_dir, 'region1_align.qza'),
                   os.path.join(self.base_dir, 'region1_counts.qza')],
                  [os.path.join(self.base_dir, 'region2_db_map.qza'),
                   os.path.join(self.base_dir, 'region2_align.qza'),
                   os.path.join(self.base_dir, 'region2_counts.qza')]],
            columns=['kmer-map', 'alignment-map', 'frequency-table'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
            ))
        count_table, summary, mapping = \
            reconstruct_counts(manifest, debug=True)
        npt.assert_array_equal(
            count_table.matrix_data.todense(),
            np.array([[100,  50,   0,  50,  50, 50],
                      [100,  25, 100,  25,  25, 25],
                      [  0, 100, 100,   0,  50, 50]]).T * 1.
           )
        npt.assert_array_equal(
            np.array(list(count_table.ids(axis='sample'))),
            np.array(['sample1', 'sample2', 'sample3'])
        )
        npt.assert_array_equal(
            np.array(list(count_table.ids(axis='observation'))),
            np.array(['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6']),
        )
        pdt.assert_series_equal(known_map, mapping)
        pdt.assert_frame_equal(known_summary, summary.to_dataframe())

    def test_build_id_set(self):
        test = _build_id_set(['seq05', 'seq06'], 
                             self.tangle.drop(columns=['region']),
                             block_size=3)
        self.assertTrue(isinstance(test, list))
        self.assertEqual(len(test), 1)
        npt.assert_array_equal(self.id_set[0], test[0])

    def test_map_id_set(self):
        tangled = pd.DataFrame(
            data=[['seq00', 0, 'seq00',  np.nan,  np.nan],
                  ['seq00', 1, 'seq00', 'seq01',  np.nan],
                  ['seq01', 1, 'seq00', 'seq01',  np.nan],
                  ['seq01', 2, 'seq01',  np.nan,  np.nan],
                  ['seq03', 0, 'seq03', 'seq04',  np.nan],
                  ['seq03', 1, 'seq03', 'seq04', 'seq05'],
                  ['seq04', 0, 'seq03', 'seq04',  np.nan],
                  ['seq04', 1, 'seq03', 'seq04', 'seq05'],
                  ['seq04', 2, 'seq04', 'seq05',  np.nan],
                  ['seq05', 1, 'seq03', 'seq04', 'seq05'],
                  ['seq05', 2, 'seq04', 'seq05',  np.nan],
                  ['seq06', 0, 'seq06',  np.nan,  np.nan],
                  ['seq06', 1, 'seq06',  np.nan,  np.nan],
                  ['seq07', 2, 'seq07',  np.nan,  np.nan],
                  ],
            columns=['db_seq', 'region', '0', '1', '2']
            )

        known = pd.Series({'seq00': 'seq00|seq01',
                           'seq01': 'seq00|seq01',
                           'seq03': 'seq03|seq04|seq05',
                           'seq04': 'seq03|seq04|seq05',
                           'seq05': 'seq03|seq04|seq05',
                           'seq06': 'seq06',
                           'seq07': 'seq07',
                           }, name='clean_name')
        known.index.set_names('db_seq', inplace=True)

        sub_map = _map_id_set(id_set=['seq00', 'seq01', 'seq03', 'seq04', 
                                      'seq05', 'seq06', 'seq07'],
                              tangle=tangled,
                              )
        pdt.assert_series_equal(known, sub_map.compute())

    def test_check_manifest_columns(self):
        manifest = Metadata(pd.DataFrame(
            data=np.array([['Bruce', 'Wayne'],
                          ['Dick', 'Grayson'],
                          ['Jason', 'Todd'],
                          ['Tim', 'Drake'],
                          ['Barbara', 'Gordon'],
                          ['Stephanie', 'Brown'],
                          ['Cassandra', 'Cain-Wayne'],
                          ['Damian', 'Wayne'],
                          ]),
            index=pd.Index(['Batman', 'Nightwing', 'Red Hood', 'Red Robin',
                            'Oracle', 'Batgirl', 'Black Bat', 'Robin'], 
                            name='id'),
            columns=['First Name', 'Last Name']
        ))
        # print(manifest)
        with self.assertRaises(ValidationError) as err:
            _check_manifest(manifest)
        self.assertEqual(str(err.exception), 
                         ('The manifest must contain the columns '
                          'kmer-map, alignment-map and frequency-table.\n'
                          'Please check the manifest and make sure all'
                          ' column names are spelled correctly')
                         )

    def test_check_manifest_missing(self):
        manifest = Metadata(pd.DataFrame(
            data=np.array([['Bruce', 'Wayne', None],
                          ['Dick', 'Grayson', '29'],
                          ['Jason', 'Todd', '24'],
                          ['Tim', 'Drake', '20'],
                          ['Barbara', 'Gordon', '32'],
                          ['Stephanie', 'Brown', '19'],
                          ['Cassandra', 'Cain-Wayne', '22'],
                          ['Damian', 'Wayne', '12'],
                          ]),
            index=pd.Index(['Batman', 'Nightwing', 'Red Hood', 'Red Robin',
                            'Oracle', 'Batgirl', 'Black Bat', 'Robin'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table']
        ))
        with self.assertRaises(ValidationError) as err:
            _check_manifest(manifest)
        self.assertEqual(str(err.exception), 
                         ('All regions must have a kmer-map, '
                          'alignment-map and frequency-table. Please '
                          'check and make sure that you have provided '
                          'all the files you need')
                         )

    def test_check_manifest_unique(self):
        manifest = Metadata(pd.DataFrame(
            data=np.array([['Bruce', 'Wayne', '45'],
                          ['Dick', 'Grayson', '29'],
                          ['Jason', 'Todd', '24'],
                          ['Tim', 'Drake', '20'],
                          ['Barbara', 'Gordon', '32'],
                          ['Stephanie', 'Brown', '19'],
                          ['Cassandra', 'Cain-Wayne', '22'],
                          ['Damian', 'Wayne', '12'],
                          ]),
            index=pd.Index(['Batman', 'Nightwing', 'Red Hood', 'Red Robin',
                            'Oracle', 'Batgirl', 'Black Bat', 'Robin'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table']
        ))
        with self.assertRaises(ValidationError) as err:
            _check_manifest(manifest)
        self.assertEqual(str(err.exception), 
                         ('All paths in the manifest must be unique.'
                          ' Please check your filepaths')
                         )

    def test_check_manifest_exists(self):
        manifest = Metadata(pd.DataFrame(
            data=np.array([['Bruce', 'Wayne', '45'],
                           ['Dick', 'Grayson', '29'],
                           ['Jason', 'Todd', '24'],
                           ]),
            index=pd.Index(['Batman', 'Nightwing', 'Red Hood'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table']
        ))
        with self.assertRaises(ValidationError) as err:
            _check_manifest(manifest)
        self.assertEqual(str(err.exception), 
                         ('All the paths in the manifest must exist.'
                          ' Please check your filepaths')
                         )

    def test_check_manifest_files(self):
        base_dir = self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        manifest = Metadata(pd.DataFrame(
            data=np.array([[self.base_dir, 
                           os.path.join(base_dir, 'full_db.qza'),
                           os.path.join(base_dir, 'region1_align.qza')]]),
            index=pd.Index(['test'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table']
            ))
        with self.assertRaises(ValidationError) as err:
            _check_manifest(manifest)
        self.assertEqual(str(err.exception), 
                         ('All the paths in the manifest must be files.'
                          ' Please check your filepaths')
                         )

    def test_read_manifest_files_path_error(self):
        base_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                                'files/little_test')
        manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(base_dir, 'region1_db_map.qza'),
                   os.path.join(base_dir, 'region1_align.qza'),
                   os.path.join(base_dir, 'region1_db_map.qza')],
                  [os.path.join(base_dir, 'region1_align.qza'),
                   os.path.join(base_dir, 'region1_align.qza'),
                   os.path.join(base_dir, 'region1_db_map.qza')],
                   ],
            columns=['kmer-map', 'alignment-map', 'frequency-table'],
            index=pd.Index(['Bludhaven', 'Gotham'], name='id')
        ))
        with self.assertRaises(TypeError) as err_:
            _read_manifest_files(manifest, 'kmer-map', 'FeatureData[KmerMap]')
        self.assertEqual(str(err_.exception),
            'Not all kmer map Artifacts are of the FeatureData[KmerMap] '
            'semantic type.\nPlease review semantic types for these'
            ' regions:\nGotham'
            )

    def test_read_manifest_artifact(self):
        base_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                                'files/little_test')
        manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(base_dir, 'region1_db_map.qza'),
                   os.path.join(base_dir, 'region1_align_map.qza'),
                   os.path.join(base_dir, 'region1_db_map.qza')],
                   ],
            columns=['kmer-map', 'alignment-map', 'frequency-table'],
            index=pd.Index(['Bludhaven'], name='id')
        ))
        test = _read_manifest_files(manifest, 'kmer-map', 
                                    'FeatureData[KmerMap]')
        self.assertEqual(len(test), 1)
        self.assertTrue(isinstance(test[0], Artifact))
        self.assertEqual(str(test[0].type), 'FeatureData[KmerMap]')
        pdt.assert_frame_equal(test[0].view(pd.DataFrame), 
                               ts.region1_db_map.view(pd.DataFrame))

    def test_read_manifest_view(self):
        base_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                                'files/little_test')
        manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(base_dir, 'region1_db_map.qza'),
                   os.path.join(base_dir, 'region1_align_map.qza'),
                   os.path.join(base_dir, 'region1_db_map.qza')],
                   ],
            columns=['kmer-map', 'alignment-map', 'frequency-table'],
            index=pd.Index(['Bludhaven'], name='id')
        ))
        test = _read_manifest_files(manifest, 'kmer-map', 
                                    'FeatureData[KmerMap]',
                                    pd.DataFrame)
        self.assertEqual(len(test), 1)
        self.assertTrue(isinstance(test[0], pd.DataFrame))
        pdt.assert_frame_equal(test[0], ts.region1_db_map.view(pd.DataFrame))

    def test_count_mapping_degen(self):
        known = pd.DataFrame.from_dict(orient='index', data={
            'seq05': {'num_regions': 2., 
                      'total_kmers_mapped': 3., 
                      'mean_kmer_per_region': 1.5,
                      'stdv_kmer_per_region': np.std([2, 3], ddof=1),
                      },
            'seq06': {'num_regions': 2, 
                      'total_kmers_mapped': 2,
                      'mean_kmer_per_region': 1,
                      'stdv_kmer_per_region': 0},
            })
        known.index.set_names('clean_name', inplace=True)
        test = _count_mapping(self.pair, count_degen=True)
        pdt.assert_frame_equal(test, known)

    def test_count_mapping_no_degen(self):
        known = pd.DataFrame.from_dict(orient='index', data={
            'seq05': {'num_regions': 2., 
                      'total_kmers_mapped': 2., 
                      'mean_kmer_per_region': 1.,
                      'stdv_kmer_per_region': 0.,
                      },
            'seq06': {'num_regions': 2, 
                      'total_kmers_mapped': 2,
                      'mean_kmer_per_region': 1,
                      'stdv_kmer_per_region': 0},
            })
        known.index.set_names('clean_name', inplace=True)
        test = _count_mapping(self.pair, count_degen=False)
        pdt.assert_frame_equal(test, known)

    def test_untangle_database_ids(self):
        matches = pd.DataFrame(
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
                           ['seq12', 'seq12', 'seq12|seq14|seq15', '1'],
                           ['seq13', 'seq13', 'seq13', '1'],
                           ['seq14', 'seq14', 'seq12|seq14|seq15', '1'],
                           ['seq15', 'seq15', 'seq12|seq14|seq15', '1'],
                           ['seq16', 'seq16', 'seq16|seq17',  '1'],
                           ['seq17', 'seq17', 'seq16|seq17', '1']
                           ], dtype=object),
            columns=['db_seq', 'seq_name', 'kmer', 'region'],
            )
        known_seq = pd.Series({'seq00': 'seq00', 
                                'seq01': 'seq01',
                                'seq02': 'seq02',
                                'seq03': 'seq03',
                                'seq04': 'seq04',
                                'seq05': 'seq05',
                                'seq06': 'seq06',
                                'seq07': 'seq07',
                                'seq08': 'seq08|seq09',
                                'seq09': 'seq08|seq09',
                                'seq10': 'seq10|seq11',
                                'seq11': 'seq10|seq11',
                                'seq12': 'seq12|seq14|seq15',
                                'seq13': 'seq13',
                                'seq14': 'seq12|seq14|seq15',
                                'seq15': 'seq12|seq14|seq15',
                                'seq16': 'seq16',
                                'seq17': 'seq17|seq18',
                                'seq18': 'seq17|seq18',
                                'seq19': 'seq19|seq20',
                                'seq20': 'seq19|seq20',
                                }, name='clean_name')
        known_seq.index.set_names('db_seq', inplace=True)
        # Generates the renaming
        seq_ = _untangle_database_ids(matches)
        pdt.assert_series_equal(known_seq, seq_[0].compute())

    def test_expand_duplicate_sequences(self):
        original = pd.DataFrame(data=[['1', '2|3', '3|4|5', '6'],
                                      [10, 20, 30, 40]],
                                index=['ids', 'counts']).T
        known = pd.DataFrame([[ '1', '2', '3', '3', '4', '5', '6'],
                              [10, 20, 20, 30, 30, 30, 40]],
                             index=['ids', 'counts']).T
        test = _expand_duplicate_sequences(original, id_col='ids')

        pdt.assert_frame_equal(test, known)

    def test_construct_align_mat(self):
        sequence_map = pd.Series({'seq1': 'seq1',
                                  'seq2': 'seq2',
                                  'seq3': 'seq3',
                                  'seq4': 'seq4',
                                  'seq5': 'seq5',
                                  'seq6': 'seq6'}, name='clean_name')
        sequence_map.index.set_names('db_seq', inplace=True)
        test_mat = _construct_align_mat(self.match1, 
                                        sequence_map,
                                        self.seq_summary).compute()

        # Rounds for the sake of sanity
        test_mat[self.float_cols] = test_mat[self.float_cols].round(4)
        test_mat[self.int_cols] = test_mat[self.int_cols].astype(int)
        pdt.assert_frame_equal(self.align1, test_mat)

    def test_solve_ml_em_iterative_1_sample(self):
        abund = sp.as_coo(np.array([0.18181818, 0.09090909, 0.09090909,
                                    0.09090909,  0.09090909, 0.09090909, 
                                    0.09090909, 0.09090909,  0.09090909, 
                                    0.09090909]))
        align = sp.as_coo(np.array([
            [0.4638, 0.4638, 0,      0,      0,      0     ],
            [0,      0,      0.4638, 0,      0,      0     ],
            [0,      0,      0,      0,      0.4638, 0.0008],
            [0,      0,      0,      0,      0.0008, 0.4638],
            [0.4638, 0,      0 ,     0,      0,      0     ],
            [0,      0.4638, 0,      0,      0,      0     ],
            [0,      0,      0.4638, 0,      0,      0     ],
            [0,      0,      0,      0.9276, 0,      0     ],
            [0,      0,      0,      0,      0.4638, 0     ],
            [0,      0,      0,      0,      0,      0.4638]]))

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


if __name__ == '__main__':
    main()