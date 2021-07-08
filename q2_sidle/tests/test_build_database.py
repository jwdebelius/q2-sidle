from unittest import TestCase, main

import copy
import os
import warnings

import biom
import dask
from dask.delayed import Delayed
import dask.dataframe as dd
import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import skbio

from qiime2 import Artifact, Metadata
from qiime2.plugin import ValidationError

from q2_sidle._build_database import (reconstruct_database,
                                      _build_region_db,
                                      _check_db_list,
                                      _check_intersection_delayed,
                                      _check_regions,
                                      _clean_kmer_list,
                                      _define_shared,
                                      _detangle_names,
                                      _filter_to_aligned,
                                      _get_regional_seqs,
                                      _get_unique_kmers,
                                      _pull_unique,
                                      )

import q2_sidle.tests.test_set as ts


class DatabaseMapTest(TestCase):
    def setUp(self):
        self.kmers_delayed = [ts.region1_db_map.view(Delayed),
                              ts.region2_db_map.view(Delayed)]
        self.kmers_df = [ts.region1_db_map.view(pd.DataFrame),
                         ts.region2_db_map.view(pd.DataFrame)]
        self.regions = ['Bludhaven', 'Gotham']
        self.reional_alignment = [ts.region1_align.view(pd.DataFrame),
                                  ts.region2_align.view(pd.DataFrame)]
        self.match0 = pd.DataFrame(
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
                           ['seq21', 'seq21', 'seq21', '0']]),
            columns=['db-seq', 'seq-name', 'kmer', 'region']
            )
        self.match1 = pd.DataFrame(
            data=np.array([['seq00', 'seq00', 'seq00', '1'],
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
            )
        self.match2 = pd.DataFrame(
            data=np.array([['seq22', 'seq22', 'seq22', '2']]),
            columns=['db-seq', 'seq-name', 'kmer', 'region'],
            )

        self.untangled = pd.Series({'seq00': 'seq00', 
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
        self.untangled.index.set_names('db-seq', inplace=True)

    def reconstruct_database_parallel(self):
        test = reconstruct_database(
            region=['0', '1', '2'],
            regional_alignment=[self.match0, self.match1, self.match2],
            kmer_map=[dask.delayed(self.match0), 
                      dask.delayed(self.match1), 
                      dask.delayed(self.match2)],
            debug=True
            )
        pdt.assert_frame_equal(self.tangled, test)

    def test_build_region_db(self):
        kmers = ['seq1', 'seq2', 'seq3']
        known = pd.DataFrame(
            data=np.array([['Bludhaven', 'seq1|seq2'],
                           ['Bludhaven', 'seq1|seq2'],
                           ['Bludhaven', 'seq3@0001'],
                           ['Bludhaven', 'seq3@0002'],
                           ]),
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq3'], 
                           name='db-seq'),
            columns=['region', 'kmer']
            ).reset_index()
        test = _build_region_db(self.kmers_df[0].reset_index(), kmers=kmers)
        pdt.assert_frame_equal(known, test)

    def test_check_db_list(self):
        ref_seqs = {'Nightwing', 'RedHood', 'RedRobin', 'Spoiler'}
        kmer = 'Nightwing|RedHood|TimDrake|DamianWayne'
        known = 'Nightwing|RedHood|X'

        test = _check_db_list(kmer, ref_seqs)
        self.assertEqual(known, test)

    def test_check_intersection_delayed(self):
        df_ = pd.DataFrame(
            data=np.array([['seq00', 'seq00', '0'],
                           ['seq01', 'seq01|seq02', '0'],
                           ['seq02', 'seq01|seq02', '0'],
                           ['seq05', 'seq05', '0'],
                           ['seq06', 'seq06', '0'],
                           ['seq07', 'seq07', '0'],
                           ['seq08', 'seq08|seq09', '0'],
                           ["seq09", 'seq08|seq09', '0'],
                           ['seq10', 'seq10|seq11', '0'],
                           ['seq11', 'seq10|seq11', '0'],
                           ['seq12', 'seq12|seq13|seq14', '0'],
                           ['seq13', 'seq12|seq13|seq14', '0'],
                           ['seq14', 'seq12|seq13|seq14', '0'],
                           ['seq16', 'seq16', '0'],
                           ['seq17', 'seq17|seq18', '0'],
                           ['seq18', 'seq17|seq18', '0'],
                           ['seq19', 'seq19|seq20', '0'],
                           ['seq20', 'seq19|seq20', '0'],
                           ['seq21', 'seq21', '0'],
                           ['seq00', 'seq00', '1'],
                           ['seq01', 'seq01', '1'],
                           ['seq02', 'seq02', '1'],
                           ['seq05', 'seq05|seq06', '1'],
                           ['seq06', 'seq05|seq06', '1'],
                           ['seq09', 'seq09', '1'],
                           ['seq10', 'seq10|seq11', '1'],
                           ['seq11', 'seq10|seq11', '1'],
                           ['seq12', 'seq12|seq14|seq15', '1'],
                           ['seq13', 'seq13|seq14|seq15', '1'],
                           ['seq14', 'seq12|seq13|seq14|seq15', '1'],
                           ['seq15', 'seq12|seq13|seq14|seq15', '1'],
                           ['seq16', 'seq16|seq17',  '1'],
                           ['seq17', 'seq16|seq17', '1'],
                           ['seq21', 'seq21|seq22|seq23', '1'],
                           ['seq22', 'seq21|seq22|seq23', '1'],
                           ['seq23', 'seq21|seq22|seq23', '1'],
                           ['seq22', 'seq22', '2'],
                           ], dtype=object),
            columns=['db-seq', 'kmer', 'region'],
            )
        df_['region'] =  df_['region'].astype(int)
        df_ = dd.from_pandas(df_, npartitions=1)
        known_kmers = pd.DataFrame(
            data=[['seq00', 'seq00'],
                  ['seq01', 'seq01'],
                  ['seq02', 'seq02'],
                  ['seq05', 'seq05'],
                  ['seq06', 'seq06'],
                  ['seq07', 'seq07'],
                  ['seq08', 'seq08|seq09'],
                  ['seq09', 'seq09'],
                  ['seq10', 'seq10|seq11'],
                  ['seq11', 'seq10|seq11'],
                  ['seq12', 'seq12|seq14'],
                  ['seq13', 'seq13|seq14'],
                  ['seq14', 'seq12|seq13|seq14'],
                  ['seq15', 'seq12|seq13|seq14|seq15'],
                  ['seq16', 'seq16'],
                  ['seq17', 'seq17'],
                  ['seq18', 'seq17|seq18'],
                  ['seq19', 'seq19|seq20'],
                  ['seq20', 'seq19|seq20'],
                  ['seq21', 'seq21'],
                  ['seq22', 'seq22'],
                  ['seq23', 'seq21|seq22|seq23'],
                  ],
            columns=['db-seq', 'kmer'],
            )
        known_check = pd.DataFrame(
            data=[['seq00', 'seq00', True],
                  ['seq01', 'seq01', True],
                  ['seq02', 'seq02', True],
                  ['seq05', 'seq05', True],
                  ['seq06', 'seq06', True],
                  ['seq07', 'seq07', True],
                  ['seq08|seq09', 'seq08', False],
                  ['seq09', 'seq09', True],
                  ['seq10|seq11', 'seq10|seq11', True],
                  ['seq12|seq13|seq14', 'seq14', False],
                  ['seq12|seq13|seq14|seq15', 'seq15', False],
                  ['seq12|seq14', 'seq12', False],
                  ['seq13|seq14', 'seq13', False],
                  ['seq16', 'seq16', True],
                  ['seq17', 'seq17', True],
                  ['seq17|seq18', 'seq18', False],
                  ['seq19|seq20', 'seq19|seq20', True],
                  ['seq21', 'seq21', True],
                  ['seq21|seq22|seq23', 'seq23', False],
                  ['seq22', 'seq22', True],
                  ],
            columns=['kmer', 'shared-name', 'tidy']
            )

        test_kmers, test_check = _check_intersection_delayed(df_)
        # print(test_check.compute())
        pdt.assert_frame_equal(known_kmers, test_kmers.compute())
        pdt.assert_frame_equal(known_check, test_check.compute())

    def test_check_regions(self):
        regions = ['Bludhaven', 'Gotham']
        k_order = {'Bludhaven': 0, 'Gotham': 1}
        k_names = {0: 'Bludhaven', 1: "Gotham"}

        t_order, t_names, t_num = _check_regions(self.regions)
        self.assertEqual(t_num, 2)
        npt.assert_array_equal(list(k_order.keys()), list(t_order.keys()))
        for k, v in k_order.items():
            self.assertEqual(v, t_order[k])
        npt.assert_array_equal(list(k_names.keys()), list(t_names.keys()))
        for k, v in k_names.items():
            self.assertEqual(v, t_names[k])

    def test_clean_kmer_list(self):
        test = _clean_kmer_list('robin@0001|robin@002|batman@001')
        self.assertEqual(test,  'batman|robin')

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

    def test_filter_to_aligned(self):
        kmers = ['seq1', 'seq2', 'seq3']
        known = pd.DataFrame(
            data=np.array([['Bludhaven', 'seq1|seq2'],
                           ['Bludhaven', 'seq1|seq2'],
                           ['Bludhaven', 'seq3@0001'],
                           ['Bludhaven', 'seq3@0002'],
                           ]),
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq3'], 
                           name='db-seq'),
            columns=['region', 'kmer']
            ).reset_index()
        test = _filter_to_aligned(self.kmers_df[0], kmers)
        pdt.assert_frame_equal(known, test)

    def test_get_unique_kmers(self):
        test = pd.Series(['seq00', 'seq00|seq01@001', 'seq01@002', 
                          'seq01|seq02'])
        known = np.array(['seq00', 'seq01', 'seq02'])
        test = _get_unique_kmers(test)
        npt.assert_array_equal(test, known)

    def test_get_regional_seqs(self):
        x = pd.Series(['robin|redhood', 
                       'robin|spoiler|batgirl',
                       'robin',
                       ])
        known = 'batgirl|redhood|robin|spoiler'
        test = _get_regional_seqs(x)
        self.assertEqual(test, known)

    def test_pull_unique(self):
        known = np.array(['seq1', 'seq2', 'seq3', 'seq5', 'seq6'])
        test = _pull_unique(self.kmers_delayed[0]).compute()
        npt.assert_array_equal(test, known)

    def test_define_shared(self):
        known = 'robin'
        g = pd.DataFrame(data=np.array([['robin|nightwing|batman',
                                         'robin|redhood', 
                                         'robin|spoiler|batgirl',
                                         'robin']]).T,
                         columns=['kmer']
                         )
        matched = {"nightwing", 'batman'}
        test = _define_shared(g, matched)
        self.assertEqual(test, known)

if __name__ == '__main__':
    main()