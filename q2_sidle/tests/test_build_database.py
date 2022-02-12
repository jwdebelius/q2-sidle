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
                                      _clean_kmer_list,
                                      _count_mapping,
                                      _define_shared,
                                      _detangle_names,
                                      _filter_to_aligned,
                                      _filter_to_defined,
                                      _get_regional_seqs,
                                      _get_unique_kmers,
                                      _map_aligned_asvs,
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
            data=np.array([['seq00', 'seq00', 'seq00', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq01', 'seq01', 'seq01|seq02', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq02', 'seq02', 'seq01|seq02', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq03', 'seq03@001', 'seq03@001', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq03', 'seq03@002', 'seq03@002', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq04', 'seq04@001', 'seq04@001', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq04', 'seq04@002', 'seq04@002', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq05', 'seq05', 'seq05', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq06', 'seq06', 'seq06', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq07', 'seq07', 'seq07', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq08', 'seq08', 'seq08|seq09', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ["seq09", 'seq09', 'seq08|seq09', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq10', 'seq10', 'seq10|seq11', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq11', 'seq11', 'seq10|seq11', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq12', 'seq12', 'seq12|seq13|seq14', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq13', 'seq13', 'seq12|seq13|seq14', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq14', 'seq14', 'seq12|seq13|seq14', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq16', 'seq16', 'seq16', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq17', 'seq17', 'seq17|seq18', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq18', 'seq18', 'seq17|seq18', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq19', 'seq19', 'seq19|seq20', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq20', 'seq20', 'seq19|seq20', '0', 'WANTCAT', 'ATGATGATG', np.nan, np.nan, 50],
                           ['seq21', 'seq21', 'seq21', '0', 'WANTCAT','ATGATGATG', np.nan, np.nan, 50],
                           ], dtype=object),
            columns=['db-seq', 'seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer', 'fwd-pos', 'rev-pos', 'kmer-length'],
            )
        self.match1 = pd.DataFrame(
            data=np.array([['seq00', 'seq00', 'seq00', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq01', 'seq01', 'seq01', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq02', 'seq02', 'seq02', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq03', 'seq03', 'seq03', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq04', 'seq04@001', 'seq04@001', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq04', 'seq04@002', 'seq04@002', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq04', 'seq04@003', 'seq04@003', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq05', 'seq05@001', 'seq05@001', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq05', 'seq05@002', 'seq05@002|seq06', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq06', 'seq06', 'seq05@002|seq06', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq09', 'seq09', 'seq09', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq10', 'seq10', 'seq10|seq11', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq11', 'seq11', 'seq10|seq11', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq12', 'seq12@0001', 'seq12@0001', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq12', 'seq12@0002', 'seq12@0002|seq14|seq15', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq13', 'seq13@0001', 'seq13@0001|seq14|seq15', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq13', 'seq13@0002', 'seq13@002', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq14', 'seq14', 'seq12@0002|seq13@0001|seq14|seq15', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq15', 'seq15', 'seq12@0002|seq13@0001|seq14|seq15', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq16', 'seq16', 'seq16|seq17',  '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq17', 'seq17', 'seq16|seq17', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq21', 'seq21', 'seq21|seq22|seq23', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq22', 'seq22', 'seq21|seq22|seq23', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ['seq23', 'seq23', 'seq21|seq22|seq23', '1', 'CACCTCGTN', 'MTGACGTG', np.nan, np.nan, 75],
                           ], dtype=object),
             columns=['db-seq', 'seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer', 'fwd-pos', 'rev-pos', 'kmer-length'],
            )
        self.match2 = pd.DataFrame(
            data=np.array([['seq22', 'seq22', 'seq22', '2', 'AACACA', 'CATAAG', np.nan, np.nan, 25]], dtype=object),
            columns=['db-seq', 'seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer', 'fwd-pos', 'rev-pos', 'kmer-length'],
            )
        for match in [self.match0, self.match1, self.match2]:
            match.set_index('db-seq', inplace=True)
            match['kmer-length'] = match['kmer-length'].astype(int)
            match[['fwd-pos', 'rev-pos']] = match[['fwd-pos', 'rev-pos']].astype(float)


        self.align0 = pd.DataFrame(
            data=np.array([['seq00', 'asv_aa', 50, 0, 2, '0'],
                           ['seq01|seq02', 'asv_ab', 50, 0, 2, '0'],
                           ['seq03@002', 'asv_ac', 50, 0, 2, '0'],
                           ['seq04@001', 'asv_ad', 50, 0, 2, '0'],
                           ['seq04@002', 'asv_ae', 50, 0, 2, '0'],
                           ['seq05', 'asv_af', 50, 0, 2, '0'],
                           ['seq06', 'asv_ag', 50, 0, 2, '0'],
                           ['seq07', 'asv_ah', 50, 0, 2, '0'],
                           ['seq08|seq09', 'asv_ai', 50, 0, 2, '0'],
                           ['seq10|seq11', 'asv_aj', 50, 0, 2, '0'],
                           ['seq10|seq11', 'asv_ak', 50, 0, 2, '0'],
                           ['seq12|seq13|seq14', 'asv_al', 50, 0, 2, '0'],
                           ['seq16', 'asv_ao', 50, 0, 2, '0'],
                           ['seq17|seq18', 'asv_ap', 50, 0, 2, '0'],
                           ['seq19|seq20', 'asv_aq', 50, 0, 2, '0'],
                           ['seq21', 'asv_ar', 50, 0, 2, '0'],
                           ], dtype=object),
            columns=['kmer', 'asv', 'length', 'mismatch', 'max-mismatch', 
                     'region'],
            )
        self.align1 = pd.DataFrame(
            data=[['seq00', 'asv_ba', 75, 0, 3, '1'],
                  ['seq01', 'asv_bb', 75, 0, 3, '1'],
                  ['seq02', 'asv_bc', 75, 0, 3, '1'],
                  ['seq03', 'asv_bd', 75, 0, 3, '1'],
                  ['seq04@003', 'asv_be', 75, 0, 3, '1'],
                  ['seq05@001', 'asv_bf', 75, 0, 3, '1'],
                  ['seq05@002|seq06', 'asv_bg', 75, 0, 3, '1'],
                  ['seq09', 'asv_bi', 75, 0, 3, '1'],
                  ['seq10|seq11', 'asv_bj', 75, 0, 3, '1'],
                  ['seq12@0002|seq14|seq15', 'asv_bl', 75, 0, 3, '1'],
                  ['seq13@002', 'asv_bm', 75, 0, 3, '1'],
                  ['seq12@0002|seq13@0001|seq14|seq15', 'asv_bn', 
                   75, 0, 3, '1'],
                  ['seq16|seq17', 'asv_bo', 75, 0, 3, '1'],
                  ['seq21|seq22|seq23', 'asv_bq', 75, 0, 3, '1'],
                  ],
            columns=['kmer', 'asv', 'length', 'mismatch', 'max-mismatch', 
                     'region'],
            )
        self.align2 = pd.DataFrame(
            data=np.array([['seq23', 'asv_cq', 25, 0, 1, '2']], dtype=object),
            columns=['kmer', 'asv', 'length', 'mismatch', 'max-mismatch',
                      'region'],
            )
        
        self.untangled = pd.DataFrame.from_dict(orient='index', data={
            'seq00': {'clean_name': 'seq00', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq01': {'clean_name': 'seq01', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq02': {'clean_name': 'seq02', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq03': {'clean_name': 'seq03', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq04': {'clean_name': 'seq04', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq05': {'clean_name': 'seq05', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq06': {'clean_name': 'seq06', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq07': {'clean_name': 'seq07', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 0, 'last-fwd-primer': 'WANTCAT', 'last-fwd-pos': np.nan, 'last-kmer-length': 50},
            'seq08': {'clean_name': 'seq08', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 0, 'last-fwd-primer': 'WANTCAT', 'last-fwd-pos': np.nan, 'last-kmer-length': 50},
            'seq09': {'clean_name': 'seq09', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq10': {'clean_name': 'seq10|seq11','first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq11': {'clean_name': 'seq10|seq11', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq12': {'clean_name': 'seq12', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq13': {'clean_name': 'seq13', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq14': {'clean_name': 'seq14', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq15': {'clean_name': 'seq15', 'first-region': 1, 'first-fwd-primer': 'CACCTCGTN', 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq16': {'clean_name': 'seq16', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq17': {'clean_name': 'seq17', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq18': {'clean_name': 'seq18', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 0, 'last-fwd-primer': 'WANTCAT', 'last-fwd-pos': np.nan, 'last-kmer-length': 50},
            'seq19': {'clean_name': 'seq19|seq20', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 0, 'last-fwd-primer': 'WANTCAT', 'last-fwd-pos': np.nan, 'last-kmer-length': 50},
            'seq20': {'clean_name': 'seq19|seq20', 'first-region': 0, 'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 0, 'last-fwd-primer': 'WANTCAT', 'last-fwd-pos': np.nan, 'last-kmer-length': 50},
            'seq21': {'clean_name': 'seq21', 'first-region': 0,  'first-fwd-primer': 'WANTCAT', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            'seq22': {'clean_name': 'seq22', 'first-region': 1, 'first-fwd-primer': 'CACCTCGTN', 'first-fwd-pos': np.nan, 'last-region': 2, 'last-fwd-primer': 'AACACA', 'last-fwd-pos': np.nan, 'last-kmer-length': 25},
            'seq23': {'clean_name': 'seq23', 'first-region': 1, 'first-fwd-primer': 'CACCTCGTN', 'first-fwd-pos': np.nan, 'last-region': 1, 'last-fwd-primer': 'CACCTCGTN', 'last-fwd-pos': np.nan, 'last-kmer-length': 75},
            })
        self.untangled.index.set_names('db-seq', inplace=True)

        self.summary = pd.DataFrame.from_dict(orient='index', data={
            'seq00': {'num-regions': 2., 'total-kmers-mapped': 2., 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_aa|asv_ba'},
            'seq01': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ab|asv_bb'},
            'seq02': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ab|asv_bc'},
            'seq03': {'num-regions': 2, 'total-kmers-mapped': 3, 'mean-kmer-per-region': 1.5, 'stdv-kmer-per-region': np.std([1, 2], ddof=1), 'mapped-asvs': 'asv_ac|asv_bd'},
            'seq04': {'num-regions': 2, 'total-kmers-mapped': 5, 'mean-kmer-per-region': 2.5, 'stdv-kmer-per-region': np.std([2, 3], ddof=1), 'mapped-asvs': 'asv_ad|asv_ae|asv_be'},
            'seq05': {'num-regions': 2, 'total-kmers-mapped': 3, 'mean-kmer-per-region': 1.5, 'stdv-kmer-per-region': np.std([1, 2], ddof=1), 'mapped-asvs': 'asv_af|asv_bf|asv_bg'},
            'seq06': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ag|asv_bg'},
            'seq07': {'num-regions': 1, 'total-kmers-mapped': 1, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ah'},
            'seq08': {'num-regions': 1, 'total-kmers-mapped': 1, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ai'},
            'seq09': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ai|asv_bi'},
            'seq10|seq11': {'num-regions': 2, 'total-kmers-mapped': 4, 'mean-kmer-per-region': 2, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_aj|asv_ak|asv_bj'},
            'seq12': {'num-regions': 2, 'total-kmers-mapped': 3, 'mean-kmer-per-region': 1.5, 'stdv-kmer-per-region': np.std([1, 2], ddof=1), 'mapped-asvs': 'asv_al|asv_bl|asv_bn'},
            'seq13': {'num-regions': 2, 'total-kmers-mapped': 3, 'mean-kmer-per-region': 1.5, 'stdv-kmer-per-region': np.std([1, 2], ddof=1), 'mapped-asvs': 'asv_al|asv_bm|asv_bn'},
            'seq14': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_al|asv_bl|asv_bn'},
            'seq15': {'num-regions': 1, 'total-kmers-mapped': 1, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_bl|asv_bn'},
            'seq16': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ao|asv_bo'},
            'seq17': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ap|asv_bo'},
            'seq18': {'num-regions': 1, 'total-kmers-mapped': 1, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ap'},
            'seq19|seq20': {'num-regions': 1, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 2, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_aq'},
            'seq21': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_ar|asv_bq'},
            'seq22': {'num-regions': 2, 'total-kmers-mapped': 2, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_bq'},
            'seq23': {'num-regions': 1, 'total-kmers-mapped': 1, 'mean-kmer-per-region': 1, 'stdv-kmer-per-region': 0, 'mapped-asvs': 'asv_bq|asv_cq'},
            })
        self.summary.index.set_names('feature-id', inplace=True)
        self.pair = pd.DataFrame(
            data=[['seq05', 0, 'seq05', 'seq05'],
                  ['seq05@001', 1, 'seq05', 'seq05'],
                  ['seq05@002', 1, 'seq05', 'seq05'],
                  ['seq06', 0, 'seq06', 'seq06'],
                  ['seq06', 1, 'seq06', 'seq06']
                  ],
            columns=['kmer', 'region', 'db-seq', 'clean_name'],
            )

    def test_reconstruct_database_parallel(self):
        t_map, t_summary = reconstruct_database(
            region=['0', '1', '2'],
            regional_alignment=[self.align0, self.align1, self.align2],
            kmer_map=[dask.delayed(self.match0), 
                      dask.delayed(self.match1), 
                      dask.delayed(self.match2)],
            debug=True
            )
        t_map.sort_index(inplace=True)
        pdt.assert_frame_equal(self.untangled, t_map)
        pdt.assert_frame_equal(self.summary, t_summary.to_dataframe())

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
                           ['seq16', 'seq16|seq17', '1'],
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

    def test_clean_kmer_list(self):
        test = _clean_kmer_list('robin@0001|robin@002|batman@001')
        self.assertEqual(test, 'batman|robin')

    # def test_count_mapping_degen(self):
    #     known = pd.DataFrame.from_dict(orient='index', data={
    #         'seq05': {'num-regions': 2., 
    #                   'total-kmers-mapped': 3., 
    #                   'mean-kmer-per-region': 1.5,
    #                   'stdv-kmer-per-region': np.std([2, 3], ddof=1),
    #                   },
    #         'seq06': {'num-regions': 2, 
    #                   'total-kmers-mapped': 2,
    #                   'mean-kmer-per-region': 1,
    #                   'stdv-kmer-per-region': 0},
    #         })
    #     known.index.set_names('clean_name', inplace=True)
    #     test = _count_mapping(self.pair, count_degen=True)
    #     pdt.assert_frame_equal(test, known)

    # def test_count_mapping_no_degen(self):
    #     known = pd.DataFrame.from_dict(orient='index', data={
    #         'seq05': {'num-regions': 2., 
    #                   'total-kmers-mapped': 2., 
    #                   'mean-kmer-per-region': 1.,
    #                   'stdv-kmer-per-region': 0.,
    #                   },
    #         'seq06': {'num-regions': 2, 
    #                   'total-kmers-mapped': 2,
    #                   'mean-kmer-per-region': 1,
    #                   'stdv-kmer-per-region': 0},
    #         })
    #     known.index.set_names('clean_name', inplace=True)
    #     test = _count_mapping(self.pair, count_degen=False)
    #     pdt.assert_frame_equal(test, known)

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
                  'seq04', 'seq05', 'seq06', 'seq07|seq08', 
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

    def test_filter_to_aligned_with_kmers(self):
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

    def test_filter_to_aligned_no_kmers(self):
        known = pd.DataFrame(
            data=np.array([['Bludhaven', 'seq1|seq2'],
                           ['Bludhaven', 'seq1|seq2'],
                           ['Bludhaven', 'seq3@0001'],
                           ['Bludhaven', 'seq3@0002'],
                           ['Bludhaven', 'seq5'],
                           ['Bludhaven', 'seq6'],
                           ]),
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq3', 'seq5', 'seq6'], 
                           name='db-seq'),
            columns=['region', 'kmer']
            ).reset_index()
        test = _filter_to_aligned(self.kmers_df[0])
        pdt.assert_frame_equal(known, test)

    def filter_to_defined(self):
        defined = ['seq1', 'seq2', 'seq4']
        known = pd.DataFrame.from_dict(orient='index', data={
            'seq1': {'region': 'Bludhaven', 
                     'fwd-primer': 'WANTCAT', 
                     'fwd-pos': 13,
                     'kmer-length': 15},
            'seq2': {'region': 'Bludhaven', 
                     'fwd-primer': 'WANTCAT', 
                     'fwd-pos': 13,
                     'kmer-length': 15},
            })
        test = _filter_to_defined(self.kmers_df[0], defined)
        pd.assert_frame_equal(known, test)

    def test_get_unique_kmers(self):
        test = pd.Series(['seq00', 'seq00|seq01@001', 'seq01@002', 
                          'seq01|seq02'])
        known = np.array(['seq00', 'seq01', 'seq02'])
        test = _get_unique_kmers(test)
        npt.assert_array_equal(test, known)

    # def test_get_regional_seqs(self):
    #     x = pd.Series(['robin|redhood', 
    #                    'robin|spoiler|batgirl',
    #                    'robin',
    #                    ])
    #     known = 'batgirl|redhood|robin|spoiler'
    #     test = _get_regional_seqs(x)
    #     self.assertEqual(test, known)

    # def test_map_aligned_kmers(self):
    #     seq_map = self.untangled['clean_name'].copy()
    #     known = pd.DataFrame(
    #        np.array([['seq00', 'asv_aa'],
    #                  ['seq01', 'asv_ab'],
    #                  ['seq02', 'asv_ab'],
    #                  ['seq03', 'asv_ac'],
    #                  ['seq04', 'asv_ad'],
    #                  ['seq04', 'asv_ae'],
    #                  ['seq05', 'asv_af'],
    #                  ['seq06', 'asv_ag'],
    #                  ['seq07', 'asv_ah'],
    #                  ['seq08', 'asv_ai'],
    #                  ['seq09', 'asv_ai'],
    #                  ['seq10|seq11', 'asv_aj'],
    #                  ['seq10|seq11', 'asv_ak'],
    #                  ['seq12', 'asv_al'],
    #                  ['seq13', 'asv_al'],
    #                  ['seq14', 'asv_al'],
    #                  ['seq16', 'asv_ao'],
    #                  ['seq17', 'asv_ap'],
    #                  ['seq18', 'asv_ap'],
    #                  ['seq19|seq20', 'asv_aq'],
    #                  ['seq21', 'asv_ar']], dtype=str),
    #        columns=['clean_name', 'asv']
    #        )
    #     known = known.set_index('clean_name')['asv']
    #     test = _map_aligned_asvs(self.align0, self.untangled['clean_name'])
    #     test.sort_index(inplace=True)
    #     pdt.assert_series_equal(known, test)

    def test_pull_unique(self):
        known = np.array(['seq1', 'seq2', 'seq3', 'seq5', 'seq6'])
        test = _pull_unique(self.kmers_delayed[0]).compute()
        npt.assert_array_equal(test, known)



if __name__ == '__main__':
    main()