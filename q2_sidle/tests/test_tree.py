from unittest import TestCase, main

import os

import dask
import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact, Metadata
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._tree import (reconstruct_fragment_rep_seqs,
                            _collapse_expanded,
                            _expand_primer,
                            _find_positions,
                            _find_multiple_region_fragments,
                            _find_single_region_fragments,
                            )
from q2_sidle.tests import test_set as ts

class TreeTest(TestCase):
    def setUp(self):
        self.base_dir = \
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                         'files/little_test')
        self.single_fragment = pd.DataFrame(
            data=[['seq01|seq02', 'seq01', 0, '[AT]A[ACGT]TCAT', 'ATGATGATG', 
                   15, 1, '-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC'],
                  ['seq01|seq02', 'seq02', 0, '[AT]A[ACGT]TCAT', 'ATGATGATG',
                   15, 1, 'ACTAGTCATGCGAAGCGGCTCAGGATGATGATGAA---']],
            columns=['clean_name', 'db-seq', 'region', 'fwd-primer', 
                     'rev-primer', 'kmer-length', 'final_coverage', 
                     'sequence'],
                )
        self.multiple_fragments = pd.DataFrame(
            data=[['seq03|seq04', 'seq03', 0, '[AT]A[ACGT]TCAT', 'ATGATGATG',
                   15, 2, ('CATAGTCATWTCCGCGTTGGAGTTATGATGATGAWACCACCTCGTCCCA'
                   'GTTCCGCGCTTCTGACGTGCA-'), 9, 60, 9, 60],
                  ['seq03|seq04', 'seq03', 1, 'CACCTCGT[ACGT]', '[AC]TGACGTG',
                   15, 2, ('CATAGTCATWTCCGCGTTGGAGTTATGATGATGAWACCACCTCGTCCC'
                    'AGTTCCGCGCTTCTGACGTGCA-'), 9, 60, 9, 60],
                  ['seq03|seq04', 'seq04', 1, 'CACCTCGT[ACGT]', '[AC]TGACGTG',
                   15, 2, ('------------------GGAGTTATGATGA--AGACCACCTCGTCCC'
                    'AGTTCCGCGCTTCTGACGTGCAC'), 46, 60, 9, 60]],
            columns=['clean_name', 'db-seq', 'region', 'fwd-primer', 
                     'rev-primer', 'kmer-length', 'final_coverage', 
                     'sequence', 'fwd-pos', 'rev-pos', 'group-fwd', 
                     'group-rev'],
                )
        self.multiple_fragments_trim = pd.Series(
            data=['WTCCGCGTTGGAGTTATGATGATGAWACCACCTCGTCCCAGTTCCGCGCTT'], 
            name='fragment', 
            index=pd.Index(['seq03|seq04'], name='clean_name')
            )

    def test_reconstruct_fragment_rep_seqs(self):
        recon_map = pd.Series(
            data=['seq01|seq02', 'seq01|seq02', 'seq03|seq04', 
                  'seq03|seq04', 'seq05'],
            index=pd.Index(['seq01', 'seq02', 'seq03', 'seq04', 'seq05'], 
                            name='db-seq'),
            name='clean_name'
            )
        recon_summary = pd.DataFrame(
            data=[[1, 2, 2, 0, 'asv01|asv02'],
                  [2, 3, 1.5, np.std([1, 2], ddof=1), 'asv03|asv04'],
                  [2, 2, 1, 0, 'asv07|asv08']],
            index=pd.Index(['seq01|seq02', 'seq03|seq04', 'seq05'], 
                            name='clean_name'),
            columns=['num-regions', 'total-kmers-mapped', 
                     'mean-kmer-per-region', 'stdv-kmer-per-region', 
                     'mapped-asvs'],
            )
        manifest = Metadata(pd.DataFrame(
            data=[[os.path.join(self.base_dir, 'frag_r1_db_map.qza'),
                   os.path.join(self.base_dir, 'region1_align.qza'),
                   os.path.join(self.base_dir, 'region1_counts.qza'), 0],
                  [os.path.join(self.base_dir, 'frag_r2_db_map.qza'),
                   os.path.join(self.base_dir, 'region2_align.qza'),
                   os.path.join(self.base_dir, 'region2_counts.qza'), 1]],
            columns=['kmer-map', 'alignment-map', 'frequency-table', 
                     'region-order'],
            index=pd.Index(['Gotham', 'Bludhaven'], name='id')
        ))
        aligned_seqs = pd.Series({
            'seq01': '-CTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC---------------------------------',
            'seq02': 'ACTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGAC---------------------------------',
            'seq03': 'CATAGTCATWTCCGCGTTGGAGTTATGATGATGAWACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCA-',
            'seq04': '------------------GGAGTTATGATGA--AGACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCAC',
            'seq05': 'CATAGTCATCGTTTATGTATGCCCATGATGATGCGAGCACCTCGTATGGATGTAGAGCCACTGACGTGCGG',
        })
        known = pd.Series(
            data=['GCGAAGCGGCTCAGG',
                  'WTCCGCGTTGGAGTTATGATGATGAWACCACCTCGTCCCAGTTCCGCGCTT'],
            index=pd.Index(['seq01|seq02', 'seq03|seq04']),
            )
        test = reconstruct_fragment_rep_seqs(
            reconstruction_map=recon_map,
            reconstruction_summary=recon_summary,
            manifest=manifest,
            aligned_sequences=aligned_seqs,
            trim_to_fragment=True,
            # debug=True,
            )
        pdt.assert_series_equal(test.view(pd.Series).astype(str), known)

    def test_collapse_expanded_drop_no_trim(self):
        self.multiple_fragments.drop_duplicates(['db-seq'], inplace=True)
        known = pd.Series({
            "seq03|seq04": 'WTCCGCGTTGGAGTTATGATGAADACCACCTCGTCCCAGTTCCGCGCTT'
            })
        test =  _collapse_expanded('seq03|seq04',
                                    self.multiple_fragments.copy(), 
                                    'drop', 
                                    trim=False)
        pdt.assert_series_equal(known, test)

    def test_collapse_expanded_keep_no_trim(self):
        self.multiple_fragments.drop_duplicates(['db-seq'], inplace=True)
        known = pd.Series({
            "seq03|seq04": 
            'WTCCGCGTTGGAGTTATGATGATGADACCACCTCGTCCCAGTTCCGCGCTT'
            })
        test = _collapse_expanded('seq03|seq04', 
                                 self.multiple_fragments.copy(), 
                                 'keep', 
                                 trim=False)
        pdt.assert_series_equal(known, test)

    def _collapse_expanded_trim_drop(self):
        self.multiple_fragments.drop_duplicates(['db-seq'], inplace=True)
        test = _collapse_expanded_trim('seq03|seq04', 
                                       self.multiple_fragments.copy(),
                                       gap_management='drop')
        pdt.assert_series_equal(test, self.multiple_fragments_trim)
    
    def test_expand_primer(self):
        primer = 'WANTCAT'
        known = '((?<=([AT])))((A[ACGT]TCAT){e<=1})'
        test = _expand_primer(primer, 1)
        self.assertEqual(known, test)

    def test_find_positions(self):
        seqs = pd.DataFrame(
            data=[['TGTACTCATGGACAATGTAATGATGATGGAT', '[AT]A[ACGT]TCAT', 'CATCATCAT', 5],
                  ['ATTACTCATTATACTCTGAATGATGATGGCG', '[AT]A[ACGT]TCAT', 'CATCATCAT', 5],
                  ['GCTAGTCATCCGCGTAGCGATGATGATGAAA', '[AT]A[ACGT]TCAT', 'CATCATCAT', 5],
                  ['GATATTCATTTGAGCTTGCATGATGATGCTA', '[AT]A[ACGT]TCAT', 'CATCATCAT', 5],
                  ['ACAAGTCATTTTACTTTGTATGATGATGCCC', '[AT]A[ACGT]TCAT', 'CATCATCAT', 5]],
            columns=['sequence', 'fwd-primer', 'rev-primer', 'kmer-length'],
            index=pd.Index(['ref1', 'ref2', 'ref3', 'ref4', 'ref5'], 
                           name='db-seq')
            )
        known = pd.DataFrame(
            data=np.vstack([np.array([9, 13])] * 5),
            columns=['fwd_pos', 'rev_pos'],
            index=seqs.index
            )
        test = seqs.apply(_find_positions, axis=1)
        pdt.assert_frame_equal(known, test)

    def test_find_multiple_region_fragments(self):
        self.multiple_fragments.drop(
            columns=['fwd-pos', 'rev-pos', 'group-fwd', 'group-rev'], 
            inplace=True)
        self.multiple_fragments.set_index('db-seq', inplace=True)
        test = _find_multiple_region_fragments(self.multiple_fragments)
        pdt.assert_series_equal(test, self.multiple_fragments_trim)

    def test_find_single_region_fragments(self):
        known = pd.DataFrame.from_dict(orient='index', data={
            "seq01|seq02": {'fragment': 'GCGAAGCGGCTCAGG'}
            })
        known.index.set_names('clean_name', inplace=True)
        test = _find_single_region_fragments(self.single_fragment)
        pdt.assert_series_equal(test, known['fragment'])


if __name__ == '__main__':
    main()