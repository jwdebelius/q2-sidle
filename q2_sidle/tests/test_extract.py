from unittest import TestCase, main
import warnings

import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._extract import (prepare_extracted_region,
                               extract_regional_database,
                               _extract_region,
                               _find_primer_end,
                               _find_primer_start,
                               _tidy_sequences,
                               _collapse_duplicates,
                               _expand_degenerates,
                               _reverse_complement,
                               _trim_masked,
                               _build_id_map,
                               )


class ExtractTest(TestCase):
    def setUp(self):
        self.seq_array = pd.DataFrame.from_dict(orient='index', data={
            0: {'id_': '0', 0: 'A', 1: 'G', 2: 'T', 3: 'C'},
            1: {'id_': '1', 0: 'A', 1: 'R', 2: 'W', 3: 'S'},
            2: {'id_': '2', 0: 'C', 1: 'T', 2: 'W', 3: 'K'},
            3: {'id_': '3', 0: 'G', 1: 'T', 2: 'C', 3: 'M'},
            4: {'id_': '4', 0: 'A', 1: 'T', 2: 'G', 3: 'N'},
            }).set_index('id_')

        self.fwd_primer = 'WANTCAT'
        self.rev_primer = 'CATCATCAT'

        self.seqs = pd.DataFrame.from_dict(orient='index', data={
            'ref1': list('TGTACTCATGGACAATGTAATGATGATGGAT'),
            'ref2': list('ATTACTCATTATACTCTGAATGATGATGGCG'),
            'ref3': list('GCTAGTCATCCGCGTAGCGATGATGATGAAA'),
            'ref4': list('GATATTCATTTGAGCTTGCATGATGATGCTA'),
            'ref5': list('ACAAGTCATTTTACTTTGTATGATGATGCCC'),
            })
        self.seqs.index.set_names('kmer', inplace=True)

        self.full_seqs = Artifact.import_data(
            'FeatureData[Sequence]',
            pd.Series({'seq1': skbio.DNA(seq1, metadata={'id': 'seq1'}), 
                       'seq2': skbio.DNA(seq2, metadata={'id': 'seq2'}),
                       'seq3': skbio.DNA(seq3, metadata={'id': 'seq3'}), 
                       'seq4': skbio.DNA(seq4, metadata={'id': 'seq4'}),
                       'seq5': skbio.DNA(seq5, metadata={'id': 'seq5'}), 
                       'seq6': skbio.DNA(seq6, metadata={'id': 'seq6'})})
            )
        self.region_1 = Artifact.import_data(
            'FeatureData[Sequence]',
            pd.Series({
                'seq1 | seq2': skbio.DNA('GCGAAGCGGCTCAGG', 
                                         metadata={'id': 'seq1 | seq2'}),
                'seq3@0001': skbio.DNA('ATCCGCGTTGGAGTT', 
                                       metadata={'id': 'seq3@0001'}),
                'seq3@0002': skbio.DNA('TTCCGCGTTGGAGTT', 
                                       metadata={'id': 'seq3@0002'}),
                'seq5': skbio.DNA('CGTTTATGTATGCCC', metadata={'id': 'seq5'}),
                'seq6': skbio.DNA('CGTTTATGTATGCCT', metadata={'id': 'seq6'}),
                })
            )
        self.region1_map =  pd.DataFrame(
            data=np.array([['seq1', 'seq1 | seq2', 'WANTCAT-CATCATCAT'],
                           ['seq2', 'seq1 | seq2', 'WANTCAT-CATCATCAT'],
                           ['seq3@0001', 'seq3@0001', 'WANTCAT-CATCATCAT'],
                           ['seq3@0002', 'seq3@0002', 'WANTCAT-CATCATCAT'],
                           ['seq5', 'seq5', 'WANTCAT-CATCATCAT'],
                           ['seq6', 'seq6', 'WANTCAT-CATCATCAT'],
                           ]),
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq3', 'seq5', 'seq6'], 
                            name='db_seq'),
            columns=['seq_name', 'kmer', 'region'],
            )

    def test_prepared_extracted_region(self):
        trimmed = Artifact.import_data('FeatureData[Sequence]',
            pd.Series({
                'seq1': skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq1'}),
                'seq2': skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq2'}),
                'seq3': skbio.DNA('WTCCGCGTTGGAGTT', metadata={'id': 'seq3'}),
                'seq5': skbio.DNA('CGTTTATGTATGCCC', metadata={'id': 'seq5'}),
                'seq6': skbio.DNA('CGTTTATGTATGCCT', metadata={'id': 'seq6'}),
                })
            )
        test_seqs, test_map = \
            prepare_extracted_region(sequences=trimmed, 
                                     region='WANTCAT-CATCATCAT',
                                     trim_length=15,
                                     debug=True,
                                     )
        pdt.assert_series_equal(test_seqs.view(pd.Series).astype(str),
                                self.region_1.view(pd.Series).astype(str))
        pdt.assert_frame_equal(test_map,  self.region1_map)


    def test_extract_regional_database(self):
        test_seqs, test_map = \
            extract_regional_database(sequences=self.full_seqs,
                                      fwd_primer=self.fwd_primer,
                                      rev_primer=self.rev_primer,
                                      trim_length=15,
                                      primer_mismatch=1,
                                      debug=True,
                                      )
        pdt.assert_series_equal(test_seqs.view(pd.Series).astype(str),
                                self.region_1.view(pd.Series).astype(str))
        pdt.assert_frame_equal(test_map,  self.region1_map)

    def test_expand_degenerates_no_degen(self):
        in_mer = pd.DataFrame(
            data=np.array([list('AGTCCATGC'), 
                           list('TACGAGTGA'),
                           list('ACTCCATGC'),
                           list('AAAAAAAGT')]),
            )
        test = _expand_degenerates(in_mer.copy())
        pdt.assert_frame_equal(
            test.compute(), in_mer.rename({i: str(i) for i in in_mer.index}))

    def test_expand_degenerates(self):
        known = pd.DataFrame(
            data=np.array([['A', 'G', 'T', 'C'], # 0-0
                           ['A', 'A', 'A', 'C'], # 1-0
                           ['A', 'A', 'A', 'G'], # 1-1
                           ['A', 'A', 'T', 'C'], # 1-2
                           ['A', 'A', 'T', 'G'], # 1-3
                           ['A', 'G', 'A', 'C'], # 1-4
                           ['A', 'G', 'A', 'G'], # 1-5
                           ['A', 'G', 'T', 'C'], # 1-6
                           ['A', 'G', 'T', 'G'], # 1-7 
                           ['C', 'T', 'A', 'G'], # 2-0
                           ['C', 'T', 'A', 'T'], # 2-1
                           ['C', 'T', 'T', 'G'], # 2-2
                           ['C', 'T', 'T', 'T'], # 2-3
                           ['G', 'T', 'C', 'A'], # 3-0
                           ['G', 'T', 'C', 'C'], # 3-1
                           ['A', 'T', 'G', 'A'], # 4-0
                           ['A', 'T', 'G', 'C'], # 4-1
                           ['A', 'T', 'G', 'G'], # 4-2
                           ['A', 'T', 'G', 'T'], # 4-3
                           ]),
            index=pd.Index([u'0', 
                            u'1@01', u'1@02', u'1@03', u'1@04', u'1@05', 
                            u'1@06', u'1@07',u'1@08', 
                            u'2@01', u'2@02', u'2@03', u'2@04', 
                            u'3@01', u'3@02',
                            u'4@01', u'4@02', u'4@03', u'4@04',
                            ]),
        ).T.reset_index(drop=True).T

        test = _expand_degenerates(self.seq_array, pad=2)
        # Gets arund some delightful casting weirdness
        pdt.assert_frame_equal(known, test.compute(), check_index_type=False)

    def test_trim_masked(self):
        known = pd.DataFrame(
            data=[list('GT'),  list('RW'),  list('TW'),  
                  list('TG')],
            index=pd.Index(['0', '1', '2', '4'], name='id_')
            )
        start =pd.Series(
            data=[1, 1, 1, 1], 
            index=pd.Index(['0', '1', '2', '4'], name='id_')
            )
        end = pd.Series(
            data=[2, 2, 2, 2], 
            index=pd.Index(['0', '1', '2', '4'], name='id_')
            )
        test = _trim_masked(self.seq_array.loc[['0', '1', '2', '4']], 
                           start, 
                           end)
        pdt.assert_frame_equal(known, test)

    def test_reverse_complement(self):
        known = pd.DataFrame(
            data=np.array([['G', 'A', 'C', 'T'],
                           ['S', 'W', 'Y', 'T'],
                           ['M', 'W', 'A', 'G'],
                           ['K', 'G', 'A', 'C'],
                           ['N', 'C', 'A', 'T']]),
            index=pd.Index(['0', '1', '2', '3', '4'], name='id_'))
        test_ = _reverse_complement(self.seq_array)

        pdt.assert_frame_equal(known, test_)

    def test_find_primer_start_match(self):
        known = pd.Series({'pos': 0, 'mis': 0})
        test = _find_primer_start('Cats are awesome', '(Cat){e<=1}', adj=0)
        pdt.assert_series_equal(known, test)

    def test_find_primer_start_no_match(self):
        known = pd.Series({'pos': np.nan, 'mis': np.nan})
        test = _find_primer_start('Dogs are awesome', '(Cat){e<=1}', adj=0)
        pdt.assert_series_equal(known, test)

    def test_find_primer_end_match(self):
        known = pd.Series({'pos': 3, 'mis': 1})
        test = _find_primer_end('Rats are awesome', '(Cat){e<=1}', )
        pdt.assert_series_equal(known, test)

    def test_find_primer_end_no_match(self):
        known = pd.Series({'pos': np.nan, 'mis': np.nan})
        test = _find_primer_end('Iguanas are awesome', '(Cat){e<=1}')
        pdt.assert_series_equal(known, test)

    def test_extract_region_trim(self):
        k_trim = pd.DataFrame(
            data=np.array([list('GGACAATGTA'),
                           list('TATACTCTGA'),
                           list('CCGCGTAGCG'),
                           list('TTGAGCTTGC'),
                           list('TTTACTTTGT')]),
            index=pd.Index(self.seqs.index, name='kmer'),
            )
        t_trim = _extract_region(self.seqs, self.fwd_primer, self.rev_primer)
        pdt.assert_frame_equal(k_trim, t_trim.compute())

    def test_extract_region_no_trim(self):
        k_trim = pd.DataFrame(
            data=np.array([list('TACTCATGGACAATGTAATGATGAT'),
                           list('TACTCATTATACTCTGAATGATGAT'),
                           list('TAGTCATCCGCGTAGCGATGATGAT'),
                           list('TATTCATTTGAGCTTGCATGATGAT'),
                           list('AAGTCATTTTACTTTGTATGATGAT'),
                           ]),
            index=self.seqs.index,
        )
        k_trim.index.set_names('kmer', inplace=True)
        t_trim = _extract_region(self.seqs, self.fwd_primer, self.rev_primer, 
                                trim_fwd=False, trim_rev=False, 
                                trim_len='min')
        pdt.assert_frame_equal(k_trim, t_trim.compute())

    def test_extract_region_too_short(self):
        k_trim = pd.DataFrame(
            columns=self.seqs.columns,
            index=pd.Index([], name='kmer')
            )
        t_trim = _extract_region(
            self.seqs, self.fwd_primer, self.rev_primer, 
            trim_fwd=False, 
            trim_rev=False, 
            trim_len=40
            )
        pdt.assert_frame_equal(k_trim, t_trim.compute())

    def test_extract_region_set_length(self):
        k_trim = pd.DataFrame(
            data=np.array([list('GGACA'),
                           list('TATAC'),
                           list('CCGCG'),
                           list('TTGAG'),
                           list('TTTAC')]),
            index=self.seqs.index
            )
        k_trim.index.set_names('kmer', inplace=True)
        t_trim = _extract_region(self.seqs, self.fwd_primer, self.rev_primer, 
                                trim_len=5)
        pdt.assert_frame_equal(k_trim, t_trim.compute())

    def test_tidy_region(self):
        sequences =   pd.DataFrame.from_dict(orient='index', data={
                '0': list('GCGAAGCGGCTCAGG'),
                '1': list('GCGAAGCGGCTCAGG'),
                '2': list('GCGAAGCGGCTCAGG'),
                '3': list('WTCCGCGTTGGAGTT'),
                '4': list('ATCCGCGTTGGAGTT'),
                '5': list('ATCCGCGTTGGAGTS'),
                '6': list('CTCCGCGTTGGAGTT'),
                '7': list('MTCCGCGTTGGAGTT'),
                '8': list('DTCCGCGTTGGAGTT'),
                '9': list('CGTTTATGTATGCCC'),
                '10': list('CGTTTATGTATGCCT'),
                '11': list('CGTTTATGTATGCCT'),
                })
        known_dups = pd.DataFrame.from_dict(orient='index', data={
            '0 | 1 | 2' : list('GCGAAGCGGCTCAGG'),
            '10 | 11': list('CGTTTATGTATGCCT'),
            '3@0001 | 4 | 7@0001 | 8@0001': list('ATCCGCGTTGGAGTT'),
            '3@0002 | 8@0003':  list('TTCCGCGTTGGAGTT'),
            '5@0001': list('ATCCGCGTTGGAGTC'),
            '5@0002': list('ATCCGCGTTGGAGTG'),
            '6 | 7@0002': list('CTCCGCGTTGGAGTT'),
            '8@0002':  list('GTCCGCGTTGGAGTT'),
            '9': list('CGTTTATGTATGCCC'),
            })
        known_dups.index.set_names('seq_name', inplace=True)
        known_map = pd.DataFrame(
            data=np.array([['0', '0', '0 | 1 | 2', 'Bludhaven'],
                           ['1', '1', '0 | 1 | 2', 'Bludhaven'],
                           ['10', '10', '10 | 11', 'Bludhaven'],
                           ['11', '11', '10 | 11', "Bludhaven"],
                           ['2', '2', '0 | 1 | 2', 'Bludhaven'],
                           ['3', '3@0001', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['3', '3@0002', '3@0002 | 8@0003', 'Bludhaven'],
                           ['4', '4', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['5', '5@0001', '5@0001', 'Bludhaven'],
                           ['5', '5@0002', '5@0002', 'Bludhaven'],
                           ['6', '6', '6 | 7@0002', 'Bludhaven'],
                           ['7', '7@0001', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['7', '7@0002', '6 | 7@0002', 'Bludhaven'],
                           ['8', '8@0001', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['8', '8@0002', '8@0002', 'Bludhaven'],
                           ['8', '8@0003', '3@0002 | 8@0003', "Bludhaven"],
                           ['9', '9', '9', 'Bludhaven']], dtype=object),
            columns=['db_seq', 'seq_name', 'kmer', 'region'],
            )
        collapsed, map_ = _tidy_sequences([sequences], 'Bludhaven', 15)
        pdt.assert_frame_equal(collapsed, known_dups)
        pdt.assert_frame_equal(map_, known_map.set_index('db_seq'))

    def test_collapse_duplicates_sequences(self):
        expand_seqs = pd.DataFrame.from_dict(orient='index', data={
                '0' : list('GCGAAGCGGCTCAGG'),
                '1' : list('GCGAAGCGGCTCAGG'),
                '2' : list('GCGAAGCGGCTCAGG'),
                '3@0001': list('ATCCGCGTTGGAGTT'),
                '3@0002': list('TTCCGCGTTGGAGTT'),
                '4': list('ATCCGCGTTGGAGTT'),
                '5@0001': list('ATCCGCGTTGGAGTC'),
                '5@0002': list('ATCCGCGTTGGAGTG'),
                '6': list('CTCCGCGTTGGAGTT'),
                '7@0001': list('ATCCGCGTTGGAGTT'),
                '7@0002': list('CTCCGCGTTGGAGTT'),
                '8@0001': list('ATCCGCGTTGGAGTT'),
                '8@0002':  list('GTCCGCGTTGGAGTT'),
                '8@0003':  list('TTCCGCGTTGGAGTT'),
                '9': list('CGTTTATGTATGCCC'),
                '10': list('CGTTTATGTATGCCT'),
                '11': list('CGTTTATGTATGCCT'),
                })

        known_dups = pd.DataFrame.from_dict(orient='index', data={
            '0 | 1 | 2' : list('GCGAAGCGGCTCAGG'),
            '10 | 11': list('CGTTTATGTATGCCT'),
            '3@0001 | 4 | 7@0001 | 8@0001': list('ATCCGCGTTGGAGTT'),
            '3@0002 | 8@0003':  list('TTCCGCGTTGGAGTT'),
            '5@0001': list('ATCCGCGTTGGAGTC'),
            '5@0002': list('ATCCGCGTTGGAGTG'),
            '6 | 7@0002': list('CTCCGCGTTGGAGTT'),
            '8@0002':  list('GTCCGCGTTGGAGTT'),
            '9': list('CGTTTATGTATGCCC'),
            })
        known_dups.index.set_names('seq_name', inplace=True)

        test_dup = _collapse_duplicates(expand_seqs).compute()
        pdt.assert_frame_equal(known_dups, test_dup)

    def test_build_id_map(self):
        kmer = ['0 | 1 | 2', '10 | 11', '3@0001 | 4 | 7@0001 | 8@0001',
                '3@0002 | 8@0003', '5@0001', '5@0002', '6 | 7@0002',
                '8@0002', '9']
        known = pd.DataFrame(
            data=np.array([['0', '0', '0 | 1 | 2', 'Bludhaven'],
                           ['1', '1', '0 | 1 | 2', 'Bludhaven'],
                           ['10', '10', '10 | 11', 'Bludhaven'],
                           ['11', '11', '10 | 11', "Bludhaven"],
                           ['2', '2', '0 | 1 | 2', 'Bludhaven'],
                           ['3', '3@0001', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['3', '3@0002', '3@0002 | 8@0003', 'Bludhaven'],
                           ['4', '4', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['5', '5@0001', '5@0001', 'Bludhaven'],
                           ['5', '5@0002', '5@0002', 'Bludhaven'],
                           ['6', '6', '6 | 7@0002', 'Bludhaven'],
                           ['7', '7@0001', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['7', '7@0002', '6 | 7@0002', 'Bludhaven'],
                           ['8', '8@0001', '3@0001 | 4 | 7@0001 | 8@0001', 
                            'Bludhaven'],
                           ['8', '8@0002', '8@0002', 'Bludhaven'],
                           ['8', '8@0003', '3@0002 | 8@0003', "Bludhaven"],
                           ['9', '9', '9', 'Bludhaven']], dtype=object),
            columns=['db_seq', 'seq_name', 'kmer', 'region'],
            )

        test = _build_id_map(kmer, 'Bludhaven')
        pdt.assert_frame_equal(test, known.set_index('db_seq'))

seq1 = 'ACTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGACACCTCGTAAGAGAGGCTGAATCCTGACGTGAAC'
# #seq2: same first region as seq1, 3 nt difference in second
seq2 = 'ACTAGTCATGCGAAGCGGCTCAGGATGATGATGAAGACACCTCGTAAGAGTTTCTGAATCCTGACGTGAAC'
# seq3: degeneracy in first region
seq3 = 'CATAGTCATWTCCGCGTTGGAGTTATGATGATGAAACCACCTCGTCCCAGTTCCGCGCTTCTGACGTGCAC'
# Skipped in 1st region
seq4 = 'GATTTTTTTATGTTCGCATGCGGAATGATGATGCCGACACCTCGTCCCGGAGACGAGAGGCTGACGTGAGC'
seq5 = 'CATAGTCATCGTTTATGTATGCCCATGATGATGCGAGCACCTCGTATGGATGTAGAGCCACTGACGTGCGG'
# # seq6: 1 nt diff in region 1, 3 nt in region 2
seq6 = 'CATAGTCATCGTTTATGTATGCCTATGATGATGCGAGCACCTCGTAAAAATGTAGAGCCACTGACGTGCGG'

if __name__ == '__main__':
    main()

