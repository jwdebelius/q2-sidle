from unittest import TestCase, main
import warnings

import dask.dataframe as dd
import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact
from q2_types.feature_data import DNAIterator, DNAFASTAFormat

from q2_sidle._extract import (prepare_extracted_region,
                               _artifical_trim,
                               _block_seqs,
                               _collapse_all_sequences,
                               _condense_seqs,
                               _expand_degenerate_gen,
                               _expand_ids,
                               _split_ids,
                               )
from q2_sidle.tests import test_set as ts


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

        self.trimmed = Artifact.import_data('FeatureData[Sequence]',
            pd.Series({
                'seq1': skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq1'}),
                'seq2': skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq2'}),
                'seq3': skbio.DNA('WTCCGCGTTGGAGTT', metadata={'id': 'seq3'}),
                'seq5': skbio.DNA('CGTTTATGTATGCCC', metadata={'id': 'seq5'}),
                'seq6': skbio.DNA('CGTTTATGTATGCCT', metadata={'id': 'seq6'}),
                })
            )
        self.region_1 = ts.region1_db_seqs
        self.region1_map = ts.region1_db_map

        self.seq_block = pd.DataFrame(
            data=[['seq1', 'GCGAAGCGGCTCAGG', 'seq1'],
                  ['seq2', 'GCGAAGCGGCTCAGG', 'seq2'],
                  ['seq3@0001', 'ATCCGCGTTGGAGTT', 'seq3'],
                  ['seq3@0002', 'TTCCGCGTTGGAGTT', 'seq3'],
                  ['seq5', 'CGTTTATGTATGCCC', 'seq5'],
                  ['seq6', 'CGTTTATGTATGCCT', 'seq6'],
                  ],
            columns=['seq-name', 'sequence', 'db-seq']
            )
        self.amplicon = pd.DataFrame(
            data=[['seq1', 'GCGAAGCGGCTCAGG'],
                  ['seq2', 'GCGAAGCGGCTCAGG'],
                  ['seq3@0001', 'ATCCGCGTTGGAGTT'],
                  ['seq3@0002', 'TTCCGCGTTGGAGTT'],
                  ['seq5', 'CGTTTATGTATGCCC'],
                  ['seq6', 'CGTTTATGTATGCCT'],
                  ],
            columns=['seq-name', 'amplicon']
            )
        self.amplicon_r = pd.DataFrame(
            data=[['seq1', 'TCAGG'], ['seq2', 'TCAGG'], ['seq3@0001', 'GAGTT'],
                  ['seq5', 'TGCCC'], ['seq6', 'TGCCT']],
            columns=['seq-name', 'amplicon'],
            index=np.array([0, 1, 2, 4, 5]),
            )
        self.reverse_seqs = pd.Series({
            'seq1|seq2': 'CCTGA',
            'seq3@0001': 'AACTC',
            'seq5': 'GGGCA',
            'seq6': 'AGGCA',
            })
        self.group_forward = \
            self.region_1.view(pd.Series).astype(str).reset_index()
        self.group_forward.columns = ['seq-name', 'seq']
        self.group_forward['seq-name'] = \
            self.group_forward['seq-name'].apply(lambda x: '>%s' % x)
        self.group_forward['amplicon'] = self.group_forward['seq']

    def test_prepared_extracted_region(self):
        test_seqs, test_map = \
            prepare_extracted_region(sequences=self.trimmed, 
                                     region='Bludhaven',
                                     trim_length=15,
                                     debug=True,
                                     fwd_primer='WANTCAT',
                                     rev_primer='CATCATCAT',
                                     reverse_complement_rev=True,
                                     reverse_complement_result=False,
                                     )

        pdt.assert_series_equal(
            test_seqs.view(pd.Series).astype(str).sort_values(),
            self.region_1.view(pd.Series).astype(str).sort_values()
        )
        test_map.sort_values(['db-seq', 'seq-name'], inplace=True)
        test_map = test_map[['seq-name', 'kmer', 'region', 
                             'fwd-primer', 'rev-primer', 'kmer-length']]
        pdt.assert_frame_equal(test_map,  
                              self.region1_map.view(pd.DataFrame))

    def test_prepared_extracted_region_rc(self):
        test_seqs, test_map = \
            prepare_extracted_region(sequences=self.trimmed, 
                                     region='Bludhaven-rev',
                                     trim_length=-5,
                                     debug=True,
                                     fwd_primer='ATGATGATG',
                                     rev_primer='ATCANTW',
                                     reverse_complement_rev=False,
                                     reverse_complement_result=True,
                                     )
        known_map = pd.DataFrame(
            data=[['seq1', 'seq1|seq2', 'Bludhaven-rev', 'ATGATGATG', 'ATCANTW', -5],
                  ['seq2', 'seq1|seq2', 'Bludhaven-rev', 'ATGATGATG', 'ATCANTW', -5],
                  ['seq3@0001', 'seq3@0001', 'Bludhaven-rev', 'ATGATGATG', 'ATCANTW', -5],
                  ['seq5', 'seq5', 'Bludhaven-rev', 'ATGATGATG', 'ATCANTW', -5],
                  ['seq6', 'seq6', 'Bludhaven-rev', 'ATGATGATG', 'ATCANTW', -5]],
            index=pd.Index(['seq1', 'seq2', 'seq3', 'seq5', 'seq6'], name='db-seq'),
            columns=['seq-name', 'kmer', 'region', 'fwd-primer', 'rev-primer', 'kmer-length']
        )
        pdt.assert_series_equal(
            test_seqs.view(pd.Series).astype(str).sort_index(), 
            self.reverse_seqs)
        test_map.sort_values(['db-seq', 'seq-name'], inplace=True)
        test_map = test_map[['seq-name', 'kmer', 'region', 
                             'fwd-primer', 'rev-primer', 'kmer-length']]
        pdt.assert_frame_equal(test_map, known_map)

    def test_artifical_trim_fwd(self):
        test = _artifical_trim(self.seq_block, 15).compute()
        pdt.assert_frame_equal(test, self.amplicon)

    def test_artifical_trim_rev(self):
        test = _artifical_trim(self.seq_block, -5).compute()
        pdt.assert_frame_equal(test, self.amplicon_r)

    def test_block_seqs(self):
        test = _block_seqs(self.trimmed.view(pd.Series).values).compute()
        pdt.assert_frame_equal(self.seq_block, test)

    def test_collapse_all_sequences_fwd(self):
        condensed = dd.from_pandas(self.amplicon, chunksize=5000)
        test_ff, test_group2 = _collapse_all_sequences(condensed, False)
        self.assertTrue(isinstance(test_ff, DNAFASTAFormat))
        pdt.assert_series_equal(
            test_ff.view(pd.Series).astype(str), 
            self.region_1.view(pd.Series).astype(str)
            )
        pdt.assert_frame_equal(
            test_group2.reset_index(drop=True), 
            self.group_forward[['amplicon', 'seq-name', 'seq']])

    def test_collapse_all_seqs_rev(self):
        condensed = dd.from_pandas(chunksize=5000, data=self.amplicon_r)
        known_grouped = self.reverse_seqs.reset_index()
        known_grouped.columns = ['seq-name', 'seq']
        known_grouped['seq-name'] = \
            known_grouped['seq-name'].apply(lambda x: '>%s' % x)
        known_grouped['amplicon'] = ['TCAGG', 'GAGTT', 'TGCCC', 'TGCCT']

        test_ff, test_group2 = _collapse_all_sequences(condensed, True)
        pdt.assert_series_equal(self.reverse_seqs, 
                                test_ff.view(pd.Series).astype(str))
        pdt.assert_frame_equal(known_grouped[['amplicon', 'seq-name', 'seq']],
                               test_group2.reset_index(drop=True))

    def test_condense_seqs(self):
        test = _condense_seqs(self.amplicon).compute()
        known = self.region_1.view(pd.Series).astype(str).reset_index()
        known.columns = ['seq-name', 'amplicon']
        known.sort_values('amplicon', inplace=True)
        known.reset_index(inplace=True, drop=True)
        pdt.assert_frame_equal(
            test.reset_index(drop=True), 
            known[['amplicon', 'seq-name']].sort_values('amplicon')
            )

    def test_expand_degenerate_gen_no_degen(self):
        seq = skbio.DNA('GCGAAGCGGCTCAGG', metadata={'id': 'seq1'})
        known = pd.Series({'seq1': 'GCGAAGCGGCTCAGG'})
        test = _expand_degenerate_gen(seq)
        pdt.assert_series_equal(test, known)

    def test_expand_degenerate_gen_degen(self):
        seq = skbio.DNA('WTCCGCGTTGGAGTT', metadata={'id': 'seq3'})
        known = pd.Series({'seq3@0001': 'ATCCGCGTTGGAGTT',
                           'seq3@0002': 'TTCCGCGTTGGAGTT'})
        test = _expand_degenerate_gen(seq)
        pdt.assert_series_equal(test, known)

    def test_expand_ids(self):
        test = _expand_ids(self.group_forward, self.fwd_primer, 'ATGATGATG',
                           'Bludhaven', 15, 1000).compute()
        test = test[['db-seq', 'seq-name', 'kmer', 'region', 
                     'fwd-primer', 'rev-primer', 'kmer-length']]
        test.set_index('db-seq', inplace=True)
        pdt.assert_frame_equal(test.sort_index(), 
                               self.region1_map.view(pd.DataFrame))

    def test_split_ids(self):
        ids = pd.Series(['seq1|seq2', 'seq3@0001', 'seq3@0002', 
                         'seq5', 'seq6'], name='kmer')
        known = pd.DataFrame(
            data=[['seq1|seq2', 'seq1'],
                  ['seq1|seq2', 'seq2'],
                  ['seq3@0001', 'seq3@0001'],
                  ['seq3@0002', 'seq3@0002'],
                  ['seq5', 'seq5'],
                  ['seq6', 'seq6'],
                  ],
            columns=['kmer', 'seq-name']
            )
        test = _split_ids(ids).compute()
        test.sort_values(['kmer', 'seq-name'], inplace=True)
        test.reset_index(drop=True, inplace=True)
        pdt.assert_frame_equal(known, test)



if __name__ == '__main__':
    main()

 