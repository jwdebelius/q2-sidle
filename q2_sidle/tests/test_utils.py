from unittest import TestCase, main

import os

import dask
import numpy as np
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact, Metadata
from qiime2.plugin import ValidationError

from q2_sidle._utils import (_convert_seq_block_to_dna_fasta_format,
                             _convert_generator_to_delayed_seq_block,
                             _convert_generator_to_seq_block,
                             _check_manifest,
                             _read_manifest_files,
                             _count_degenerates,
                             _find_primer_end,
                             _find_primer_start,
                             )
import q2_sidle.tests.test_set as ts
from q2_types.feature_data import DNAIterator, DNAFASTAFormat



class UtilTest(TestCase):
    def setUp(self):
        self.seq_block = [pd.DataFrame(
            data=[list('CATS'), list('WANT'), list("CANS")],
            index=['0', '1', '2']
        )]
        self.skbio_series = pd.Series(data={
            "0": skbio.DNA('CATS', metadata={'id': '0'}),
            "1": skbio.DNA('WANT', metadata={'id': '1'}),
            "2": skbio.DNA('CANS', metadata={'id': '2'}),
             })
        self.seq_artifact = Artifact.import_data('FeatureData[Sequence]', 
                                                 self.skbio_series, 
                                                 pd.Series)

    def test_convert_seq_block_to_dna_fasta_format(self):
        ff = _convert_seq_block_to_dna_fasta_format(self.seq_block)
        test = ff.view(pd.Series)
        pdt.assert_series_equal(self.skbio_series.astype(str), 
                                self.skbio_series.astype(str))

    def test_convert_generator_to_delayed_seq_block(self):
        gen = self.seq_artifact.view(DNAIterator)
        test = _convert_generator_to_delayed_seq_block(gen, chunksize=5)
        self.assertTrue(test, list)
        self.assertEqual(len(test), 1)
        pdt.assert_frame_equal(self.seq_block[0], dask.compute(test[0])[0])

    def test_convert_generator_to_seq_block(self):
        gen = self.seq_artifact.view(DNAIterator)
        test = _convert_generator_to_seq_block(gen, chunksize=5)
        self.assertTrue(test, list)
        self.assertEqual(len(test), 1)
        pdt.assert_frame_equal(self.seq_block[0], test[0])

    def test_count_degenerates(self):
        seq_array = pd.DataFrame.from_dict(orient='index', data={
            0: {'id_': '0', 0: 'A', 1: 'G', 2: 'T', 3: 'C'},
            1: {'id_': '1', 0: 'A', 1: 'R', 2: 'W', 3: 'S'},
            3: {'id_': '3', 0: 'G', 1: 'T', 2: 'C', 3: 'M'},
            4: {'id_': '4', 0: 'A', 1: 'T', 2: 'G', 3: 'N'},
            }).set_index('id_')
        known = pd.Series(data=[0, 3, 1, 1], 
                          index=pd.Index(['0', '1', '3', '4'], name='id_'), 
                          dtype=int)
        test = _count_degenerates(seq_array)
        pdt.assert_series_equal(known, test)

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
            data=np.array([['Bruce', 'Wayne', None, '0'],
                          ['Dick', 'Grayson', '29', '0'],
                          ['Jason', 'Todd', '24', '0'],
                          ['Tim', 'Drake', '20', '0'],
                          ['Barbara', 'Gordon', '32', '1'],
                          ['Stephanie', 'Brown', '19', '1'],
                          ['Cassandra', 'Cain-Wayne', '22', '1'],
                          ['Damian', 'Wayne', '12', '1'],
                          ]),
            index=pd.Index(['Batman', 'Nightwing', 'Red Hood', 'Red Robin',
                            'Oracle', 'Batgirl', 'Black Bat', 'Robin'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table', 'region-order']
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
            data=np.array([['Bruce', 'Wayne', '45', '0'],
                          ['Dick', 'Grayson', '29', '0'],
                          ['Jason', 'Todd', '24', '0'],
                          ['Tim', 'Drake', '20', '0'],
                          ['Barbara', 'Gordon', '32', '1'],
                          ['Stephanie', 'Brown', '19', '1'],
                          ['Cassandra', 'Cain-Wayne', '22', '1'],
                          ['Damian', 'Wayne', '12', '1'],
                          ]),
            index=pd.Index(['Batman', 'Nightwing', 'Red Hood', 'Red Robin',
                            'Oracle', 'Batgirl', 'Black Bat', 'Robin'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table',  
                     'region-order'],
        ))
        with self.assertRaises(ValidationError) as err:
            _check_manifest(manifest)
        self.assertEqual(str(err.exception), 
                         ('All paths in the manifest must be unique.'
                          ' Please check your filepaths')
                         )

    def test_check_manifest_exists(self):
        manifest = Metadata(pd.DataFrame(
            data=np.array([['Bruce', 'Wayne', '45', '0'],
                           ['Dick', 'Grayson', '29', '0'],
                           ['Jason', 'Todd', '24', '0'],
                           ]),
            index=pd.Index(['Batman', 'Nightwing', 'Red Hood'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table', 
                     'region-order'],
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
                           os.path.join(base_dir, 'region1_align.qza'),
                           'Bludhaven',
                           ]]),
            index=pd.Index(['test'], 
                            name='id'),
            columns=['kmer-map', 'alignment-map', 'frequency-table', 
                     'region-order'],
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


if __name__ == '__main__':
    main()