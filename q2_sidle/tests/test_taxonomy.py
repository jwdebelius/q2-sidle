from unittest import TestCase, main

import warnings

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.util.testing as pdt
import skbio

from qiime2 import Artifact

from q2_sidle._taxonomy import (reconstruct_taxonomy,
                                # _combine_taxonomic_annotation,   
                                # _parse_generic,
                                # _parse_greengenes,
                                )
import q2_sidle.tests.test_set as ts

class TaxonomyTest(TestCase):
    def setUp(self):
        self.gg_taxonomy = pd.Series(
            data={'seq01': ('k__DCU; p__Superhero; c__Gotham; o__Civillian; ''f__; g__; s__'),
                  'seq02': ('k__DCU; p__Superhero; c__Gotham; o__Civillian; ''f__Doctor; g__Thompson; s__Leslie'),                   
                  'seq03': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Batgirl; g__Gordon; s__Barbara'),
                  'seq04': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Batgirl; g__Gordon; s__Barbara'),
                  'seq05': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__; g__; s__'),
                  'seq06': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Oracle; g__Gordon; s__'),
                  'seq07': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Oracle; g__; s__'),                   
                  'seq08': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Batman; g__; s__'),
                  'seq09': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Batman; g__Wayne; s__Bruce'),
                  'seq10': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Batman; g__Grayson; s__'),
                  'seq11': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Batman; g__; s__'),
                  'seq12': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Robin; g__Grayson; s__Dick'),
                  'seq13': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Robin; g__Todd; s__'),
                  'seq14': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Robin; g__Drake; s__'),
                  'seq15': ('k__DCU; p__Superhero; c__Gotham; o__Batfamily; ''f__Robin; g__Brown; s__')},
            name='Taxon'
        )
        self.gg_taxonomy.index.set_names('Feature ID', inplace=True)
        
        self.silva_taxonomy = pd.Series(
            data={'seq01': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Civillian;D_4__ambigious taxa;D_5__ambigious taxa;D_6__ambigious taxa',
                  'seq02': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Civillian;D_4__Doctor;D_5__Thompson;D_6__Leslie',
                  'seq03': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batgirl;D_5__Gordon;D_6__Barbara',
                  'seq04': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batgirl;D_5__Gordon;D_6__Barbara',
                  'seq05': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__uncultured_child;D_5__uncultured_child;D_6__uncultured_child',
                  'seq06': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Oracle;D_5__Gordon;D_6__bat_metagenome',
                  'seq07': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Oracle;D_5__bat_metagenome;D_6__bat_metagenome',
                  'seq08': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batman;D_5__ambigious batman;D_6__ambigious batman',
                  'seq09': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batman;D_5__Wayne;D_6__Bruce',
                  'seq10': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batman;D_5__Grayson;D_6__ambigious grayson',
                  'seq11': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Oracle;D_5__bat_metagenome;D_6__ambigious_child',
                  'seq12': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Grayson;D_6__Dick',
                  'seq13': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Todd;D_6__uncultured_child',
                  'seq14': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Drake;D_6__uncultured_child',
                  'seq15': 'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Brown;D_6__ambigious_child'},
            name='Taxon'
        )
        self.silva_taxonomy.index.set_names('Feature ID', inplace=True)
        self.seq_map = pd.Series({"seq01": "seq01", 
                                  'seq03': 'seq03|seq04', 
                                  'seq04': 'seq03|seq04', 
                                  'seq06': 'seq06|seq07',
                                  'seq07': 'seq06|seq07',
                                  'seq08': 'seq08|seq09|seq10',
                                  'seq09': 'seq08|seq09|seq10',
                                  'seq10': 'seq08|seq09|seq10',
                                  'seq12': 'seq12',
                                  'seq13': 'seq13|seq14',
                                  'seq14': 'seq13|seq14'},
                                name='db_seq')
        self.seq_map.index.set_names('clean_name', inplace=True)

    def test_reconstruct_taxonomy_none_ignore_ignore(self):
        known = pd.Series(
            data=['k__DCU;p__Superhero;c__Gotham;o__Civillian;f__;g__;s__',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Batgirl;g__Gordon;s__Barbara',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Oracle;g__Gordon|g__;g__Gordon|g__',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Batman;g__|g__Wayne|g__Grayson;g__|g__Wayne|g__Grayson',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Robin;g__Grayson;s__Dick',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Robin;g__Todd|g__Drake;g__Todd|g__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
            )
        test = reconstruct_taxonomy(
            taxonomy=self.gg_taxonomy.apply(lambda x: x.replace(" ", "")), 
            reconstruction_map=self.seq_map,
            database='none',
            define_missing='ignore',
            ambiguity_handling='ignore',
            )
        pdt.assert_series_equal(known, test)

    def test_reconstruct_taxonomy_none_warning(self):
        known = pd.Series(
            data=['k__DCU;p__Superhero;c__Gotham;o__Civillian;f__;g__;s__',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Batgirl;g__Gordon;s__Barbara',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Oracle;g__Gordon|g__;g__Gordon|g__',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Batman;g__|g__Wayne|g__Grayson;g__|g__Wayne|g__Grayson',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Robin;g__Grayson;s__Dick',
                  'k__DCU;p__Superhero;c__Gotham;o__Batfamily;f__Robin;g__Todd|g__Drake;g__Todd|g__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
            )
        with warnings.catch_warnings(record=True) as w:
            test = reconstruct_taxonomy(
                taxonomy=self.gg_taxonomy.apply(lambda x: x.replace(" ", "")), 
                reconstruction_map=self.seq_map,
                database='none',
                define_missing='merge',
                ambiguity_handling='missing',
            )
            self.assertEqual(len(w), 2)  
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertEqual(str(w[0].message), 
                              'When no database is specified, missing values '
                              'are ignored by default')
            self.assertTrue(issubclass(w[1].category, UserWarning))
            self.assertEqual(str(w[1].message), 
                              'When no database is specified, '
                              'ambiguious values are ignored by default')
        pdt.assert_series_equal(known, test)

    def test_reconstruct_taxonomy_delim_error(self):
        with self.assertRaises(ValueError) as err:
            test = reconstruct_taxonomy(
                taxonomy=self.gg_taxonomy.apply(lambda x: x.replace(' ', '')), 
                reconstruction_map=self.seq_map,
                database='greengenes',
                )
        self.assertEqual(
            str(err.exception), 
            'Only one taxonomic level was found. Please check your database '
            'and delimiter.'
        )

    def test_reconstruct_taxonomy_gg_ignore_ignore(self):
        known = pd.Series(
            data=['k__DCU; p__Superhero; c__Gotham; o__Civillian; f__; g__; s__',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batgirl; g__Gordon; s__Barbara',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Oracle; g__Gordon|g__; g__Gordon|g__',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batman; g__|g__Wayne|g__Grayson; g__|g__Wayne|g__Grayson',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Grayson; s__Dick',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Todd|g__Drake; g__Todd|g__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
            )
        test = reconstruct_taxonomy(taxonomy=self.gg_taxonomy, 
                                    reconstruction_map=self.seq_map,
                                    database='greengenes',
                                    define_missing='ignore',
                                    ambiguity_handling='ignore',
                                    )
        pdt.assert_series_equal(known, test)

    def test_reconstruct_taxonomy_gg_inheriet_missing(self):
        known = pd.Series(
            data=['k__DCU; p__Superhero; c__Gotham; o__Civillian; o__Civillian; o__Civillian; o__Civillian',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batgirl; g__Gordon; s__Barbara',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Oracle; g__Gordon|f__Oracle; g__Gordon|f__Oracle',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batman; f__Batman|g__Wayne|g__Grayson; f__Batman|g__Wayne|g__Grayson',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Grayson; s__Dick',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Todd|g__Drake; g__Todd|g__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
            )
        with warnings.catch_warnings(record=True) as w:
            test = reconstruct_taxonomy(
                taxonomy=self.gg_taxonomy, 
                reconstruction_map=self.seq_map,
                database='greengenes',
                define_missing='inherit',
                ambiguity_handling='missing',
            )
            self.assertEqual(len(w), 1)  
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertEqual(str(w[0].message), 
                             'Greengenes does not include ambigious taxa. '
                             'The ambiguity handling will be ignored.')
        pdt.assert_series_equal(known, test)

    def test_reconstruct_taxonomy_gg_merge_ignore(self):
        known = pd.Series(
            data=['k__DCU; p__Superhero; c__Gotham; o__Civillian; o__Civillian; o__Civillian; o__Civillian',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batgirl; g__Gordon; s__Barbara',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Oracle; g__Gordon; g__Gordon',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Batman; g__Wayne|g__Grayson; g__Wayne|g__Grayson',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Grayson; s__Dick',
                  'k__DCU; p__Superhero; c__Gotham; o__Batfamily; f__Robin; g__Todd|g__Drake; g__Todd|g__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
            )
        test = reconstruct_taxonomy(
            taxonomy=self.gg_taxonomy, 
            reconstruction_map=self.seq_map,
            database='greengenes',
        )
        pdt.assert_series_equal(known, test)

    def test_test_reconstruct_taxonomy_silva_ignore_ignore(self):
        known = pd.Series(
            data=['D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Civillian;D_4__ambigious taxa;D_5__ambigious taxa;D_6__ambigious taxa',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batgirl;D_5__Gordon;D_6__Barbara',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Oracle;D_5__Gordon|D_5__bat_metagenome;D_5__Gordon|D_5__bat_metagenome',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batman;D_5__ambigious batman|D_5__Wayne|D_5__Grayson;D_5__ambigious batman|D_5__Wayne|D_5__Grayson',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Grayson;D_6__Dick',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Todd|D_5__Drake;D_5__Todd|D_5__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
        )
        test = reconstruct_taxonomy(
            taxonomy=self.silva_taxonomy, 
            reconstruction_map=self.seq_map,
            database='silva',
            define_missing='ignore',
            ambiguity_handling='ignore',
        )
        pdt.assert_series_equal(known, test)

    def test_reconstruct_taxonomy_silva_inherit_missing(self):
        known = pd.Series(
            data=['D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Civillian;D_3__Civillian;D_3__Civillian;D_3__Civillian',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batgirl;D_5__Gordon;D_6__Barbara',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Oracle;D_5__Gordon|D_4__Oracle;D_5__Gordon|D_4__Oracle',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batman;D_4__Batman|D_5__Wayne|D_5__Grayson;D_4__Batman|D_5__Wayne|D_5__Grayson',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Grayson;D_6__Dick',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Todd|D_5__Drake;D_5__Todd|D_5__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
        )
        test = reconstruct_taxonomy(
            taxonomy=self.silva_taxonomy, 
            reconstruction_map=self.seq_map,
            database='silva',
            define_missing='inherit',
            ambiguity_handling='missing',
        )
        pdt.assert_series_equal(known, test)

    def test_reconstruct_taxonomy_silva_merge_missing(self):
        known = pd.Series(
            data=['D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Civillian;D_4__ambigious taxa;D_5__ambigious taxa;D_6__ambigious taxa',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batgirl;D_5__Gordon;D_6__Barbara',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Oracle;D_5__Gordon;D_5__Gordon',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Batman;D_5__ambigious batman|D_5__Wayne|D_5__Grayson;D_5__ambigious batman|D_5__Wayne|D_5__Grayson',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Grayson;D_6__Dick',
                  'D_0__DCU;D_1__Superhero;D_2__Gotham;D_3__Batfamily;D_4__Robin;D_5__Todd|D_5__Drake;D_5__Todd|D_5__Drake',
                  ],
            index=pd.Index(['seq01', 'seq03|seq04', 'seq06|seq07', 'seq08|seq09|seq10', 'seq12', 'seq13|seq14'], name='Feature ID'),
            name='Taxon',
        )
        test = reconstruct_taxonomy(
            taxonomy=self.silva_taxonomy, 
            reconstruction_map=self.seq_map,
            database='silva',
        )
        pdt.assert_series_equal(known, test)

if __name__ == '__main__':
    main()

