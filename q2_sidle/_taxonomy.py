import os
import warnings

import numpy as np
import pandas as pd

from qiime2 import Metadata, Artifact
from qiime2.plugin import ValidationError
from q2_sidle._utils import (database_params)



def reconstruct_taxonomy(reconstruction_map: pd.Series, 
                         taxonomy: pd.Series, 
                         database: str='none', 
                         define_missing: str='merge',
                         ambiguity_handling: str='ignore',
                         ) -> pd.Series:
    """
    Reconstructs the taxonomic annotation based on a sidle database by 
    identifying the lowest taxonomic level  where the taxonomic annotation 
    diverges for two sequences

    Parameters
    ----------
    reconstruction_map : pd.Series
        The relationship between raw  sequences and the reconstructed sidle
        database
    taxonomy: pd.Series
        A taxonomic description of each sequence
    database: {'greengenes', 'silva', 'none'}
        The database used taxonomy. This is important for selecting the 
        correct taxonomic delimiter and for removing missing sequences.
    define_missing: {'merge'; 'inherit'; 'ignore'}
        Taxonomic strings may be missing information (for example  `g__` in 
        greengenes  or `D_5__uncultured bacteria` in Silva). These can be 
        ignored  (`"ignore"`) and treated like any other taxonomic 
        designation; they can be first inherited in merged sequences 
        (`"merge`"), where, when there are two strings being merged and 
        one has a missing level, the missing level is taken form the 
        defined one, or they can be inherited from the previous level 
        (`"inherit"`) first, and then merged.
    ambiguity_handling: {'missing', 'ignore'}
        whether "ambigious taxa" (Silva-specific) should be treated as 
        missing values (`"missing`") or ignored (`"ignore"`)

    Returns
    -------
    pd.Series
        A series describing the new taxonomy 
    """
    if (database == 'none') & (define_missing != 'ignore'):
        warnings.warn('When no database is specified, '
                      'missing values are ignored by default', UserWarning)
    if (database == 'none') & (ambiguity_handling != 'ignore'):
        warnings.warn('When no database is specified, '
                      'ambiguious values are ignored by default', UserWarning)
    if (database == 'greengenes') and (ambiguity_handling != 'ignore'):
        warnings.warn('Greengenes does not include ambigious taxa. The '
                       'ambiguity handling will be ignored.', UserWarning)

    # Filters the taxonomy and converts to levels
    db_lookup = database_params.get(database, 'none')
    delim = db_lookup['delim']
    def split_taxonomy(x):
        return pd.Series([s.strip(' ') for s in  x.split(delim)])
    taxonomy = taxonomy.loc[reconstruction_map.index]
    taxonomy = taxonomy.apply(split_taxonomy)
    taxonomy.index.set_names('Feature ID', inplace=True)

    if len(taxonomy.columns) == 1:
        raise ValueError('Only one taxonomic level was found. Please check '
                         'your database and delimiter.')

    # Finds the undefined levels
    defined_f = db_lookup['defined']
    undefined_levels = ~pd.concat(axis=1, objs=[
        taxonomy[c].apply(defined_f) for c in taxonomy.columns
    ])
    ambigious_levels = pd.concat(axis=1, objs=[
        taxonomy[c].apply(lambda x: 'ambig' in x) for c in taxonomy.columns
        ])
    ambigious_levels = ambigious_levels.cummax(axis=1)
    undefined = (undefined_levels | 
                 (ambigious_levels & (ambiguity_handling == 'missing'))
                 ).astype(bool)
    

    # Filters missing taxa and  hanldes initial inherietence.
    if define_missing != 'ignore':
        taxonomy.mask(undefined, np.nan, inplace=True)
    
    if define_missing == 'inherit':
        taxonomy.fillna(method='ffill', axis=1, inplace=True)

    # Combines the taxonomy across multiple  levels
    def _combine_f(x):
        if pd.isnull(x).all():
            return np.nan
        else:
            return '|'.join(np.sort(x.dropna().unique()))

    def _combine_taxa(g):
        """Help function ot tidy taxonomy"""
        if len(g) == 1:
            return g.iloc[0]
        else:
            return g.apply(_combine_f)

    taxonomy['clean_name'] = reconstruction_map
    collapsed  = taxonomy.groupby('clean_name').apply(_combine_taxa)
    collapsed.drop(columns=['clean_name'], inplace=True)

    # Finds splits in the data
    disjoint = pd.concat(axis=1, objs=[
        collapsed[c].apply(lambda x: True if pd.isnull(x) else '|' in x) 
        for c in collapsed.columns
        ]).cummax(axis=1)
    # Set up inherietence so you inheriet the first split in each row 
    # of the data
    disjoint_inheriet = (disjoint.cummax(axis=1) & 
                         ~((disjoint.cumsum(axis=1) == 1) & disjoint))
    collapsed.mask(disjoint_inheriet, np.nan, inplace=True)

    # Does nan inherietence
    collapsed.fillna(method='ffill', axis=1,  inplace=True)
    # Returns  the summarized taxonomy
    new_taxa = collapsed.apply(lambda x: delim.join(list(x.values)), axis=1)

    new_taxa.name = 'Taxon'
    new_taxa.index.set_names('Feature ID', inplace=True)

    return new_taxa


