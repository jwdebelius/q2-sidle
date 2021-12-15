import os
import warnings

import biom
import numpy as np
import pandas as pd

from qiime2 import Metadata, Artifact


def filter_with_aligned(alignment: pd.DataFrame,
                        table: biom.Table,
                        max_mismatch: int = None,
                        discarded: bool = True,
                        ) -> biom.Table:
    """
    Filters to the sequences retained during a sidle alignmeent

    Parameters
    ----------
    alignment: pd.DataFrame
        The alignment matrix for the dataset
    table : biom.Table
        The counts table for the region
    max_mismatch: int, optional
        The maximum mismatch to use, if not using the value from the alignment
    discard: bool, optional
        Whether the features that were discarded (discard) should be in the 
        table

    Returns
    -------
    biom.Table
        The filtered biom table
    """
    if max_mismatch is None:
        aligned_asvs = alignment['asv'].unique()
    else:
        aligned_asvs = \
            alignment.loc[alignment['mismatch'] <= max_mismatch, 'asv']
        aligned_asvs = aligned_asvs.unique()
    filt_table = table.copy().filter(aligned_asvs,
                                     axis='observation',
                                     invert=discarded,
                                     )
    return filt_table


def track_aligned_counts(
    region: str,
    regional_alignment: pd.DataFrame,
    regional_table: biom.Table,
    ) -> Metadata:
    """
    Tracks where sample counts have gone.
    """
    # Maps the regioional counts
    regional_counts = [
        _alignment_accounting(r, a, t)
        for (r, a, t) 
        in zip(*(region, regional_alignment, regional_table))
        ]
    total_table = regional_table[0].copy().merge(*regional_table[1:])
    total_counts = _alignment_accounting('total', 
                                         pd.concat(regional_alignment),
                                         total_table)

    counts = pd.concat(axis=1, objs=[total_counts, *regional_counts])
    counts.index.set_names('sample-id', inplace=True)

    return Metadata(counts)


def _alignment_accounting(region, alignment, table):
    keep_asvs = alignment['asv'].unique()
    table_n = table.copy().norm(axis='sample')
    regional_account = pd.DataFrame(
        data=np.vstack([
            table.sum(axis='sample'),
            table.filter(keep_asvs, axis='observation', 
                         inplace=False).sum(axis='sample'),
            table_n.filter(keep_asvs, axis='observation', 
                           inplace=False).sum(axis='sample')
        ]),
        columns=table.ids(axis='sample'),
        index=['starting counts', 'aligned counts', 'aligned percentage']
    ).T
    return regional_account.add_prefix(f'{region} ')
