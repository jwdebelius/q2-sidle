import os
import warnings

import biom
import numpy as np
import pandas as pd

from qiime2 import Metadata, Artifact
from q2_sidle._utils import (_setup_dask_client, 
                             _convert_generator_to_seq_block,
                             _convert_generator_to_delayed_seq_block, 
                             _convert_seq_block_to_dna_fasta_format,
                             degen_reps,
                             database_params,
                             )


def filter_to_alignment_discard(alignment: pd.DataFrame,
                                table: biom.Table,
                                max_mismatch: int=None,
                                discarded: bool=True,
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
    regions: str,
    regional_alignments: pd.DataFrame,
    regional_tables: biom.Table,
    ) -> Metadata:
    """
    Tracks where sample counts have gone.
    """
    print(regions)
    # Maps the regioional counts
    regional_counts = [
        _alignment_accounting(region, align, table)
        for (region, align, table) 
        in zip(*(regions, regional_alignments, regional_tables))
        ]
    total_table = regional_tables[0].copy().merge(*regional_tables[1:])
    total_counts = _alignment_accounting('total', 
                                         pd.concat(regional_alignments),
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
    return regional_account.add_prefix('%s ' % region)
