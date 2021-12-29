import os
import warnings

import biom
import numpy as np
import pandas as pd
import skbio

from qiime2 import Metadata


def reverse_complement_sequence(sequence: pd.Series) -> pd.Series:
    """
    Reverse complements the set of sequences
    """
    reverse_complement = sequence.apply(lambda x: x.reverse_complement())
    
    return reverse_complement


def find_first_alignment_position(alignment: pd.Series, 
                                  representative_sequences: pd.Series,
                                  table: biom.Table=None,
                                  direction: str='fwd',
                                  ) -> Metadata:
    """
    Finds the starting position for each sequence
    """
    ids_both = representative_sequences.index
    if set(ids_both).issubset(set(alignment.index)) == False:
        raise ValueError('All the representative sequences must be present '
                         'in the alignment. Please make sure that you have '
                         'the correct representative sequences and '
                         'alignment.')
    alignment = alignment.loc[ids_both].astype(str)
    if table is not None:
        table = table.filter(ids_both, axis='observation')

    # Calculates the places the sequences cover
    coverage = _generate_align_mask(alignment, ids_both)

    # Identifies the starting position for the alignment
    first_pos = _get_first_align_pos(coverage, table)
    if direction == 'rev':
        first_pos['starting-position'] = \
            len(coverage.columns) - first_pos['starting-position']

    # Cleans up the result to be presented as metadata to make filtering 
    # and alignment easier
    first_pos['direction'] = direction
    first_pos['starting-position'] = \
        first_pos['starting-position'].astype(int).astype(str)
    first_pos.index.set_names('feature-id', inplace=True)

    return Metadata(first_pos)


def find_span_positions(alignment: pd.Series, 
                        representative_sequences: pd.Series,
                        region_name: str='region',
                        ) -> Metadata:
    """
    Finds the minimum and maximum alignment positions where the
    reference sequences are within the alignment
    """
    ids_both = representative_sequences.index
    if set(ids_both).issubset(set(alignment.index)) == False:
        raise ValueError('All the representative sequences must be present '
                         'in the alignment. Please make sure that you have '
                         'the correct representative sequences and '
                         'alignment.')
    alignment = alignment.loc[ids_both].astype(str)
    
    coverage = _generate_align_mask(alignment, ids_both)
    coverage.replace({0: np.nan}, inplace=True)
    coverage.dropna(how='all', axis=1, inplace=True)
    left = coverage.columns.min()
    right = coverage.columns.max()

    summary = pd.DataFrame(
        data=np.array([[left, right]], dtype=float),
        index=pd.Index([region_name], name='id'),
        columns=['left', 'right'],
        )

    return Metadata(summary)


def _generate_align_mask(alignment, seq_ids):
    """
    Gets the masked positions
    """
    expanded = alignment.loc[seq_ids].astype(str)
    expanded = expanded.apply(lambda x: pd.Series(list(x)))
    expanded.replace({'-': np.nan}, inplace=True)
    coverage = (expanded.notna() * 1)

    return coverage


def _get_first_align_pos(coverage, table=None, min_median_freq=100):
    """
    Finds the starting position for the alignment
    """
    # Gets the first position in hte coverage matrix
    pos = np.vstack([coverage.columns.values] * len(coverage)) + 1
    first_pos = ((coverage.cumsum(axis=1) == 1) * pos).max(axis=1)
    first_pos = first_pos[first_pos > 0] - 1
    first_pos.index.set_names('feature-id', inplace=True)
    first_pos = pd.DataFrame(first_pos, columns=['starting-position'])
    
    # Appends the first position with the sequenced counts
    if table is not None:
        first_pos['sequence-counts'] = \
            pd.Series(table.sum(axis='observation'),
                      index=list(table.ids(axis='observation')),
                      )
    else:
        first_pos['sequence-counts'] = min_median_freq + 1

    return first_pos
