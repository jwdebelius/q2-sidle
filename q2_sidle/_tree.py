import re

import dask
import numpy as np
import pandas as pd
from skbio import DNA, TabularMSA
from skbio.alignment import local_pairwise_align_ssw

from qiime2 import Metadata
from q2_feature_classifier._cutter import (_local_aln)
from q2_sidle._utils import (degen_sub2,
                             _check_regions,
                            )


def reconstruct_fragment_rep_seqs(
    reconstruction_map: pd.DataFrame, 
    reconstruction_summary:pd.DataFrame,
    aligned_sequences:pd.Series,
    ) -> pd.Series:
    """
    Assembles fragments from assemblies that include multiple sequences

    Parameters
    ----------
    reconstruction_map: Series
        A map between the final kmer name and the original database sequence.
    reconstruction_summary: DataFrame
        A summary of the statitics for the regional map describing the number
        of regions mapped to each reference sequence and the number of kmers.
        The kmer mapping estimate can account for degeneracy when the 
        `count_degenerates` flag is used or can ignore degenrate sequences in
        mapping
    manifest: Metadata
         A manifest file describing the relationship between regions and their
        alignment mapping. The manifest must have at least three columns
        (`kmer-map`, `alignment-map` and `frequency-table`) each of which
        contains a unique filepath. The regions should be sorted and the 
        region labels must match between the kmer map and the alignment map.
    aligned_sequences: Series
        The aligned representative sequences for the data

    Returns
    -------
    pd.Series
        The sequence fragments to be used for fragment insertion and multiple
        sequence alignment

    """
    # Gets the total number of database sequences mapped in the summary
    reconstruction_summary['num_seqs_mapped'] = \
        pd.Series({i: i.count("|") + 1 for i in reconstruction_summary.index})
    if (reconstruction_summary['num_seqs_mapped'] == 1).all():
        return pd.Series({'1': DNA('ATCG')})

    # Filters down to only clusters that include multiple sequences 
    multiple_mapped = reconstruction_summary['num_seqs_mapped'] > 1
    multiple_mapped = multiple_mapped.index[multiple_mapped].values
    sub_map = \
        reconstruction_map.loc[reconstruction_map['clean_name'].isin(
            multiple_mapped)].copy()

    target_seqs = aligned_sequences.loc[sub_map.index].copy()
 
    ### Gets primer positions in the alignment
    # Finds primer sequences and positions if the sequence positions haven't
    # been identified
    first_ = sub_map[['first-fwd-pos', 'first-fwd-primer']].copy()
    last_ = sub_map[['last-fwd-pos', 'last-fwd-primer']].copy()
    primers = pd.concat(axis=0, objs=[
        primer.rename(columns={c: c.split('-')[-1] for c in primer.columns})
        for primer in [first_, last_]
        ])
    primers.drop_duplicates(inplace=True)
    primers.sort_values(by=['pos'], ascending=False, inplace=True)
    primers.drop_duplicates('primer', inplace=True)
    primers = primers.reset_index().set_index('primer')

    # Finds the primer position in the alignment
    for primer, [seq_name, pos] in primers.iterrows():
        if pd.isnull(pos):
            seq_ = target_seqs.loc[seq_name]
            pos = _find_approx_forward([seq_, primer])
            primers.loc[primer, 'pos'] = pos
    pos_map = primers['pos'].to_dict()
    print(pos_map)

    # Sets the position for all the primers and finds the span positions
    miss_1st = sub_map['first-fwd-pos'].isna()
    miss_end = sub_map['last-fwd-pos'].isna()

    sub_map.loc[miss_1st, 'first-fwd-pos'] = \
        sub_map.loc[miss_1st, 'first-fwd-primer'].replace(pos_map)
    sub_map.loc[miss_end, 'last-fwd-pos'] = \
        sub_map.loc[miss_end, 'last-fwd-primer'].replace(pos_map)

    sub_map['start'] = sub_map['first-fwd-pos']
    sub_map['end'] = sub_map['last-fwd-pos'] + sub_map['last-kmer-length']

    # Extracts the fragments
    sub_map['seq'] = target_seqs

    concensus = pd.DataFrame({
        'seq': sub_map.groupby('clean_name')['seq'].apply(_group_concensus),
        'start': sub_map.groupby('clean_name')['start'].min(),
        'end': sub_map.groupby('clean_name')['end'].max(),
        })
    concensus_sub = concensus.apply(
        lambda x: x['seq'][x['start']:x['end']], axis=1
        )

    return concensus_sub


def _find_approx_forward(args):
    """
    Finds an approximate match for a forward primer
    """
    [sequence, primer] = args
    primer = DNA(primer)
    align_ = _local_aln(sequence, primer)
    return align_[2][0][1] + 1


def _group_concensus(seqs):
    """
    Collapses grouped sequences to a concensus sequence
    """
    concensus = TabularMSA(seqs).consensus().degap()
    return concensus
