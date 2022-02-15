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
    # num_search: int=10,
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
    # Temp arg TODO: Move into pluginb setup

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

    num_search = int(np.max([10, np.ceil(len(target_seqs) * 0.1)]))


    # If the positions are already defined, life is fantastic. Otherwise, 
    # we need to find them.
    no_pos = pd.isnull(sub_map[['first-fwd-pos', 'last-fwd-pos']])

    ### We don't need to map any positions. Whoo!
    if (no_pos == False).all().all():
        sub_map['start'] = sub_map['first-fwd-pos']
        sub_map['end'] = sub_map['last-fwd-pos'] + sub_map['last-kmer-length']

    ### We need to map all the positions
    else:
        # Gets the primers to map
        primer_seqs = pd.concat(axis=0, objs=[
            sub_map['first-fwd-primer'], sub_map['last-fwd-primer']
            ])
        primer_seqs.index.set_names('db-seq', inplace=True)
        primer_seqs.name = 'primer'
        primer_seqs = primer_seqs.reset_index()
        primer_seqs.drop_duplicates(inplace=True)

        # Extracts the primer positions from the alignment
        primer_pos = dict()
        for primer, sub_ids in primer_seqs.groupby('primer')['db-seq']:
            # Extracts the sequences and the primer
            exp_primer = _expand_primer(primer)
            sub_seqs = target_seqs[sub_ids[:num_search]].astype(str).copy()

            # Checks for an exact match
            exact_match = \
                sub_seqs.apply(_find_exact_forward, primer=exp_primer)
            position = exact_match.dropna().unique()[0]
            primer_pos[primer] = position

        # Updates the data map with the positions
        sub_map['start'] = sub_map['first-fwd-primer'].replace(primer_pos)
        sub_map['end'] = sub_map['last-fwd-primer'].replace(primer_pos) + \
            sub_map['last-kmer-length']

    # Extracts the fragments
    sub_map[['start', 'end']] = sub_map[['start', 'end']].astype(int)
    sub_map['seq'] = target_seqs

    concensus = pd.DataFrame({
        'seq': sub_map.groupby('clean_name')['seq'].apply(_group_concensus),
        'start': sub_map.groupby('clean_name')['start'].min(),
        'end': sub_map.groupby('clean_name')['end'].max(),
        })
    concensus['fragment'] = \
        concensus.apply(lambda x: x['seq'][x['start']:x['end']], axis=1)

    concensus_sub = pd.Series({
        name: DNA(str(seq), metadata={'id': name})
        for name, seq in concensus['fragment'].items()
        })
    concensus_sub.index.set_names('clean_name', inplace=True)

    return concensus_sub


def _expand_primer(primer):
    """
    Expands degenerates in a primer to allow alignment
    """
    for bp, rep in degen_sub2.items():
        primer = primer.replace(bp, rep)

    return primer


def _find_exact_forward(sequence, primer):
    """
    Finds the exact match for a forward primer
    """
    exp = _expand_primer(primer)
    match = re.search(exp, sequence)
    if match is None:
        return np.nan
    else:
        return match.span()[1]

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
