import re

import dask
import numpy as np
import pandas as pd
from skbio import DNA, TabularMSA
from skbio.alignment import local_pairwise_align_ssw

from qiime2 import Metadata
from q2_feature_classifier._cutter import (_local_aln)
from q2_sidle._utils import (degen_sub,
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
    trim_to_fragment: bool, optional
        During alignmnt, if multiple sequences are mapped to the same region
        with one having shorter fragments, when `trim_to_fragment` is true,
        the alignmnent will be handled using that shorter sequence. When 
        False, the fragment will be reconstructed using the full length
        sequences for the joined table.
    gap_handling: str {'keep_nt', 'drop'}
        Whether gaps in the middle of alignments should inheriet the defined
        sequence or simply drop out.
    primer_mismatch: int, optional
        The number of mismatches between the primer and sequence allowed.
    debug: bool
        Whether the function should be run in debug mode (without a client)
        or not. `debug` superceeds all options
    n_workers: int, optional
        The number of jobs to initiate. When `n_workers` is 0, the cluster 
        will be able to access all avalaibel resources.

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
    reconstruction_map = \
        reconstruction_map.loc[reconstruction_map['clean_name'].isin(
            multiple_mapped)].copy()

    # print(aligned_sequences)
    reconstruction_map['sequence'] = \
        aligned_sequences.loc[reconstruction_map.index]
    reconstruction_map.reset_index(inplace=True)
    reconstruction_map.sort_values(['clean_name', 'first-region'], 
                                   inplace=True)

    # # print(reconstruction_map)

    # Calculates and degaps the concensus sequence for the region using skbio
    # to calculate the concensus and then fiters the concensus down to a 
    # single sequence
    concensus_map = pd.concat(axis=1, objs=[
        reconstruction_map.groupby('clean_name', sort=False)['sequence'].apply(
            _group_concensus),
        reconstruction_map.groupby('clean_name', sort=False)[
        'first-fwd-primer'].first(),
        reconstruction_map.sort_values('last-region').groupby('clean_name', 
            sort=False)[['last-fwd-primer', 'last-kmer-length']].last()
        ])

    # Determines the positions for each of the primer regions. We first look
    # for an exact match using regex
    fwd_rep = {primer: _expand_primer(primer, None) 
               for primer in concensus_map['first-fwd-primer'].unique()}
    rev_rep = {primer: _expand_primer(primer, None) 
               for primer in concensus_map['last-fwd-primer'].unique()}
    concensus_map['f_exact'] = \
        concensus_map['first-fwd-primer'].replace(fwd_rep)
    concensus_map['r_exact'] = \
        concensus_map['last-fwd-primer'].replace(rev_rep)

    # Finds sequences with exact matches
    concensus_map['start'] = \
        concensus_map[['sequence', 'f_exact']].astype(str).apply(
        _find_exact_forward, axis=1
    )
    concensus_map['end'] = \
        concensus_map[['sequence', 'r_exact']].astype(str).apply(
        _find_exact_forward, axis=1
    )

    # Finds fuzzy matches. We assume that since we know sequences from this
    # concensus were amplified using this primer that we let the fuzzy
    # match go
    missing_start = concensus_map['start'].isna()
    concensus_map.loc[missing_start, 'start'] = \
        concensus_map.loc[missing_start, ['sequence', 'first-fwd-primer']
        ].apply(_find_approx_forward, axis=1)
    missing_end = concensus_map['end'].isna()
    concensus_map.loc[missing_start, 'end'] = \
        concensus_map.loc[missing_start, ['sequence', 'last-fwd-primer']
        ].apply(_find_approx_forward, axis=1)
    concensus_map['end'] = \
        concensus_map[['end', 'last-kmer-length']].sum(axis=1)

    concensus_map[['start', 'end']] = \
        concensus_map[['start', 'end']].astype(int)
    
    # Gets that fragment based on the most extreme forward and reverse 
    # positions in the amplified regions
    concensus_sub = concensus_map.apply(
        lambda x: DNA(x['sequence'][x['start']:x['end']]), axis=1)

    return concensus_sub


def _expand_primer(primer, error=None):
    """
    Expands degenerates in a primer to allow alignment
    """
    if error is not None:
        new_ = '((?<=(%s)))((%s){e<=%s})' % (primer[0], primer[1:], error)
    else:
        new_ = primer
    for bp, rep in degen_sub.items():
        if bp in {'A', 'T', 'G', 'C'}:
            continue
        new_ = new_.replace(bp, '[%s]' % rep)

    return new_


def _find_exact_forward(args):
    """
    Finds the exact match for a forward primer
    """
    [sequence, primer] = args.values.flatten()
    match = re.search(primer, sequence)
    if match is None:
        return np.nan
    else:
        return match.span()[1]
    

def _find_exact_reverse(args):
    """
    Finds the exact match for a reverse primer
    """
    [sequence, primer] = args.values.flatten()
    match = re.search(primer, sequence)
    if match is None:
        return np.nan
    else:
        return match.span()[0]


def _find_approx_forward(args):
    """
    Finds an approximate match for a forward primer
    """
    [sequence, primer] = args
    primer = DNA(primer)
    align_ = _local_aln(sequence, primer)
    return align_[2][0][1] + 1


def _find_approx_reverse(args):
    """
    Finds an approximate match for a reverse primer
    """
    [sequence, primer] = args
    primer = DNA(primer)
    align_ = _local_aln(sequence, primer)
    return align_[2][0][0]


def _group_concensus(seqs):
    """
    Collapses grouped sequences to a concensus sequence
    """
    concensus = TabularMSA(seqs).consensus().degap()
    return concensus
