import re

import dask
import numpy as np
import pandas as pd
from skbio import DNA, TabularMSA
from skbio.alignment import local_pairwise_align_ssw

from qiime2 import Metadata
from q2_feature_classifier._cutter import (_local_aln)
from q2_sidle._utils import (_setup_dask_client, 
                             _check_manifest,
                             _read_manifest_files,
                             degen_sub,
                             )


def reconstruct_fragment_rep_seqs(
    reconstruction_map: pd.Series, 
    reconstruction_summary:pd.DataFrame,
    aligned_sequences:pd.Series,
    manifest:Metadata,
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
    # Checks the  manifest
    _check_manifest(manifest)
    region_order = manifest.get_column('region-order').to_series().astype(int)
    region_names = {i: r for r, i in region_order.items()}

    # Gets the total number of database sequences mapped in the summary
    reconstruction_summary['num_seqs_mapped'] = \
        pd.Series({i: i.count("|") + 1 for i in reconstruction_summary.index})

    # Filters down to only clusters that include multiple sequences 
    multiple_mapped = reconstruction_summary['num_seqs_mapped'] > 1
    multiple_mapped = multiple_mapped.index[multiple_mapped].values
    reconstruction_map = \
        reconstruction_map.loc[reconstruction_map.isin(multiple_mapped)].copy()
    reconstruction_summary = reconstruction_summary.loc[multiple_mapped]

    # Loads the mapping of the regional kmers, adds. the clean name
    # and tidies the information
    kmer_map = pd.concat(
        axis=0, 
        sort=False,
        objs=_read_manifest_files(manifest, 'kmer-map', 
                                 'FeatureData[KmerMap]', pd.DataFrame)
        )
    kmer_map = kmer_map.loc[reconstruction_map.index].copy()
    kmer_map['clean_name'] = reconstruction_map
    kmer_map.replace({'region': region_order}, inplace=True)
    kmer_map.reset_index(inplace=True)
    kmer_map['kmer-length'] = kmer_map['kmer-length'].astype(int)
    
    # Gets the first kmer and the associated information
    kmer_map.sort_values(by=['clean_name', 'db-seq', 'seq-name', 'region'], 
                         ascending=True, 
                         inplace=True)
    kmer_map.drop_duplicates(['clean_name', 'db-seq', 'region'], 
                             inplace=True)
    kmer_map.set_index('db-seq', inplace=True)
    kmer_map['sequence'] = aligned_sequences.loc[kmer_map.index]
    kmer_map.reset_index(inplace=True)

    # We'll drop kmer-specific info
    kmer_map.drop(columns=['seq-name', 'kmer'], 
                  inplace=True)

    # Calculates and degaps the concensus sequence for the region using skbio
    # to calculate the concensus and then fiters the concensus down to a 
    # single sequence
    concensus = \
        kmer_map.drop_duplicates('db-seq').copy().set_index('db-seq')
    concensus = \
        concensus.groupby('clean_name')['sequence'].apply(_group_concensus)


    concensus_cols = ['clean_name', 'fwd-primer', 
                      'rev-primer', 'kmer-length']
    concensus_map = kmer_map[concensus_cols].copy()
    concensus_map.drop_duplicates(inplace=True)
    concensus_map.set_index('clean_name', inplace=True)
    concensus_map['sequence'] = concensus

    # Determines the positions for each of the primer regions. We first look
    # for an exact match using regex
    fwd_rep = {primer: _expand_primer(primer, None) 
               for primer in concensus_map['fwd-primer'].unique()}
    rev_rep = {primer: _expand_primer(primer, None) 
               for primer in concensus_map['rev-primer'].unique()}
    concensus_map['f_exact'] = concensus_map['fwd-primer'].replace(fwd_rep)
    concensus_map['r_exact'] = concensus_map['rev-primer'].replace(fwd_rep)

    # Finds sequences with exact matches
    concensus_map['start'] = \
        concensus_map[['sequence', 'f_exact']].astype(str).apply(
        _find_exact_forward, axis=1
    )

    # Finds fuzzy matches. We assume that since we know sequences from this
    # concensus were amplified using this primer that we let the fuzzy
    # match go
    missing_start = concensus_map['start'].isna()
    concensus_map.loc[missing_start, 'start'] = \
    concensus_map.loc[missing_start, ['sequence', 'fwd-primer']].apply(
        _find_approx_forward, axis=1)

    # Adds the kmer offset
    concensus_map['start_off'] = \
        concensus_map['start'] + concensus_map['kmer-length']

    # Finds the most extreme positions based on the overlapping positions at
    # each region
    concensus_map['fwd-pos'] = \
        concensus_map[['start', 'start_off']].min(axis=1)
    concensus_map['rev-pos'] = \
        concensus_map[['start', 'start_off']].max(axis=1)
    
    # Gets that fragment based on the most extreme forward and reverse 
    # positions in the amplified regions
    concensus_map.reset_index(inplace=True)
    concensus_map['sequence'] = concensus_map['sequence'].astype(str)

    con_group = concensus_map.groupby(['clean_name', 'sequence'])
    concensus_sub = pd.concat(axis=1, objs=[
        con_group['fwd-pos'].min().astype(int), 
        con_group['rev-pos'].max().astype(int)
    ]).astype(int).reset_index(level='sequence')
    concensus_sub = concensus_sub.apply(
        lambda x: DNA(x['sequence'][x['fwd-pos']:x['rev-pos']]), axis=1)

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
