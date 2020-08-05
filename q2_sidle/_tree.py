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
    
    # We assume that all the sequences here amplified so we're going to ignore
    # the reverse primer and just barrel forward from there.
    kmer_map.drop(columns=['seq-name', 'kmer', 'rev-primer'], 
                  inplace=True)

    # Finds the start and end position for each region
    fwd_rep = {primer: _expand_primer(primer, None) 
               for primer in kmer_map['fwd-primer'].unique()}
    kmer_map['exact'] = kmer_map['fwd-primer'].replace(fwd_rep)
    kmer_map['start'] = kmer_map[['sequence', 'exact']].astype(str).apply(
        _find_exact_forward, axis=1) 

    missing_start = kmer_map['start'].isna()
    kmer_map.loc[missing_start, 'start'] = \
        kmer_map.loc[missing_start, ['sequence', 'fwd-primer']].apply(
            _find_approx_forward, axis=1)

    kmer_map['end'] = kmer_map['start'] + kmer_map['kmer-length']

    # Since the forward and reverse positions may be flipped, we pull out 
    # the forward and reverse only
    kmer_map['fwd-pos'] = kmer_map[['start', 'end']].min(axis=1)
    kmer_map['rev-pos'] = kmer_map[['start', 'end']].max(axis=1)

    # Gets the concensus sequence
    concensus = kmer_map.groupby('clean_name').apply(_get_concensus_seq)

    return concensus


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
    Finds the exact match for sequences
    """
    [sequence, primer] = args.values.flatten()
    match = re.search(primer, sequence)
    if match is None:
        return np.nan
    else:
        return match.span()[1]

def _find_approx_forward(args):
    [sequence, primer] = args
    primer = DNA(primer)
    align_ = _local_aln(sequence, primer)
    return align_[2][0][1] + 1


def _get_concensus_seq(g):
    """
    Gets the concensus sequence for a multiple sequence alignment
    """
    fwd = int(g['fwd-pos'].min())
    rev = int(g['rev-pos'].max())
    g['sequence'] = g['sequence'].apply(lambda x: DNA(x[fwd:rev]))
    align = TabularMSA(g['sequence'].values)
    concensus = align.consensus().degap()
    
    return concensus
