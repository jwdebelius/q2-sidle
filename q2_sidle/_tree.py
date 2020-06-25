from functools import partial

import dask
import numpy as np
import pandas as pd
from skbio import DNA

from qiime2 import Metadata
from q2_types.feature_data import DNAFASTAFormat
from q2_sidle._utils import (_setup_dask_client, 
                             _find_primer_end,
                             _find_primer_start,
                             _check_manifest,
                             _read_manifest_files,
                             degen_sub,
                             degen_sub2,
                             degen_undo,
                             )
from q2_sidle._extract import _trim_masked


def reconstruct_fragment_rep_seqs(reconstruct_map:pd.Series, 
    reconstruction_summary:pd.DataFrame,
    reconstruction_manifest:Metadata,
    aligned_sequences:pd.Series,
    trim_to_fragment:bool=True,
    gap_handling:str='keep',
    primer_mismatch:int=2,
    debug: bool=False, 
    n_workers: int=0
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
    reconstruct_map = \
        reconstruct_map.loc[reconstruct_map.isin(multiple_mapped)].copy()
    reconstruction_summary = reconstruction_summary.loc[multiple_mapped]

    # Loads the mapping of the regional kmers
    kmer_map = pd.concat(
        axis=0, 
        sort=False,
        objs=_read_manifest_files(manifest, 'kmer-map', 
                                 'FeatureData[KmerMap]', pd.DataFrame)
        )
    kmer_map = kmer_map.loc[reconstruct_map.index]
    kmer_map.sort_values(by=['clean_name', 'db-seq', 'seq-name', 'region'], 
                         ascending=True, 
                         inplace=True)
    kmer_map.drop_duplicates(['clean_name', 'db-seq', 'region'], inplace=True)
    kmer_map.drop(columns=['seq-name', 'kmer', 'region'], inplace=True)

    # Updates the kmer map with information about coverage and kmer length
    kmer_map['clean_name'] = reconstruct_map
    kmer_map.reset_index(inplace=True)
    kmer_map.set_index('clean_name', inplace=True)
    kmer_map['final_coverage'] = reconstruction_summary['num-regions']
    kmer_map.reset_index(inplace=True)
    kmer_map['region'] = kmer_map['region'].replace(region_order)
    kmer_map['kmer-length'] = kmer_map['kmer-length'].astype(int)
    kmer_map.set_index('db-seq', inplace=True)


    # Switches kmer primers to expanded primers
    fwd_rep = {primer: _expand_primer(primer, primer_mismatch) 
               for primer in kmer_map['fwd-primer'].unique()}
    rev_rep = {primer: _expand_primer(primer, primer_mismatch) 
               for primer in kmer_map['rev-primer'].unique()}
    kmer_map.replace({'fwd-primer': fwd_rep, 'rev-primer': rev_rep}, 
                     inplace=True)

    # Adds in the full length aligned sequence
    kmer_map['sequence'] = aligned_sequences[kmer_map.index]

    # Splits the data based on whether it covers a single region or 
    # multiple regions
    single_region = kmer_map.loc[kmer_map['final_coverage'] == 1].copy()
    multiple_regions = kmer_map.loc[kmer_map['final_coverage'] > 1].copy()

    # Sets up the function for multiple mapping
    if trim_to_fragment:
        collapse_f = partial(_collapse_expanded_trim,
                             gap_management=gap_handling)
    else:
        collapse_f = partial(_collapse_expanded_all, 
                             gap_management=gap_handling)
    fragment = pd.concat(axis=0, objs=[
        _find_single_region_fragments(single_region),
        _find_multiple_region_fragments(multiple_regions),
        ])
    fragment = pd.Series({id_: DNA(seq_, metadata={'id': id_})
                          for id_, seq_ in fragment})
    return fragment


def _expand_primer(primer, error):
    """
    Expands degenerates in a primer to allow alignment
    """
    new_ = '((?<=(%s)))((%s){e<=%s})' % (primer[0], primer[1:], error)
    for bp, rep in degen_sub.items():
        if bp in {'A', 'T', 'G', 'C'}:
            continue
        new_ = new_.replace(bp, '[%s]' % rep)

    return new_


def _collapse_expanded(id_, g, gap_management='drop', trim=True, 
                       sub=degen_undo):
    """
    Collapses two squences using the full length sequences
    
    Parameters
    ----------
    id_: str
        The name of the collapsed sequence
    g: DataFrame
        information for all the sequences associated with the 
        collapsed sequence id
    gap_management: str
        The way gaps in the aligmnet should handled. "keep_nt" will
        ignore the gap and use the nucleotides from the original
        sequnce. "drop" will remove all gaps that don't appear at the
        begining or end of an alignment. Any other string will be 
        used to replace the gap.
    trim: bool
        Whether shorter sequences should trimmed before alignment (True) or 
        not (False)
    sub: dict
        A look up to convert degenerate nuclotides to a single letter 
        code
    
    Returns
    -------
    pd.Series
        The collapsed aligned sequence
    """
    expanded2 = g['sequence'].apply(lambda x: pd.Series(list(x)))
    expanded2.replace(degen_sub2, inplace=True)
    
    # Expands the list of sequnces
    nt_pos = expanded2.columns.values + 1
    num_nt = len(nt_pos)
    num_seq = len(expanded2)
    position = np.vstack([nt_pos] * num_seq)

    if trim:
        fwd_ = g[['fwd-pos']].values * np.ones((1, num_nt))
        rev_ = g[['rev-pos']].values * np.ones((1, num_nt))
    else:
        fwd_ = g[['group-fwd']].values * np.ones((1, num_nt))
        rev_ = g[['group-rev']].values * np.ones((1, num_nt))
    trim_mask = (position < (fwd_ + 1)) | (position > (rev_))

    # Covers skips at the begining or end of the sequence
    start_check = (expanded2.isin(['-', ''])).cumsum(axis=1)
    start_check = start_check == position
    end_check = \
        expanded2.isin(['-', ''])[expanded2.columns[::-1]].cumsum(axis=1)
    end_check = end_check[expanded2.columns] == (num_nt + 1 - position)
    flank_check = start_check | end_check
    expanded2.mask((flank_check | trim_mask), '', inplace=True)

    skips = expanded2 == '-'
    # Handles gaps in the sequneces
    if gap_management == 'keep':
        expanded2.mask(skips, '', inplace=True)
    elif gap_management == 'drop':
        expanded2.mask(skips, np.nan, inplace=True)
    else:
        expanded2.replace({"-", gap_management}, inplace=True)
    expanded2.dropna(axis=1, inplace=True, how='any')

    # Replaces the seequence and adds the degenerates
    def _fix_nt(x):
        x1 = sorted(''.join(x).replace('[', '').replace(']', ''))
        x2 = np.unique(x1)
        return ''.join(x2)
    collapse = expanded2.apply(_fix_nt)
    collapse.replace(sub, inplace=True)

    return pd.Series({id_: ''.join(collapse.values)})


def _find_positions(x):
    [sequence, fwd_primer, rev_primer, kmer_length] = \
        x[['sequence', 'fwd-primer', 'rev-primer', 'kmer-length']]
    if kmer_length > 0:
        fwd_pos = _find_primer_end(sequence, fwd_primer)['pos']
        rev_pos = fwd_pos + kmer_length - 1
    else:
        rev_pos = _find_primer_start(sequence, rev_primer)['pos']
        fwd_pos = rev_pos - kmer_length - 1
    return pd.Series({"fwd_pos": fwd_pos, 'rev_pos': rev_pos})


def _find_multiple_region_fragments(seq_map, gap_management='keep', 
    trim=True):
    """
    Collapses sequence fragments which cover multiple regions
    """
    # Gets the expanded sequences and tidies the gapped alignments
    expand = seq_map.set_index('clean_name', append=True)['sequence'].apply(
        lambda x: pd.Series(list(x)))
    expand.replace({'-': np.nan}, inplace=True)
    expand.dropna(how='all', axis=1, inplace=True)
    expand.reset_index(level='clean_name', inplace=True)

    # Gets the sub alignments associated with the sequence cluster
    def _get_sub_aligned(g):
        g.dropna(axis=1, how='all', inplace=True)
        g.fillna('-', inplace=True)
        g.drop(columns=['clean_name'], inplace=True)
        return g.apply(lambda x: ''.join(x), axis=1)
    seq_map['sequence'] = pd.concat([
        _get_sub_aligned(g.copy()) for _, g in expand.groupby('clean_name')
        ])

    # Finds the primer positions for each region
    seq_map[['fwd-pos', 'rev-pos']] = seq_map.apply(_find_positions, axis=1)

    # Gets the individual and regional limits
    seq_map.reset_index(inplace=True)
    seq_map.sort_values(['db-seq', 'fwd-pos'], ascending=True, inplace=True)
    seq_map['fwd-pos'] = seq_map.groupby('db-seq')['fwd-pos'].cummin()
    seq_map.sort_values(['db-seq', 'rev-pos'], ascending=False, inplace=True)
    seq_map['rev-pos'] = seq_map.groupby('db-seq')['rev-pos'].cummax()

    # Gets the grouped starting position
    seq_map.sort_values(['clean_name', 'rev-pos'], ascending=False, 
                        inplace=True)
    seq_map['group-rev'] = seq_map.groupby('clean_name')['rev-pos'].cummax()
    seq_map.sort_values(['clean_name', 'fwd-pos'], ascending=True, 
                        inplace=True)
    seq_map['group-fwd'] = seq_map.groupby('clean_name')['fwd-pos'].cummin()

    # Gets the grouped sequence
    fragment = pd.concat(axis=0, objs=[
        _collapse_expanded(id_, g, gap_management, trim) 
        for id_, g in seq_map.groupby('clean_name')
        ])
    fragment.name = 'fragment'
    fragment.index.set_names('clean_name', inplace=True)

    return fragment.sort_index(ascending=True)


def _find_single_region_fragments(seq_map):
    """
    Collapes sequences that cover a single region into a fragment over that 
    region
    """
    # Because we only have a single region, all the kmers in the region are 
    # the same and therefore we only have to work with it once
    seq_map.drop_duplicates(['clean_name'], inplace=True)
    seq_map.set_index('clean_name', inplace=True)
    
    # We return to an unaligned sequence   
    seq_map['sequence'] = \
        seq_map['sequence'].apply(lambda x: x.replace("-", ''))

    # Extracts the sequence
    seq_map[['fwd-pos', 'rev-pos']] = seq_map.apply(_find_positions, axis=1)
    seq_map['fragment'] = _trim_masked(
        seq_map[['sequence']], 
        seq_map, 
        start_col='fwd-pos', 
        end_col='rev-pos'
        ).apply(lambda x: ''.join(x), axis=1)
    # print(seq_map)

    return seq_map['fragment']


def _degen_replacement(sequence):
    seq = pd.Series(list(x))
    seq.replace(degen_sub2, inplace=True)
    return ''.join(x)

