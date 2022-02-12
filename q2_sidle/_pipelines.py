import copy
import re

import numpy as np
import pandas as pd
from skbio import DNA

from qiime2 import Metadata

from q2_sidle._utils import degen, degenerate_map

def sidle_reconstruction(ctx, 
                         region, 
                         kmer_map, 
                         regional_alignment, 
                         regional_table, 
                         reference_taxonomy, 
                         min_counts=1000,
                         database='none', 
                         define_missing='merge',
                         block_size=10000,
                         debug=False, 
                         n_workers=1,
                         client_address=None,
                         ):
    """
    A pipeline to reconstruct the full data
    """
    reconstruct_database = ctx.get_action('sidle', 'reconstruct_database')
    reconstruct_counts = ctx.get_action('sidle', 'reconstruct_counts')
    reconstruct_taxonomy = ctx.get_action('sidle', 'reconstruct_taxonomy')
    reconstruct_fragment = \
        ctx.get_action('sidle', 'reconstruct_fragment_rep_seqs')

    print("Reconstructing database")
    db_map, db_summary = reconstruct_database(
        region=region,
        kmer_map=kmer_map,
        regional_alignment=regional_alignment,
        count_degenerates=True,
        block_size=block_size,
        n_workers=n_workers,
        client_address=client_address,
        debug=debug,
        )

    print("Reconstructing counts")
    counts, = reconstruct_counts(
        region=region,
        regional_alignment=regional_alignment,
        regional_table=regional_table,
        database_map=db_map,
        database_summary=db_summary,
        region_normalize='average',
        min_counts=min_counts,
        block_size=block_size,
        n_workers=n_workers,
        client_address=client_address,
        debug=debug,
        )

    print("Reconstructing Taxonomy")
    taxonomy, = reconstruct_taxonomy(
        reconstruction_map=db_map,
        taxonomy=reference_taxonomy,
        database=database,
        define_missing=define_missing,
        )

    results = [db_map, db_summary, counts, taxonomy]

    return tuple(results)


def reconstruct_tree(ctx,
                     reconstruction_summary,
                     reconstruction_map,
                     aligned_sequences,
                     sepp_reference_database,
                     n_threads=1,
                     ):
    """A pipeline for tree building, because why not?
    """
    reconstruct_fragments = \
        ctx.get_action('sidle', 'reconstruct_fragment_rep_seqs')
    insert_fragments = ctx.get_action('fragment_insertion', 'sepp')
    prune_tree = ctx.get_action('phylogeny', 'filter_tree')
    
    print('Reconstructing Fragments')
    rep_fragments, = reconstruct_fragments(
        reconstruction_summary=reconstruction_summary,
        reconstruction_map=reconstruction_map,
        aligned_sequences=aligned_sequences,
        )

    print('Reconstructing Tree')
    tree, placements = insert_fragments(
        representative_sequences=rep_fragments,
        referencce_database=sepp_reference_database,
        n_threads=n_threads,
        )

    return (rep_fragments, tree, placements)


def map_alignment_positions(ctx, 
                            alignment,
                            sequences,
                            direction,
                            table=None,
                            reverse_complement_ref=True,
                            n_threads=1,
                            add_fragments=False,
                            colormap=None,
                            ):
    rc_seqs = \
        ctx.get_action('sidle', 'reverse_complement_sequence')
    expand_alignment = \
        ctx.get_action('alignment', 'mafft_add')
    first_position = \
        ctx.get_action('sidle', 'find_first_alignment_position')
    ballet_recital = \
        ctx.get_action('sidle', 'summarize_alignment_positions')




    if (direction == 'rev') & reverse_complement_ref:
        alignment, = rc_seqs(alignment)

    expanded, = expand_alignment(alignment=alignment,
                                 sequences=sequences,
                                 addfragments=add_fragments,
                                 n_threads=n_threads,
                                 )
    print('alignment expanded')

    starts, = first_position(alignment=expanded,
                             representative_sequences=sequences,
                             table=table,
                             direction=direction,
                             )
    print('got starting position')

    viz, = ballet_recital(
        alignment=expanded,
        position_summary=starts,
        colormap=colormap,
        weight_by_abundance=(table is not None),
        )

    return (expanded, starts, viz)


def find_and_prepare_regional_seqs(ctx,
                                   alignment,
                                   region,
                                   amplicons=None,
                                   fwd_primer=None,
                                   rev_primer=None,
                                   length=None,
                                   subset_size=0.1,
                                   subset_seed=None,
                                   add_fragments=False,
                                   reverse_complement_rev_primer=True,
                                   chunk_size=10000, 
                                   debug=False, 
                                   n_workers=1,
                                   client_address=None,
                                   ):
    subsample_alignment = ctx.get_action('rescript', 'subsample_fasta')
    expand_alignment = ctx.get_action('alignment', 'mafft_add')
    get_span = ctx.get_action('sidle', 'find_alignment_span_positions')
    extract_positions = ctx.get_action('rescript', 'trim_alignment')
    degap_seqs = ctx.get_action('rescript', 'degap_seqs')
    filter_seqs = ctx.get_action('rescript', 'filter_seqs_length')
    prep_seqs = ctx.get_action('sidle', 'prepare_extracted_region')

    if pd.isnull([amplicons, fwd_primer, rev_primer, length]).all():
        raise ValueError('You must provide sequences, or a forward'
                         ' and reverse primer pair.')
    elif (pd.isnull([fwd_primer, rev_primer, length]).sum() in {0, 3}) == False:
        raise ValueError("You must either specify both a forward "
                         "primer and a reverse primer, or skip "
                         "the primer pair.")
    elif pd.isnull(amplicons) == pd.isnull(fwd_primer):
        raise ValueError('This function accepts a pair of primers and length'
                         ' or a set of references sequences. '
                         'Please only provide one.')

    # Determines if we should use the amplicons or if we should just trim the primers
    trim_primers = pd.isnull(amplicons)

    if (trim_primers == False):
        rep_seq_lens = \
            amplicons.view(pd.Series).astype(str).apply(lambda x: len(x))
        rep_seq_lens = rep_seq_lens.value_counts()
        if len(rep_seq_lens.index) > 1:
            raise ValueError('The representative sequences are of multiple '
                             'lengths. Please trim the representative sequences'
                             ' to a consistent length before using them for '
                             'extraction.')
        else:
            length = rep_seq_lens.index[0]

        print('Subsampling reference alignment')
        if subset_seed is None:
            subset_seed = np.random.randint(0, 9999, 1)
        small_align, = subsample_alignment(sequences=alignment,
                                           random_seed=subset_seed,
                                           subsample_size=subset_size,
                                           )

        print('Aligning representative sequences to reference')
        expanded, = expand_alignment(alignment=small_align,
                                     sequences=amplicons,
                                     addfragments=add_fragments,
                                     n_threads=n_workers,
                                     )
     
        print('Finding alignment positions')
        span_summary, = get_span(alignment=expanded, 
                                 representative_sequences=amplicons, 
                                 region_name=region
                                 )
        span2 =  span_summary.view(Metadata).to_dataframe().astype(int)
        [[left, right]] = span2.values

        print('Extracting and degapping reference sequences')
        sub_ref_align, = extract_positions(aligned_sequences=alignment,
                                           position_start=left,
                                           position_end=right,
                                           )
    else:
        if reverse_complement_rev_primer:
            rev_primer = str(DNA(rev_primer).reverse_complement())

        sub_ref_align, = extract_positions(aligned_sequences=alignment,
                                           primer_fwd=fwd_primer,
                                           primer_rev=rev_primer,
                                           )
        [left, right] = [None, None]

    sub_ref_seq, = degap_seqs(aligned_sequences=sub_ref_align)
    sub_ref_len, sub_ref_discard = \
        filter_seqs(sequences=sub_ref_seq, 
                    global_min=int(np.floor(length * 0.75)), 
                   )

    print('Preparing regional sidle database')
    kmer_seqs, kmer_map = \
        prep_seqs(sequences=sub_ref_len,
                 trim_length=length,
                 region=region,
                 fwd_primer=fwd_primer,
                 rev_primer=rev_primer,
                 fwd_pos=left,
                 rev_pos=right,
                 trim_primers=trim_primers,
                 reverse_complement_rev=reverse_complement_rev_primer,
                 chunk_size=chunk_size,
                 n_workers=n_workers,
                 client_address=client_address,
                 debug=debug,
                 )

    output_arts = (kmer_seqs, kmer_map)

    return output_arts

