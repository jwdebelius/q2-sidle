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




    if direction == 'rev':
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




