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
    
    rep_fragments, = reconstruct_fragments(
        reconstruction_summary=reconstruction_summary,
        reconstruction_map=reconstruction_map,
        aligned_sequences=aligned_sequences,
        )

    tree, placements = insert_fragments(
        representative_sequences=rep_fragments,
        referencce_database=sepp_reference_database,
        n_threads=n_threads,
        )

    return (rep_fragments, tree, placements)
