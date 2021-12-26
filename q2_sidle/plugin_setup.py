import importlib

from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch, Visualization)
from q2_types.feature_data import (FeatureData,
                                   Sequence,
                                   Taxonomy,
                                   AlignedSequence,
                                   ) 
from q2_types.sample_data import (SampleData,
                                  )
from q2_types.feature_table import (FeatureTable, 
                                    Frequency, 
                                    )
from q2_types.tree import (Phylogeny, 
                           Rooted
                           )
from q2_fragment_insertion._type import (SeppReferenceDatabase, 
                                         Placements
                                         )
from q2_feature_table import heatmap_choices

from q2_sidle import (KmerMap, 
                      KmerMapFormat, 
                      KmerMapDirFmt, 
                      KmerAlignment, 
                      KmerAlignFormat, 
                      KmerAlignDirFmt,
                      SidleReconstruction, 
                      SidleReconFormat, 
                      SidleReconDirFormat,
                      ReconstructionSummary,
                      ReconSummaryFormat,
                      ReconSummaryDirFormat,  
                      AlignmentPosSummary,
                      AlignmentPosFormat, 
                      AlignmentPosDirFmt,                
                      )
import q2_sidle

citations = Citations.load('citations.bib', package='q2_sidle')

plugin = Plugin(
    name='sidle',
    version='2022.2-dev',
    website='https://github.com/jwdebelius/q2-sidle',
    package='q2_sidle',
    description=('This plugin reconstructs a full 16s sequence from short '
                 'reads over a marker gene region using the Short MUltiple '
                 'Read Framework (SMURF) algorithm.'),
    short_description='Plugin for kmer-based marker gene reconstruction.',
    citations=[citations['Debelius2021']],
)


# Functions
plugin.methods.register_function(function=q2_sidle.align_regional_kmers,
    name='Aligns ASV representative sequences to a regional kmer database.',
    description=('This takes an "amplified" region of the database and '
                 'performs alignment with representative ASV sequences. The '
                 'alignment assumes the ASVs and kmers start at the same '
                 'position in the sequence and that they are the same length.'
                 ),
    inputs={
        'kmers': FeatureData[Sequence],
        'rep_seq': FeatureData[Sequence],
    },
    outputs=[
        ('regional_alignment', FeatureData[KmerAlignment]),
        # ('discarded_sequences', FeatureData[Sequence]),

    ],
    parameters={
        'region': Str,
        'max_mismatch': Int % Range(0, None),
        'chunk_size':  (Int % Range(1, None)),
        'client_address': Str,
        'n_workers': Int % Range(1, None),
        'debug': Bool,
    },
    input_descriptions={
        'kmers': ('The reference kmer sequences from the database which have '
                  'had degenerate sequences expanded and duplicate sequences '
                  'identified.'),
        'rep_seq': ('The representative sequences for the ASVs being aligned. '
                    'These must be a consistent length.'),
    },
    output_descriptions={
        'regional_alignment': ('A mapping between the database kmer name and '
                               'the ASV.'),
    },
     parameter_descriptions={
        'region': ('A unique description of the hypervariable region being '
                   'aligned. Ideally, this matches the unique identifier '
                   'used during the regional extraction.'),
        'max_mismatch': ('The maximum number of mismatched nucleotides '
                         'allowed in mapping between a sequence and kmer.'),
        'chunk_size': ('The number of sequences to be analyzed in parallel '
                       'blocks. It is highly recommended this number stay '
                       'relatively small (<1000) in combination with parallel '
                       'processing (`n_workers`>1) for the best performance '
                       'and memory optimization.'),
        'n_workers': ('The number of jobs to initiate.'),
        'client_address': ('The IP address for an existing cluster. '
                          'Please see the dask client documentation for more '
                          'information: '
                          'https://distributed.dask.org/en/latest/client.html'
                          ),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options.'),
    },
    citations=[citations['Fuks2018']],
)


plugin.methods.register_function(function=q2_sidle.find_first_alignment_position,
    name='Finds the first position of a sequence in an alignment',
    description=('The function uses an alignment between regional ASV '
                 'representative sequences and a larger refernece alignment'
                 ' to map the representative sequences to a starting position'
                 ' on the representative sequence.'),
    inputs={
        'alignment': FeatureData[AlignedSequence],
        'representative_sequences': FeatureData[Sequence],
        'table': FeatureTable[Frequency],
    },
    outputs=[('position_summary', FeatureData[AlignmentPosSummary])],
    parameters={
        'direction': Str % Choices('fwd', 'rev'),
    },
    input_descriptions={
        'alignment': ('A muliple sequence alignment between the reference '
                      'database and the sequences represented in '
                      '`representative-sequences` to identify the starting '
                      'position of each sequence'),
        'representative_sequences': ("ASV sequences from mixed regions with "
                                     "unknown starting positions"),
        'table': ('The ASV table corresponding to `representative-sequences`.'
                  ' Some users may find this useful to set threshholds for '
                  'filtering or summarizing their data.'),
    },
    output_descriptions={
        'position_summary': ('The starting position, and if provided, total '
                              'counts, for each representaative sequence in'
                              ' the alignment.'),
    },
    parameter_descriptions={
        'direction': ('The direction of the read')
    },
)


plugin.methods.register_function(function=q2_sidle.prepare_extracted_region,
    name='Prepares an already extracted region to be a kmer database.',
    description=('This function takes an amplified region of the database, '
                 'expands the degenerate sequences and collapses the '
                 'duplicated sequences under a single id that can be '
                 'untangled later.'),
    inputs={
        'sequences': FeatureData[Sequence]
    },
    outputs=[
        ('collapsed_kmers', FeatureData[Sequence]), 
        ('kmer_map', FeatureData[KmerMap]),
    ],
    parameters={
        'trim_length': Int,
        'region': Str,
        'fwd_primer': Str,
        'rev_primer': Str,
        'reverse_complement_rev': Bool,
        'reverse_complement_result': Bool,
        'chunk_size':  (Int % Range(1, None)),
        'n_workers': Int % Range(1, None),
        'client_address': Str,
        'debug': Bool,

    },
    input_descriptions={
        'sequences': 'The full length sequences from the reference database.',
    },
    output_descriptions={
        'collapsed_kmers': ('Reference kmer sequences for the region with '
                            'the degenerate sequences expanded and '
                            'duplicated sequences identified.'
                            ),
        'kmer_map': ('A mapping relationship between the name of the '
                     'sequence in the database and the kmer identifier used '
                     'in this region.'),
    },
    parameter_descriptions={
        'region': ('A unique description of the hypervariable region being '
                   'extracted.'),
        'trim_length': ('The length of the extracted regional kmers.'),
        'fwd_primer': ('The forward primer used to amplify the region of '
                       'interest.'),
        'rev_primer': ('The reverse primer used to amplify the region of '
                       'interest.'),
        'reverse_complement_rev': ('If the reverse primer was reverse '
                                   'complemented during sequence extraction. '
                                   'This is used to later generate fragments '
                                   'for the phylogenetic tree.'),
        'reverse_complement_result': ('Whether the sequences for alignment '
                                      'should be reverse complemented for '
                                      'alignment, for example, in cases where '
                                      'the forward and reverse primers do '
                                      'not overlap and you want to align with '
                                      'the reverse sequence.'),
        'chunk_size': ('The number of sequences to be analyzed in parallel '
                       'blocks.'),
        'n_workers': ('The number of jobs to initiate.'),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options.'),
        'client_address': ('The IP address for an existing cluster. '
                           'Please see the dask client documentation for more '
                           'information: '
                           'https://distributed.dask.org/en/latest/client.html'
                           ),
    },
    citations=[citations['Fuks2018']],
)


plugin.methods.register_function(function=q2_sidle.reconstruct_counts,
    name='Reconstructs multiple aligned regions into a count table.',
    description=('This takes multiple regional alignments and regional counts '
                 'and reconstructs them into a single table with region-'
                 'normalized abundance counts.'),
    inputs={
        'regional_alignment': List[FeatureData[KmerAlignment]],
        'regional_table': List[FeatureTable[Frequency]],
        'database_map': FeatureData[SidleReconstruction],
        'database_summary': FeatureData[ReconstructionSummary],
    },
    outputs=[
        ('reconstructed_table', FeatureTable[Frequency]),
    ],
    parameters={
        'region': List[Str],
        'per_nucleotide_error': Float % Range(0, 1),
        'min_abund': Float % Range(0, 1),
        'region_normalize': Str % Choices('average', 'weighted', 'unweighted'),
        'min_counts': Int % Range(0, None),
        'block_size': Int,
        'n_workers': Int % Range(1, None),
        'client_address': Str,
        'debug': Bool,
    },
    input_descriptions={
        'regional_alignment': ('A mapping between the kmer names (in the kmer'
                               ' map) and the features (found in the regional'
                               ' table)'),
        'regional_table': ('A feature-table for each region, where  the '
                           'features in the table correspond to the ASVs '
                           'which were aligned in the regional alignment '
                           'artifact'),
        'database_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
         'database_summary': ('A summary of the statitics for the '
                              'regional map describing the number of '
                              'regions mapped to each reference sequence'
                              ' and the number of kmers. The kmer '
                              'mapping estimate can account for '
                              'degeneracy when the `--count-degenerates`'
                              ' flag is used or can ignore degenrate '
                              'sequences in mapping'),
    },
    output_descriptions={
        'reconstructed_table': ('The feature table with the reconstructed '
                                'abundance using ASVs from all regions mapped'
                                ' to database sequences.'),
    },
    parameter_descriptions={
        'region': ('The name of the sub region used in alignment. The region'
                   ' names do not matter, however, the region order must '
                   'match the order along the hypervariable region.'),
        'per_nucleotide_error': ('The assumed per-nucelotide error rate '
                                 'throughout amplification and sequencing.'),
        'min_abund': ('The minimum frequency for a feature to be retained '
                      'during reconstruction. The default from the original '
                      'smurf algorithm is 1e-10, the number here may depend '
                      'on the total number of sequences in your sample.'),
        'min_counts': ('The mininum depth across all regions after alignment '
                       'for a sample to be included in a reconstruction.'),
        'region_normalize': ('Whether the relative abundance should be '
                             'normalized by region during reconstruction. '
                             'When using kmer-based alignment to '
                             'reconstruct multiple regions within a single '
                             'sample, this is most appropriate. In meta-'
                             'analysis, region normalization may increase '
                             'count retention.'),
        'block_size': ('The number of sequences to use in parallel '
                       'computation. The larger the block_size, the faster '
                       'processing can happen, but more the memory that will '
                       'be required.'),
        'n_workers': ('The number of jobs to initiate.'),
        'client_address': ('The IP address for an existing cluster. '
                          'Please see the dask client documentation for more '
                          'information: '
                          'https://distributed.dask.org/en/latest/client.html'
                          ),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options.'),
    },
    citations=[citations['Fuks2018']],
)


plugin.methods.register_function(function=q2_sidle.reconstruct_database,
    name='Reconstructs regional kmers into full database names',
    description=('Uses the kmer maps and regional alignments to solve the '
                 'regional sidle database'),
    inputs={'regional_alignment': List[FeatureData[KmerAlignment]],
            'kmer_map': List[FeatureData[KmerMap]],
            },
    outputs=[('database_map', FeatureData[SidleReconstruction]),
             ('database_summary', FeatureData[ReconstructionSummary]),
             ],
    parameters={'region': List[Str],
                'count_degenerates': Bool,
                'block_size': Int,
                'n_workers': Int % Range(1, None),
                'client_address': Str,
                 'debug': Bool,
                 },
    input_descriptions={
        'regional_alignment': ('A mapping between the kmer names (in the kmer'
                               ' map) and the features (found in the regional'
                               ' table)'),
        'kmer_map': ('A mapping relationship between the name of the '
                     'sequence in the database and the kmer identifier used'
                     ' in this region. The kmer map should correspond to the '
                     'kmers used in regional alignment'),
        },
    output_descriptions={
        'database_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
        'database_summary': ('A summary of the statitics for the '
                                   'regional map describing the number of '
                                   'regions mapped to each reference sequence'
                                   ' and the number of kmers. The kmer '
                                   'mapping estimate can account for '
                                   'degeneracy when the `--count-degenerates`'
                                   ' flag is used or can ignore degenrate '
                                   'sequences in mapping'),
        },
    parameter_descriptions={
        'region': ('The name of the sub region used in alignment. The region'
                   ' names do not matter, however, the region order must '
                   'match the order along the hypervariable region.'),
        'count_degenerates': ("Whether sequences which contain degenerate "
                              "nucleotides should be counted as unqiue kmers"
                              " or whether the number of original database "
                              "sequences before degenerate expansion should "
                              "be used."),
        'block_size': ('The number of sequences to use in parallel '
                       'computation. The larger the block_size, the faster '
                       'processing can happen, but the more memory that will'
                       ' be required.'),
        'n_workers': ('The number of jobs to initiate.'),
        'client_address': ('The IP address for an existing cluster. '
                          'Please see the dask client documentation for more'
                          ' information: '
                          'https://distributed.dask.org/en/latest/client.html'
                          ),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
        },
)


plugin.methods.register_function(function=q2_sidle.reconstruct_fragment_rep_seqs,
    name='Reconstract representative sequences for shared fragments.',
    description=('EXPERIMENTAL!!!\n'
                 'This function simulates a represenative sequence for '
                 'reference regions that are derived from multiple sequences '
                 'to allow tree building via fragment insertion. The function '
                 'will find the consensus sequence for all the database '
                 'regions covered between the amplicons.'
                 ),
    inputs={
        'reconstruction_summary': FeatureData[ReconstructionSummary],
        'reconstruction_map': FeatureData[SidleReconstruction],
        'aligned_sequences': FeatureData[AlignedSequence]
    },
    outputs=[
        ('representative_fragments', FeatureData[Sequence]),
    ],
    parameters={},
    input_descriptions={
        'reconstruction_summary': ('A summary of the statitics for the '
                                   'regional map describing the number of '
                                   'regions mapped to each reference sequence '
                                   'and the number of kmers. The kmer '
                                   'mapping estimate can account for '
                                   'degeneracy when the `--count-degenerates` '
                                   'flag is used or can ignore degenrate '
                                   'sequences in mapping.'),
        'reconstruction_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
        'aligned_sequences': ('The aligned representative sequences '
                              'corresponding to the database used in '
                              'reconstruction.'),
    },
    output_descriptions={
        'representative_fragments': ('The consensus sequence fragments '
                                     'to be used for fragment insertion.')
    },
    parameter_descriptions={},
)


plugin.methods.register_function(function=q2_sidle.reconstruct_taxonomy,
    name='Reconstructs taxonomic strings for a reconstructed sidle table.',
    description=('Reconstructs the taxonomic annotation based on a sidle '
                 'database by identifying the lowest taxonomic level  where'
                 ' the taxonomic annotation diverges for two sequences.'),
    inputs={
        'reconstruction_map': FeatureData[SidleReconstruction],
        'taxonomy': FeatureData[Taxonomy]
    },
    outputs=[
        ('reconstructed_taxonomy', FeatureData[Taxonomy])
    ],
    parameters={
        'database': Str % Choices('none', 'greengenes', 'silva'),
        'define_missing': Str % Choices('merge', 'inherit', 'ignore'),
        'ambiguity_handling': Str % Choices('missing', 'ignore'),
    },
    input_descriptions={
        'reconstruction_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
        'taxonomy': ('The taxonomic strings from the database.'),
    },
    output_descriptions={
        'reconstructed_taxonomy': ('Taxonomy which addresses the ambiguity '
                                   'in assignment associated with mapping '
                                   'to multiple sequences in reconstruction.')
    },
    parameter_descriptions={
        'database': ('The taxonomic database being used. Currently, the only '
                    'two supported databases are Greengenes and Silva. '
                    'The database choice influences the handling of missing '
                    'and ambiguous taxa.'),
        'define_missing':  ('Taxonomic strings may be missing information '
                             '(for example  `g__` in greengenes  or '
                             '`D_5__uncultured bacteria` in Silva). These '
                             'can be ignored (`"ignore"`) and treated like '
                             'any other taxonomic designation; they can be '
                             'first inherited in merged sequences '
                             '(`"merge"`), where, when there are two strings '
                             ' being merged and one has a missing level, '
                             'the missing level is taken form the defined '
                             'one, or they can be inherited from the '
                             'previous level (`"inherit"`) first, and then '
                             'merged.'),
        'ambiguity_handling': ('whether "ambiguous taxa" (Silva-specific) '
                               'should be treated as missing values '
                               '(`"missing"`) or ignored (`"ignore"`).')
    },
    citations=[citations['Fuks2018']],
)


seq_match = TypeMatch([Sequence, AlignedSequence])
plugin.methods.register_function(function=q2_sidle.reverse_complement_sequence,
    name='Reverse Complements a sequence',
    description=('This function reverse complements a sequence'),
    inputs={'sequence': FeatureData[seq_match]},
    outputs=[('reverse_complement', FeatureData[seq_match])],
    parameters={},
    input_descriptions={
        'sequence': ('The sequences to be reverse complemented'),
        },
    output_descriptions={
        'reverse_complement': 'The reverse complement of the input sequences',
        },
    parameter_descriptions={},
)


plugin.methods.register_function(function=q2_sidle.trim_dada2_posthoc,
    name='Trim a dada2 ASV table and rep set to a consistent length.',
    description=('This function trims ASVs generated by DADA2 to a '
                 'consistent length and collapses them into a single '
                 'sequence as appropriate. This is necessary to be able to '
                 'use a DADA2 ASV table with the SMURF algorithm, which '
                 'requires a consistent sequence length.\nNOTE: This should '
                 'be done before anything else that requires ASV sequences '
                 'because it may change ASV IDs.'
                 ),
    inputs={
        'table': FeatureTable[Frequency],
        'representative_sequences': FeatureData[Sequence]
    },
    outputs=[
        ('trimmed_table', FeatureTable[Frequency]),
        ('trimmed_representative_sequences', FeatureData[Sequence]),
    ],
    parameters={
        'trim_length': Int % Range(0, None),
        'hashed_feature_ids': Bool,
    },
    input_descriptions={
        'table': ('A feature table generated by a denoising algorithm with '
                  'variable length ASVs.'),
        'representative_sequences': ('The corresponding representative '
                                     'sequences for the ASV table.'),
    },
    output_descriptions={
        'trimmed_table': ('ASV table with consistent length ASVs.'),
        'trimmed_representative_sequences': ('ASV sequences of a consistent '
                                             'length.')
    },
    parameter_descriptions={
        'trim_length': ('The length ASVs should be trimmed to. If the trim '
                        'length is 0, the shortest sequence will be used. '
                        'Alternatively, if a trim length is specified which '
                        'is longer than some of the sequences, those ASVs '
                        'will be discarded.'),
        'hashed_feature_ids': ('If true, the feature ids in the resulting '
                               'table will be presented as hashes of the '
                               'sequences defining each feature. The hash '
                               'will always be the same for the same sequence '
                               'so this allows feature tables to be merged '
                               'across runs of this method. You should only '
                               'merge tables if the exact same parameters '
                               'are used for each run.')
    }
)

# Visualizations
plugin.visualizers.register_function(function=q2_sidle.summarize_alignment_positions,
    name='Generate a heatmap summarizing the starting position of sequences',
    description=('Summarizes the number of sequences at a starting count and '
                 'visualizes the aligned representative sequences'),
    inputs={
        'alignment': FeatureData[AlignedSequence],
        'position_summary': FeatureData[AlignmentPosSummary],
    },
    parameters={
        'sort_cols': Str,
        'weight_by_abundance': Bool,
        'colormap': Str % Choices(heatmap_choices['color_scheme']),
        'heatmap_maskcolor': Str,
        'heatmap_grid': Bool,
        'tick_interval': Int,
    },
    input_descriptions={
        'alignment': ('A muliple sequence alignment between the reference '
                      'database and represetnative sequences of interest'),
        'position_summary': ('The starting position, and if provided, total '
                             'counts, for each representaative sequence in'
                             ' the alignment.'),
    },
    parameter_descriptions={
        'sort_cols': ('The column in the position summary to be used to sort'
                      ' the aligned sequence'),
        'weight_by_abundance': ('If the abundance information is present in '
                                'the position summary (if, for example, a '
                                'table was based when the starting position '
                                'was identified), this will color the heatmap'
                                ' by the total sequence abundance'),
        'colormap': 'The matplotlib colorscheme to generate the heatmap with',
        'heatmap_maskcolor': ('The color to use as a background to highlight'
                              ' specifically where sequences are present. '
                              'Otherwwise, they will be displayed based on '
                              'the colormap specified'),
        'heatmap_grid': 'When true, a grid will be displayed on the heatmap',
        'tick_interval': ('How frequently xticks, indicating position in the'
                          ' alignment, should be displayed on the heatmap'),
    },
    citations=[citations['Hunter2007Matplotlib']]
)

# Pipelines
plugin.pipelines.register_function(function=q2_sidle.map_alignment_positions,
    name='Finds the starting positions of denoised amplicons in an alignment',
    description=('For studies without primers, this will take amplicons, '
                 'align them against a provided reference to identify '
                 'starting positions, and provide a map to split the '
                 'data into regions.'),
    inputs={
        'alignment': FeatureData[AlignedSequence],
        'sequences': FeatureData[Sequence],
        'table': FeatureTable[Frequency],
        },
    outputs=[
        ('expanded_alignment', FeatureData[AlignedSequence]),
        ('position_summary', FeatureData[AlignmentPosSummary]),
        ('position_map', Visualization),
        ],
    parameters={
        'direction': Str % Choices('fwd', 'rev'),
        'reverse_complement_ref': Bool,
        'n_threads': Int % Range(1, None),
        'add_fragments': Bool,
        'colormap': Str % Choices(heatmap_choices['color_scheme']),
        },
    input_descriptions={
        'alignment': ('The reference multiple sequence alignment which '
                      'representatative ASVs should be aligned against'),
        'sequences': ('The representative ASV sequences that need to be '
                      'positioned to find their starting regions'),
        'table': ('The feature table corresponding to the sequences. This is '
                  'optional and can be used to identify starting positions'
                  ' in the alignment with more abundance reads.')
    },
    output_descriptions={
        'expanded_alignment': ("The multiple sequence alignment including the"
                               " representative sequences"),
        'position_summary': ("Describes the relationship between an ASV, "
                             "and its starting position in the multiple "
                             "sequence alignment. The per-ASV information can"
                             " be viewed by tabulating the metadata"),
        'position_map': ("A summary of the number of sequences mapped to "
                         "different starting positions along the multiple "
                         "sequence alignment"),
    },
    parameter_descriptions={
        'direction': ("The orientation of the sequences being aligned"),
        'reverse_complement_ref': ('If the sequences are reverse complemented'
                                   ', this will also reverse complement the'
                                   ' reference sequences within the pipeline'
                                   ),
        'n_threads': ('The number of threads to use in multiple sequence '
                      'alignment.'),
        'add_fragments': ('Whether the ASV sequences should be added to the'
                          ' alignment as full sequences, or if they should '
                          'be allowed to be inserted with gaps. This is '
                          'most useful for short sequences like primers.'),
        'colormap': ('The colormap used to render the heatmap of aligned '
                     'sequences'),
        },
    )



plugin.pipelines.register_function(function=q2_sidle.sidle_reconstruction,
    name="A pipeline to reconstruct the database, count table, and taxonomy",
    description=("A pipeline to reconstruct a database, count table, and "
                 "taxonomy"),
    inputs={'kmer_map': List[FeatureData[KmerMap]],
            'regional_alignment': List[FeatureData[KmerAlignment]],
            'regional_table': List[FeatureTable[Frequency]],
            'reference_taxonomy': FeatureData[Taxonomy],
            },
    outputs=[('database_map', FeatureData[SidleReconstruction]),
             ('database_summary', FeatureData[ReconstructionSummary]),
             ('reconstructed_table', FeatureTable[Frequency]),
             ('reconstructed_taxonomy', FeatureData[Taxonomy]),
             ],
    parameters={'region': List[Str],
                'min_counts': Int % Range(0, None),
                'database': Str % Choices('none', 'greengenes', 'silva'),
                'define_missing': Str % Choices('merge', 'inherit', 'ignore'),
                'block_size': Int,
                'n_workers': Int % Range(1, None),
                'client_address': Str,
                'debug': Bool,
                },
    input_descriptions={
        'kmer_map': ('A mapping relationship between the name of the '
                     'sequence in the database and the kmer identifier used'
                     ' in this region. The kmer map should correspond to the '
                     'kmers used in regional alignment'),
        'regional_alignment': ('A mapping between the kmer names (in the kmer'
                               ' map) and the features (found in the regional'
                               ' table)'),
        'regional_table': ('A feature-table for each region, where  the '
                           'features in the table correspond to the ASVs '
                           'which were aligned in the regional alignment '
                           'artifact'),
        'reference_taxonomy': ('The taxonomic strings from the database'),
        },
    output_descriptions={
        'database_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
        'database_summary': ('A summary of the statitics for the '
                                   'regional map describing the number of '
                                   'regions mapped to each reference sequence'
                                   ' and the number of kmers. The kmer '
                                   'mapping estimate can account for '
                                   'degeneracy when the `--count-degenerates`'
                                   ' flag is used or can ignore degenrate '
                                   'sequences in mapping'),
        'reconstructed_table': ('The feature table with the reconstructed '
                                'abundance using ASVs from all regions mapped'
                                ' to database sequences.'),
        'reconstructed_taxonomy': ('Taxonomy which addresses the ambiguity'
                                   ' in assignment associated with mapping '
                                   'to multiple sequences in reconstruction.')
        },
    parameter_descriptions={
        'region': ('The name of the sub region used in alignment. The region'
                   ' names do not matter, however, the region order must '
                   'match the order along the hypervariable region.'),
        'min_counts': ('The mininum depth across all regions after alignment '
                       'for a sample to be included in a reconstruction'),
        'database': ('The taxonomic database being used. Currently, the only'
                    ' two supported databases are Greengenes and Silva. '
                    'The database choice influences the handling of missing '
                    'and ambigious taxa.'),
        'define_missing':  ('Taxonomic strings may be missing information '
                             '(for example  `g__` in greengenes  or '
                             '`D_5__uncultured bacteria` in Silva). These '
                             'can be ignored  (`"ignore"`) and treated like'
                             ' any other taxonomic designation; they can be'
                             ' first inherited in merged sequences '
                             '(`"merge`"), where, when there are two strings'
                             ' being merged and one has a missing level, '
                             'the missing level is taken form the defined'
                             ' one, or they can be inherited from the '
                             'previous level (`"inherit"`) first, and then '
                             'merged.'),
        'block_size': ('The number of sequences to use in parallel '
                       'computation. The larger the block_size, the faster '
                       'processing can happen, but the more memory that will'
                       ' be required.'),
        'n_workers': ('The number of jobs to initiate.'),
        'client_address': ('The IP address for an existing cluster. '
                          'Please see the dask client documentation for more'
                          ' information: '
                          'https://distributed.dask.org/en/latest/client.html'
                          ),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
        },
)


plugin.pipelines.register_function(function=q2_sidle.reconstruct_tree,
    name=("A pipeline to build a phylogenetic tree based on reconstructed "
          "sequences"),
    description=("A pipeline to reconstruct a database, count table, and "
                 "taxonomy\nThis pipeline is somewhat experimental.\n"
                 "Representative sequences for non-unique sequences are "
                 "assigned based on the concensus sequence for the "
                 "amplicons"),
    inputs={
        'reconstruction_summary': FeatureData[ReconstructionSummary],
        'reconstruction_map': FeatureData[SidleReconstruction],
        'aligned_sequences': FeatureData[AlignedSequence],
        'sepp_reference_database': SeppReferenceDatabase,
    },
    outputs=[('representative_fragments', FeatureData[Sequence]),
             ('tree', Phylogeny[Rooted]),
             ('placements', Placements)
             ],
    parameters={'n_threads': Int % Range(1, None)},
    input_descriptions={
        'reconstruction_summary': ('A summary of the statitics for the '
                                   'regional map describing the number of '
                                   'regions mapped to each reference sequence '
                                   'and the number of kmers. The kmer '
                                   'mapping estimate can account for '
                                   'degeneracy when the `--count-degenerates` '
                                   'flag is used or can ignore degenrate '
                                   'sequences in mapping.'),
        'reconstruction_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
        'aligned_sequences': ('The aligned representative sequences '
                              'corresponding to the database used in '
                              'reconstruction.'),
        'sepp_reference_database': ('A SEPP reference database object '
                                    'corresponding to the database you used '
                                    'for reconstruction.\n'
                                    'IMPORTANT: The database versions must '
                                    'match for the table construction and '
                                    'tree.'),
        },
    output_descriptions={
        'representative_fragments': ('The consensus sequence fragments '
                                     'to be used for fragment insertion.'),
        'tree': 'The tree with all the reconstructed features',
        'placements': ('The location of the representative fragments '
                       'in the phylogenetic tree')
        },
    parameter_descriptions={
        'n_threads': 'the number of threads to use during fragment insertion',
        }, 
)



# Registers semantic types
plugin.register_formats(KmerMapFormat, 
                        KmerMapDirFmt, 
                        KmerAlignFormat, 
                        KmerAlignDirFmt,
                        SidleReconFormat, 
                        SidleReconDirFormat,
                        ReconSummaryFormat,
                        ReconSummaryDirFormat,
                        AlignmentPosFormat,
                        AlignmentPosDirFmt,
                        )


plugin.register_semantic_types(KmerMap, 
                               KmerAlignment,
                               SidleReconstruction,
                               ReconstructionSummary,
                               AlignmentPosSummary,
                               )

plugin.register_semantic_type_to_format(FeatureData[KmerMap], 
                                        KmerMapDirFmt)


plugin.register_semantic_type_to_format(FeatureData[KmerAlignment], 
                                        KmerAlignDirFmt)


plugin.register_semantic_type_to_format(FeatureData[SidleReconstruction], 
                                        SidleReconDirFormat)


plugin.register_semantic_type_to_format(FeatureData[ReconstructionSummary], 
                                        ReconSummaryDirFormat)


plugin.register_semantic_type_to_format(FeatureData[AlignmentPosSummary],
                                        AlignmentPosDirFmt)


importlib.import_module('q2_sidle._transformer')
