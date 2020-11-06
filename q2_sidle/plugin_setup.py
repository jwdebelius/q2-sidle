import importlib

from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch)
from q2_types.feature_data import (FeatureData,
                                   Sequence,
                                   Taxonomy,
                                   AlignedSequence,
                                   ) 
from q2_types.feature_table import (
                      FeatureTable, 
                      Frequency, 
                      )
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
                      )
import q2_sidle

citations = Citations.load('citations.bib', package='q2_sidle')

plugin = Plugin(
    name='sidle',
    version='2020.08',
    website='https://github.com/jwdebelius/q2-sidle',
    package='q2_sidle',
    description=('This plugin reconstructs a full 16s sequence from short '
                 'reads over a marker gene region using the Short MUltiple '
                 'Read Framework (SMURF) algorithm'),
    short_description='Plugin for kmer-based marker gene reconstruction.',
)

plugin.methods.register_function(
    function=q2_sidle.filter_degenerate_sequences,
    name='Filter degenerate sequences',
    description=('Prefiltering a sequence database to remove sequences with '
                 'too many degenerate nucleotides'),
    inputs={
        'sequences': FeatureData[Sequence]
    },
    outputs=[
        ('filtered_sequences', FeatureData[Sequence]),
    ],
    parameters={
        'max_degen': (Int % Range(0, None)),
        'chunk_size':  (Int % Range(1, None)),
        'n_workers': Int % Range(0, None),
        'debug': Bool,
        'client_address': Str,
    },
    input_descriptions={
        'sequences': 'The sequences to be filtered.'
    },
    output_descriptions={
        'filtered_sequences': ('The sequences filtered to remove full length'
                               ' sequences with excess degenerate '
                               'nucleotides')
    },
    parameter_descriptions={
        'max_degen': ('the maximum number of degenerate sequences for a '
                      'sequence to be retained.'),
        'chunk_size': ('The number of sequences to be analyzed in parallel '
                       'blocks'),
        'n_workers': ('The number of jobs to initiate. When `n_workers` is 0,'
                      ' the cluster will be able to access all avaliable'
                      ' resources.'),
        'client_address': ('The IP address for an existing cluster. '
                          'Please see the dask client documentation for more'
                          ' information: '
                          'https://distributed.dask.org/en/latest/client.html'
                          ),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
    },
)


plugin.methods.register_function(
    function=q2_sidle.prepare_extracted_region,
    name='Prepares an already extracted region to be a kmer database',
    description=('This function takes an amplified region of the database, '
                 'expands the degenerate sequences and collapses the '
                 'duplciated sequences under a single id that can be '
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
        'n_workers': Int % Range(0, None),
        'client_address': Str,
        'debug': Bool,

    },
    input_descriptions={
        'sequences': 'The full length sequences from the reference database',
    },
    output_descriptions={
        'collapsed_kmers': ('Reference kmer sequences for the region with'
                            ' the degenerate sequences expanded and '
                            'duplicated sequences identified'
                            ),
        'kmer_map': ('A mapping relationship between the name of the '
                     'sequence in the database and the kmer identifier used'
                     ' in this region.'),
    },
    parameter_descriptions={
        'region': ('A unique description of the hypervariable region being '
                   'extracted.'),
        'trim_length': ('The length of the extracted regional kmers.'),
        'fwd_primer': ('The forward primer used to amplify the region of '
                       'interest'),
        'rev_primer': ('The reverse primer used to amplify the region of '
                       "interest"),
        'reverse_complement_rev': ('If the reverse primer was reverse '
                                   'complemented during sequence extraction. '
                                   'This is used to later generate fragments '
                                   'for the phylogenetic tree.'),
        'reverse_complement_result': ('Whether the sequences for alignment '
                                      'should be reverse complemented for '
                                      'alignment, for example, in cases where'
                                      ' the forward and reverse primers do '
                                      'not overlap and you want to align with'
                                      ' the reverse sequence.'),
        'chunk_size': ('The number of sequences to be analyzed in parallel '
                       'blocks'),
        'n_workers': ('The number of jobs to initiate. When `n_workers` is 0,'
                      ' the cluster will be able to access all avaliable'
                      ' resources.'),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
        'client_address': ('The IP address for an existing cluster. '
                           'Please see the dask client documentation for more'
                           ' information: '
                           'https://distributed.dask.org/en/latest/client.html'
                           ),
    },
    citations=[citations['Fuks2018']],

)

plugin.methods.register_function(
    function=q2_sidle.align_regional_kmers,
    name='Aligns ASV representative sequences to a regional kmer database',
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
        'n_workers': Int % Range(0, None),
        'debug': Bool,
    },
    input_descriptions={
        'kmers': ('The reference kmer sequences from the database which have'
                  ' had degenerate sequences expanded and duplicate sequences'
                  ' identified'),
        'rep_seq': ('The representative sequences for the ASVs being aligned.'
                    'These must be a consistent length.'),
    },
    output_descriptions={
        'regional_alignment': ('A mapping between the database kmer name and'
                               ' the asv'),
    },
     parameter_descriptions={
        'region': ('A unique description of the hypervariable region being '
                   'aligned. Ideally, this matches the unique identifier '
                   'used during the regional extraction.'),
        'max_mismatch': ('the maximum number of mismatched nucleotides '
                         'allowed in mapping between a sequence and kmer'),
        'chunk_size': ('The number of sequences to be analyzed in parallel '
                       'blocks. It is highly recommend this number stay '
                       'relatively small (>1000) in combination with parallel'
                       ' processing (`n_workers`>1) for the best performance'
                       ' and memory optimization.'),
        'n_workers': ('The number of jobs to initiate. When `n_workers` is 0,'
                      ' the cluster will be able to access all avaliable'
                      ' resources.'),
        'client_address': ('The IP address for an existing cluster. '
                          'Please see the dask client documentation for more'
                          ' information: '
                          'https://distributed.dask.org/en/latest/client.html'
                          ),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
    },
    citations=[citations['Fuks2018']],
)


plugin.methods.register_function(
    function=q2_sidle.reconstruct_counts,
    name='Reconstructs multiple aligned regions into a count table',
    description=('This takes multiple regional alignments and regional counts'
                 ' and reconstructs them into a single table with region-'
                 'normalized abundance counts.'),
    inputs={
        'regional_alignment': List[FeatureData[KmerAlignment]],
        'kmer_map': List[FeatureData[KmerMap]],
        'regional_table': List[FeatureTable[Frequency]],
    },
    outputs=[
        ('reconstructed_table', FeatureTable[Frequency]),
        ('reconstruction_summary', FeatureData[ReconstructionSummary]),
        ('reconstruction_map', FeatureData[SidleReconstruction])
    ],
    parameters={
        'region': List[Str],
        'per_nucleotide_error': Float % Range(0, 1),
        'min_abund': Float % Range(0, 1),
        'count_degenerates': Bool,
        'region_normalize': Str % Choices('average', 'weighted', 'unweighted'),
        'min_counts': Int % Range(0, None),
        'block_size': Int,
        'n_workers': Int % Range(0, None),
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
                     'kmers used in regional alignmeent'),
        'regional_table': ('A feature-table for each region, where  the '
                           'features in the table correspond to the ASVs '
                           'which were aligned in the regional alignment '
                           'artifact')
    },
    output_descriptions={
        'reconstructed_table': ('The feature table with the reconstructed '
                                'abundance using ASVs from all regions mapped'
                                ' to database sequences.'),
        'reconstruction_summary': ('A summary of the statitics for the '
                                   'regional map describing the number of '
                                   'regions mapped to each reference sequence'
                                   ' and the number of kmers. The kmer '
                                   'mapping estimate can account for '
                                   'degeneracy when the `--count-degenerates`'
                                   ' flag is used or can ignore degenrate '
                                   'sequences in mapping'),
        'reconstruction_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
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
        'per_nucleotide_error': ('The assumed per-nucelotide error rate '
                                 'through out amplification and sequencing'),
        'min_abund': ('The minimum frequency for a feature to be retained '
                      'during reconstruction. The default from the original '
                      'smurf algorithm is 1e-10, the number here may depend'
                      ' on the total number of sequences in your sample'),
        'min_counts': ('The mininum depth across all regions after alignment '
                       'for a sample to be included in a reconstruction'),
        'region_normalize': ('Whether the relative abundance should be '
                             'normalized by region during reconstruction. '
                             'When using kmer-based alignment to '
                             'reeconstruct multiple regions within a single '
                             'sample, this is most appropriate. In meta-'
                             'analysis, region normalization may increase'
                             ' count retention.'),
        'block_size': ('The number of sequences to use in parallel '
                       'computation. The larger the block_size, the faster '
                       'processing can happen, but the more memory that will'
                       ' be required.'),
        'n_workers': ('The number of jobs to initiate. When `n_workers` is 0,'
                      ' the cluster will be able to access all avaliable'
                      ' resources.'),
        'client_address': ('The IP address for an existing cluster. '
                          'Please see the dask client documentation for more'
                          ' information: '
                          'https://distributed.dask.org/en/latest/client.html'
                          ),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
    },
    citations=[citations['Fuks2018']],
)


plugin.methods.register_function(
    function=q2_sidle.reconstruct_taxonomy,
    name='Reconstructs taxonomic strings for a reconstructed sidle table',
    description=('Reconstructs the taxonomic annotation based on a sidle '
                 'database by identifying the lowest taxonomic level  where'
                 ' the taxonomic annotation diverges for two sequences'),
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
        'taxonomy': ('The taxonomic strings from the database'),
    },
    output_descriptions={
        'reconstructed_taxonomy': ('Taxonomy which addresses the ambiguity'
                                   ' in assignment associated with mapping '
                                   'to multiple sequences in reconstruction.')
    },
    parameter_descriptions={
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
        'ambiguity_handling': ('whether "ambigious taxa" (Silva-specific) '
                               'should be treated as missing values '
                               '(`"missing`") or ignored (`"ignore"`)')
    },
    citations=[citations['Fuks2018']],
)


plugin.methods.register_function(
    function=q2_sidle.trim_dada2_posthoc,
    name='Trim a dada2 ASV table and rep set to a consistent length',
    description=('This function trims ASVs generated by DADA2 to a '
                 'consistent length and collapses them into a single '
                 'sequence as appropriate. This is necessary to be able to '
                 'use a DADA2 ASV table with the SMURF algorithm, which '
                 'requires a consistent sequence length.\nNOTE: This shoud'
                 ' be done before anything else that requires ASV sequences'
                 ' because it may change ASV IDs'
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
                  'variable length ASVs'),
        'representative_sequences': ('The corresponding representative '
                                     'sequences for the ASV table'),
    },
    output_descriptions={
        'trimmed_table': ('ASV table with consistent length ASVs'),
        'trimmed_representative_sequences': ('ASV sequences of a consistent '
                                             'length')
    },
    parameter_descriptions={
        'trim_length': ('The length ASVs should be trimmed to. If the trim '
                        'length is 0, the shortest sequence will be used.'
                        'Alternatively, if a trim length is specified which'
                        ' is longer than some of the sequences, those ASVs'
                        ' will be discarded.'),
        'hashed_feature_ids': ('If true, the feature ids in the resulting '
                               'table will be presented as hashes of the '
                               'sequences defining each feature. The hash '
                               'will always be the same for the same sequence'
                               ' so this allows feature tables to be merged '
                               'across runs of this method. You should only'
                               ' merge tables if the exact same parameters'
                               ' are used for each run.')
    }
)


plugin.methods.register_function(
    function=q2_sidle.reconstruct_fragment_rep_seqs,
    name='Reconstract representative sequences for shared fragments',
    description=('EXPERIMENTAL!!!\n'
                 'This function simulates a represenative sequence for '
                 'reference regions that are dervied from multiple sequences'
                 'to allow tree building via fragment insertion. The function'
                 'will find the concensus sequence for all the database '
                 'regions covered between the amplicons'
                 ),
    inputs={
        'kmer_map': List[FeatureData[KmerMap]],
        'reconstruction_summary': FeatureData[ReconstructionSummary],
        'reconstruction_map': FeatureData[SidleReconstruction],
        'aligned_sequences': FeatureData[AlignedSequence]
    },
    outputs=[
        ('representative_fragments', FeatureData[Sequence]),
    ],
    parameters={
        'region': List[Str],
    },
    input_descriptions={
        'kmer_map': ('A mapping relationship between the name of the '
                     'sequence in the database and the kmer identifier used'
                     ' in this region. The kmer map should correspond to the '
                     'kmers used in regional alignment'),
        'reconstruction_summary': ('A summary of the statitics for the '
                                   'regional map describing the number of '
                                   'regions mapped to each reference sequence'
                                   ' and the number of kmers. The kmer '
                                   'mapping estimate can account for '
                                   'degeneracy when the `--count-degenerates`'
                                   ' flag is used or can ignore degenrate '
                                   'sequences in mapping'),
        'reconstruction_map': ('A map between the final kmer name and the '
                               'original database sequence. Useful for '
                               'reconstructing taxonomy and trees.'),
        'aligned_sequences': ('The aligned representative sequences '
                              'corresponding to the database used in '
                              'reconstruction.'),
    },
    output_descriptions={
        'representative_fragments': ('The concensus sequence fragments '
                                     'to be used for fragment insertion')
    },
    parameter_descriptions={
        'region': ('The name of the sub region used in alignment. The region'
                   ' names do not matter, however, the region order must '
                   'match the order along the hypervariable region.'),

    }
)


plugin.register_formats(KmerMapFormat, 
                        KmerMapDirFmt, 
                        KmerAlignFormat, 
                        KmerAlignDirFmt,
                        SidleReconFormat, 
                        SidleReconDirFormat,
                        ReconSummaryFormat,
                        ReconSummaryDirFormat
                        )


plugin.register_semantic_types(KmerMap, 
                               KmerAlignment,
                               SidleReconstruction,
                               ReconstructionSummary
                               )


plugin.register_semantic_type_to_format(FeatureData[KmerMap], 
                                        KmerMapDirFmt)


plugin.register_semantic_type_to_format(FeatureData[KmerAlignment], 
                                        KmerAlignDirFmt)


plugin.register_semantic_type_to_format(FeatureData[SidleReconstruction], 
                                        SidleReconDirFormat)


plugin.register_semantic_type_to_format(FeatureData[ReconstructionSummary], 
                                        ReconSummaryDirFormat)


importlib.import_module('q2_sidle._transformer')
