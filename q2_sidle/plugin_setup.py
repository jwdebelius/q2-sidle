from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch)
from q2_types.feature_data import (FeatureData,
                                   Sequence,
                                   ) 
from q2_types.feature_table import (
                      FeatureTable, 
                      Frequency, 
                      )
#                       SampleData)
from q2_sidle import (KmerMap, 
                      KmerMapFormat, 
                      KmerMapDirFmt, 
                      KmerAlignment, 
                      KmerAlignFormat, 
                      KmerAlignDirFmt
                      )
import q2_sidle

plugin = Plugin(
    name='sidle',
    version='0.0.1',
    website='https://github.com/jwdebelius/q2-sidle',
    package='q2_sidle',
    description=('This plugin reconstructs a full 16s sequence from short '
                 'reads over a marker gene region using the Short MUltiple '
                 'Read Framework (SMURF) algorithm'),
    short_description='Plugin for taxonomic classification.',
    # citations=[citations['bokulich2018optimizing']]
)
plugin.register_formats(KmerMapFormat, 
                        KmerMapDirFmt, 
                        KmerAlignFormat, 
                        KmerAlignDirFmt
                        )
plugin.register_semantic_types(KmerMap, 
                               KmerAlignment)
plugin.register_semantic_type_to_format(FeatureData[KmerMap], 
                                        KmerMapDirFmt)
plugin.register_semantic_type_to_format(FeatureData[KmerAlignment], 
                                        KmerAlignDirFmt)

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
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
    }
)

plugin.methods.register_function(
    function=extract_regional_database,
    name='Extract and expand regional kmer database',
    description=('Performs in silico PCR to extract a region of the full '
                 'length sequence, and then expands degenerate sequences '
                 'and collapses duplicated sequences.'),
    inputs={
        'sequences': FeatureData[Sequence]
    },
    outputs=[
        ('collapsed_kmers', FeatureData[Sequence]), 
        ('kmer_map', FeatureData[KmerMap]),
    ],
    parameters={
        'fwd_primer': Str,
        'rev_primer': Str,
        'trim_length': Int,
        'region': Str,
        'primer_mismatch': Int,
        'trim_primers': Bool,
        'reverse_complement_rev': Bool,
        'chunk_size':  (Int % Range(1, None)),
        'n_workers': Int % Range(0, None),
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
        'fwd_primer': ('The forward primer used to amplify the region of '
                       'interest'),
        'rev_primer': ('The reverse primer used to amplify the region of '
                       "interest. If this is 3'-5' anti-sense primer, then "
                       'the `reverse-complement-rev` parameter should be '
                       'used'),
        'region': ('A unique description of the hypervariable region being '
                   'extracted. If no region is provided, the primers will'
                   ' be used'),
        'trim_length': ('The length of the extracted regional kmers.'),
        'primer_mismatch': ('The allowed mismatch between the database '
                            'sequence and the primer'),
        'trim_primers': ('Whether the primer should be trimmed from the '
                        'region. This is removed before `trim_length` is '
                        'applied'),
        'reverse_complement_rev': ('Indicates the reverse primer should be '
                                   'reverse complemented'),
        'chunk_size': ('The number of sequences to be analyzed in parallel '
                       'blocks'),
        'n_workers': ('The number of jobs to initiate. When `n_workers` is 0,'
                      ' the cluster will be able to access all avaliable'
                      ' resources.'),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
    },
)

plugin.methods.register_function(
    function=prepare_extracted_region,
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
        'chunk_size':  (Int % Range(1, None)),
        'n_workers': Int % Range(0, None),
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
        'chunk_size': ('The number of sequences to be analyzed in parallel '
                       'blocks'),
        'n_workers': ('The number of jobs to initiate. When `n_workers` is 0,'
                      ' the cluster will be able to access all avaliable'
                      ' resources.'),
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
    },
)

plugin.methods.register_function(
    function=align_regional_kmers,
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
        ('sequence_alignment', FeatureData[KmerAlignment]),
        ('discarded_sequences', FeatureData[Sequence]),

    ],
    parameters={
        'region': Str,
        'max_mismatch': Int % Range(1, None),
        'chunk_size':  (Int % Range(1, 5000)),
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
        'sequence_alignment': ('A mapping between the database kmer name and'
                               ' the asv'),
        'discarded_sequences': ('The sequences which could not be aligned to'
                                ' the database at the matching threshhold.'),
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
        'debug': ('Whether the function should be run in debug mode (without '
                  'a client) or not. `debug` superceeds all options'),
    },
)

