from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations, TypeMatch)
from q2_types.feature_data import (FeatureData,
                                   Sequence,
                                   ) 
# from q2_types.feature_table import (
#                       FeatureTable, 
#                       Frequency, 
#                       SampleData)
import q2_sidle

plugin = Plugin(
    name='sidle',
    version='0.0.1',
    website='https://github.com/jwdebelius/q2-sidle',
    package='q2_sidle',
    description=('This plugin reconstructs a full 16s sequence from short '
                 'reads over a marker gene region using the Short MUltiple '
                 'Read Framework (SMURF) algorith'),
    short_description='Plugin for taxonomic classification.',
    # citations=[citations['bokulich2018optimizing']]
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
        'max_degen': (Int % Range(1, None)),
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