import csv

from qiime2.plugin import SemanticType
from q2_types.feature_data import FeatureData

KmerMap = SemanticType('KmerMap', variant_of=FeatureData.field['type'])

KmerAlignment = SemanticType('KmerAlignment', 
                             variant_of=FeatureData.field['type'])

SidleReconstruction = SemanticType('SidleReconstruction', 
                                   variant_of=FeatureData.field['type'])

ReconstructionSummary = SemanticType('ReconstructionSummary',
                                     variant_of=FeatureData.field['type'])

AlignmentPosSummary = SemanticType('AlignmentPosSummary',
                                   variant_of=FeatureData.field['type'])