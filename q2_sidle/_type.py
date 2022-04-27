import csv

from qiime2.plugin import SemanticType
from q2_types.feature_data import FeatureData
from q2_types.sample_data import SampleData

KmerMap = SemanticType('KmerMap', variant_of=FeatureData.field['type'])

KmerAlignment = SemanticType('KmerAlignment', 
                             variant_of=FeatureData.field['type'])

SidleReconstruction = SemanticType('SidleReconstruction', 
                                   variant_of=FeatureData.field['type'])

ReconstructionSummary = SemanticType('ReconstructionSummary',
                                     variant_of=FeatureData.field['type'])

AlignmentLedger = SemanticType('AlignmentLedger', 
                               variant_of=SampleData.field['type'])