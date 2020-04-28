from qiime2.plugin import SemanticType, model
from q2_types.feature_data import FeatureData

KmerMap = SemanticType('KmerMap', variant_of=FeatureData.field['type'])

class KmerMapFormat(model.TextFileFormat):
    def validate(*args):
        pass

KmerMapsDirFmt = model.SingleFileDirectoryFormat(
    'KmerMapDirFmt', 'map.tsv', KmerMapFormat)

KmerAlignment = SemanticType('FeatureData', variant_of=FeatureData.field['type'])

class KmerAlignFormat(model.TextFileFormat):
    def validate(*args):
        pass

KmerAlignDirFmt = model.SingleFileDirectoryFormat(
    'KmerAlignDirFmt', 'align.tsv', KmerAlignFormat)
