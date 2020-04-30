import csv

from qiime2.plugin import model, ValidationError
from q2_types.feature_data import FeatureData


class KmerMapFormat(model.TextFileFormat):
    def validate(*args):
        pass

KmerMapDirFmt = model.SingleFileDirectoryFormat(
    'KmerMapDirFmt', 'map.tsv', KmerMapFormat)


class KmerAlignFormat(model.TextFileFormat):
    def validate(*args):
        pass

KmerAlignDirFmt = model.SingleFileDirectoryFormat(
    'KmerAlignDirFmt', 'align.tsv', KmerAlignFormat)


class SidleReconFormat(model.TextFileFormat):
    def validate(*args):
        pass

SidleReconDirFormat = model.SingleFileDirectoryFormat(
    'SidleReconDirFormat', 'mapping.tsv', SidleReconFormat
    )

class ReconSummaryFormat(model.TextFileFormat):
    def validate(*args):
        pass

ReconSummaryDirFormat = model.SingleFileDirectoryFormat(
    'ReconSummaryDirFormat', 'mapping.tsv', ReconSummaryFormat
    )