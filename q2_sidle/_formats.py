import csv
import pandas as pd

from qiime2.plugin import model, ValidationError
from q2_types.feature_data import FeatureData


class KmerMapFormat(model.TextFileFormat):
    def validate(self, *args):
        pass

KmerMapDirFmt = model.SingleFileDirectoryFormat(
    'KmerMapDirFmt', 'kmer-map.tsv', KmerMapFormat)


class KmerAlignFormat(model.TextFileFormat):
    def validate(*args):
        pass

KmerAlignDirFmt = model.SingleFileDirectoryFormat(
    'KmerAlignDirFmt', 'kmer-align.tsv', KmerAlignFormat)


class SidleReconFormat(model.TextFileFormat):
    def validate(*args):
        pass

SidleReconDirFormat = model.SingleFileDirectoryFormat(
    'SidleReconDirFormat', 'sidle-reconstruction-mapping.tsv', 
    SidleReconFormat
    )

class ReconSummaryFormat(model.TextFileFormat):
    def validate(*args):
        pass

ReconSummaryDirFormat = model.SingleFileDirectoryFormat(
    'ReconSummaryDirFormat', 'sidle-summary.tsv', ReconSummaryFormat
    )

class AlignmentPosFormat(model.TextFileFormat):
    def validate(*args):
        pass

AlignmentPosDirFmt = model.SingleFileDirectoryFormat(
    'AlignmentPosDirFmt', 'position-summary.tsv', AlignmentPosFormat
    )
