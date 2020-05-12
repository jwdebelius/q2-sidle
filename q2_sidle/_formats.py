import csv
import pandas as pd

from qiime2.plugin import model, ValidationError
from q2_types.feature_data import FeatureData


class KmerMapFormat(model.TextFileFormat):
    def validate(self, *args):
        col_set = set(['db-seq', 'seq-name', 'kmer', 'region', 
                      'fwd-primer', 'rev-primer', 'kmer-length'])
        map_ = pd.read_csv(str(self), dtype=str, sep='\t')
        if set(map_.columns) != col_set:
            raise ValidationError('The KmerMap does not contain '
                                  'the correct columns')
        try:
            map_['kmer-length'].astype(float)
        except:
            raise ValidationError('The kmer-length column must be numeric')

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