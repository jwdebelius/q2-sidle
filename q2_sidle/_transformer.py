import pandas as pd
from qiime2 import Metadata

from q2_sidle import (KmerMapFormat, 
                      KmerAlignFormat, 
                      SidleReconFormat,
                      ReconSummaryFormat,
                      )
from q2_sidle.plugin_setup import plugin

@plugin.register_transformer
def _1(ff:KmerMapFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    return df.set_index('db_seq')

@plugin.register_transformer
def _2(obj: pd.DataFrame) -> KmerMapFormat:
    ff = KmerMapFormat()
    obj.to_csv(str(ff), sep='\t')
    return ff

@plugin.register_transformer
def _3(ff:KmerAlignFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    return df

@plugin.register_transformer
def _4(obj: pd.DataFrame) -> KmerAlignFormat:
    ff = KmerAlignFormat()
    obj.to_csv(str(ff), sep='\t')
    return ff

@plugin.register_transformer
def _5(ff: SidleReconFormat) -> pd.Series:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df.set_index('clean_name', inplace=True)
    return df['db_seq']

@plugin.register_transformer
def _6(obj: pd.Series) -> SidleReconFormat:
    ff = SidleReconFormat()
    obj.to_csv(str(ff), sep='\t', header=True)
    return ff

@plugin.register_transformer
def _7(ff:ReconSummaryFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df.set_index('feature-id', inplace=True)
    df.drop('#q2:types', inplace=True)
    df[['num_regions', 'total_kmers_mapped']] = \
        df[['num_regions', 'total_kmers_mapped']].astype(float).astype(int)
    df[['mean_kmer_per_region', 'stdv_kmer_per_region']] = \
        df[['mean_kmer_per_region', 'stdv_kmer_per_region']].astype(float)
    return df

@plugin.register_transformer
def _8(obj:pd.DataFrame) -> ReconSummaryFormat:
    ff = ReconSummaryFormat()
    obj.to_csv(ff, sep='\t')
    return ff

@plugin.register_transformer
def _9(ff:ReconSummaryFormat) -> Metadata:
    return qiime2.Metadata.load(str(ff))

@plugin.register_transformer
def _10(obj:Metadata) -> ReconSummaryFormat:
    ff = ReconSummaryFormat()
    obj.save(str(ff))
    return ff
