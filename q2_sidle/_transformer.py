import pandas as pd
import dask.dataframe as dd
from qiime2 import Metadata

from q2_sidle import (KmerMapFormat, 
                      KmerAlignFormat, 
                      SidleReconFormat,
                      ReconSummaryFormat,
                      )
from q2_types.feature_data import  AlignedDNAFASTAFormat, DNAFASTAFormat
from q2_types.feature_data._transformer import _dnafastaformats_to_series
from q2_sidle.plugin_setup import plugin

@plugin.register_transformer
def _1(ff:KmerMapFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df = df[['db-seq', 'seq-name', 'kmer', 'region', 'fwd-primer', 
             'rev-primer', 'kmer-length']]
    df[['kmer-length']] = df[['kmer-length']].astype(int)
    df.sort_values(['db-seq', 'seq-name', 'kmer'], inplace=True)
    return df.set_index('db-seq')


@plugin.register_transformer
def _2(ff:KmerMapFormat) -> Metadata:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df = df[['db-seq', 'seq-name', 'kmer', 'region',
             'fwd-primer', 'rev-primer',  'kmer-length']]
    df[['kmer-length']] = df[['kmer-length']].astype(int)
    df.sort_values(['db-seq', 'seq-name', 'kmer'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.index = df.index.astype(str)
    df.index.set_names('id', inplace=True)
    return Metadata(df)


@plugin.register_transformer
def _3(ff:KmerMapFormat) -> dd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df = df[['db-seq', 'seq-name', 'kmer', 'region', 'fwd-primer', 
             'rev-primer', 'kmer-length']]
    df[['kmer-length']] = df[['kmer-length']].astype(int)
    df.sort_values(['db-seq', 'seq-name', 'kmer'], inplace=True)
    return dd.from_pandas(df.set_index('db-seq'), chunksize=50000)


@plugin.register_transformer
def _4(obj: pd.DataFrame) -> KmerMapFormat:
    ff = KmerMapFormat()
    obj.to_csv(str(ff), sep='\t')
    return ff


@plugin.register_transformer
def _5(ff:KmerAlignFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df = df[['kmer', 'asv', 'length', 'mismatch', 'region']]
    df[['mismatch', 'length']] = df[['mismatch', 'length']].astype(int)
    return df


@plugin.register_transformer
def _6(ff:KmerAlignFormat) -> Metadata:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df = df[['kmer', 'asv', 'length', 'mismatch', 'region']]
    df.index.set_names('feature-id', inplace=True)
    df.index = df.index.astype(str)
    df[['mismatch', 'length']] = df[['mismatch', 'length']].astype(int)
    return Metadata(df)

@plugin.register_transformer
def _7(ff:KmerAlignFormat) -> dd.DataFrame:
    df = dd.read_csv(str(ff), sep='\t', dtype=str)
    df = df[['kmer', 'asv', 'length', 'mismatch', 'region']]
    df[['mismatch', 'length']] = df[['mismatch', 'length']].astype(int)
    return df

@plugin.register_transformer
def _8(obj: pd.DataFrame) -> KmerAlignFormat:
    ff = KmerAlignFormat()
    obj.to_csv(str(ff), sep='\t', index=False)
    return ff

@plugin.register_transformer
def _15(obj: dd.DataFrame) -> KmerAlignFormat:
    ff = KmerAlignFormat()
    obj.to_csv(str(ff), sep='\t', index=False, single_file=True)
    return ff

@plugin.register_transformer
def _9(ff: SidleReconFormat) -> pd.Series:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    if df.index.name is None:
        df.set_index('db-seq', inplace=True)
    return df['clean_name']

@plugin.register_transformer
def _10(ff: SidleReconFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    if df.index.name is None:
        df.set_index('db-seq', inplace=True)
    df[[c for c in df.columns if ('length' in c)]] = \
        df[[c for c in df.columns if ('length' in c)]].astype(float)
    return df

@plugin.register_transformer
def _11(obj: pd.DataFrame) -> SidleReconFormat:
    ff = SidleReconFormat()
    obj.to_csv(str(ff), sep='\t', header=True)
    return ff

@plugin.register_transformer
def _12(ff:ReconSummaryFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', dtype=str)
    df.set_index('feature-id', inplace=True)
    df.drop('#q2:types', inplace=True)
    df[['num-regions', 'total-kmers-mapped']] = \
        df[['num-regions', 'total-kmers-mapped']].astype(float).astype(int)
    df[['mean-kmer-per-region', 'stdv-kmer-per-region']] = \
        df[['mean-kmer-per-region', 'stdv-kmer-per-region']].astype(float)
    return df

@plugin.register_transformer
def _13(ff:ReconSummaryFormat) -> Metadata:
    return Metadata.load(str(ff))

@plugin.register_transformer
def _14(obj:Metadata) -> ReconSummaryFormat:
    ff = ReconSummaryFormat()
    obj.save(str(ff))
    return ff
