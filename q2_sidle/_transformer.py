from qiime2 import Metadata

from q2_sidle import KmerMapFormat, KmerAlignDirFmt
from q2_sidle.plugin_setup import plugin

@plugin.register_transformer
def _1(ff:KmerMapFormat) -> Metadata:
    return Metadata.load(str(ff))

@plugin.register_transformer
def _2(ff:KmerMapFormat) -> pd.DataFrame:
    return Metadata.load(str(ff)).to_dataframe()

@plugin.register_transformer
def _3(obj: Metadata) -> KmerMapFormat:
    ff = KmerMapFormat()
    obj.save(str(ff))
    return ff

@plugin.register_transformer
def _4(ff:KmerAlignDirFmt) -> Metadata:
    return Metadata.load(str(ff))