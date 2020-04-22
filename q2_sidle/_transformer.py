from qiime2 import Metadata

from q2_sidle import KmerMapFormat
from q2_sidle.plugin_setup import plugin

@plugin.register_transformer
def _1(ff:KmerMapFormat) -> Metadata:
	return Metadata.load(str(ff))

@plugin.register_transformer
def _2(obj: Metadata) -> KmerMapFormat:
	ff = KmerMapFormat()
	obj.save(str(ff))
	return ff
