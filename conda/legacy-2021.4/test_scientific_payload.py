"""Run with the built package environment and the public IUC SIDLE fixtures."""
import argparse,hashlib,json
from pathlib import Path
import biom
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins.sidle.methods import trim_dada2_posthoc
p=argparse.ArgumentParser();p.add_argument('fixtures',type=Path);p.add_argument('evidence',type=Path);a=p.parse_args();a.evidence.mkdir(parents=True,exist_ok=True)
table=Artifact.load(str(a.fixtures/'3_table.qza'));seqs=Artifact.load(str(a.fixtures/'3_rep-seqs.qza'))
result=trim_dada2_posthoc(table=table,representative_sequences=seqs,trim_length=200,hashed_feature_ids=True)
actual=result.trimmed_table.view(biom.Table);expected=Artifact.load(str(a.fixtures/'3_200_table.qza')).view(biom.Table)
# Fixture table is the independently retained result used by align/count tests.
features=sorted(expected.ids(axis='observation'));samples=sorted(expected.ids());assert set(actual.ids(axis='observation'))==set(features);assert set(actual.ids())==set(samples)
actual=actual.sort_order(features,axis='observation').sort_order(samples,axis='sample');expected=expected.sort_order(features,axis='observation').sort_order(samples,axis='sample');assert (actual.matrix_data.toarray()==expected.matrix_data.toarray()).all()
sequences=result.trimmed_representative_sequences.view(pd.Series);assert set(sequences.index)==set(features)
assert all(len(str(seq))==200 and hashlib.md5(str(seq).encode()).hexdigest()==name for name,seq in sequences.items())
proof={'status':'PASS','action':'sidle trim-dada2-posthoc','features':len(features),'samples':len(samples),'count_cells_matching':len(features)*len(samples),'sequences_verified_by_content_md5_and_length':len(features),'sequence_length':200,'inputs':{name:hashlib.sha256((a.fixtures/name).read_bytes()).hexdigest() for name in ['3_table.qza','3_rep-seqs.qza','3_200_table.qza']},'scope':'Fresh Conda package environment, not the prior nf-core container. All output count cells compared with the retained fixture; every sequence identity and length checked.'}
(a.evidence/'scientific-smoke.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
