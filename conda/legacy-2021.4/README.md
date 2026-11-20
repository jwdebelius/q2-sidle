# Legacy QIIME 2021.4 Conda recipe

This recipe packages the existing `0.1.0-beta` tag at commit `35589ef2df684df597369ed722dc31e84208a3b2`, preserving the runtime used by nf-core/ampliseq 2.11.0. Distribution version: `2021.2.dev0`; QIIME plugin version: `2020.08`. It does not change the current source or release version.

```
conda build conda/legacy-2021.4 --override-channels -c qiime2 -c conda-forge -c bioconda
```

The QIIME packages are provided by the `qiime2` channel; this is not a standard Bioconda recipe and does not add any publishing automation. It provides a tested package for a possible QIIME-channel distribution path, respecting the channel concern in #34.

Local build and fresh-environment plugin tests pass. The scientific smoke test uses the public fixtures in [Galaxy tools-iuc PR #8831](https://github.com/galaxyproject/tools-iuc/pull/8831):

```
python conda/legacy-2021.4/test_scientific_payload.py /path/to/tools-iuc/tools/qiime2_sidle/test-data /path/to/evidence
```

All 27 count cells match the retained fixture; all nine representative sequences have the expected 200-base length and content-derived identifiers. [Build logs, package SHA256, explicit environment and scientific comparison](https://github.com/nekrut/foundry-paper/tree/b651d08b/workflows/ampliseq/maintained-tools/2026-10-29/sidle-wrapper/conda-build-evidence) are public. No package has been uploaded to a channel.
