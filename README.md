![Build Passing](https://travis-ci.com/jwdebelius/q2-sidle.svg?branch=2021.2-release) 
[![Coverage Status](https://coveralls.io/repos/github/jwdebelius/q2-sidle/badge.svg)](https://coveralls.io/github/jwdebelius/q2-sidle)
# q2-sidle
A QIIME 2 implementation of SMURF. 

Full documentation on [read the docs](https://q2-sidle.readthedocs.io/)

## Installation

q2-sidle requires a QIIME 2 enviroment. [Install QIIME2](https://docs.qiime2.org/2021.4/install/) according to the method that works best for your system.

Sidle depends on two conda libraries, `dask` and `regex`, as well as the [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt/) qiime2 plugin. To install the plugin, start with adding the addtional conda libraries:

```bash
conda install dask regex
conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict
```

Next, you will need the RESCRIPt and sidle plugins.

```bash
pip install git+https://github.com/bokulich-lab/RESCRIPt.git
pip install git+https://github.com/jwdebelius/q2-sidle
```

Finally, update your qiime libraries. 

```bash
qiime dev refresh-cache
```

## Tutorial

Please see our tutorial on [read the docs](https://q2-sidle.readthedocs.io/)


## Getting Help

Problem? Suggestion? Technical errors and user support requests can be filed on the [QIIME 2 Forum](https://forum.qiime2.org/).


## Citation

If you use Sidle in your reserach, please cite:

1. Debelius, J.W.; Robeson, M.; Lhugerth, L.W.; Boulund, F.; Ye, W.; Engstrand, L. "A comparison of approaches to scaffolding multiple regions along the 16S rRNA gene for improved resolution." Preprint in BioRxiv. doi: 10.1101/2021.03.23.436606

1. Fuks, G.; Elgart, M.; Amir, A.; Zeisel, A.; Turnbaugh, P.J., Soen, Y.; and Shental, N. (2018). "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**: 17. doi: 10.1186/s40168-017-0396-x
