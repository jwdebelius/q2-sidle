<!--![Build Passing](https://travis-ci.com/jwdebelius/q2-sidle.svg?branch=main) 
[![Coverage Status](https://coveralls.io/repos/github/jwdebelius/q2-sidle/badge.svg)](https://coveralls.io/github/jwdebelius/q2-sidle)
-->
# q2-sidle
A QIIME 2 implementation of SMURF. 

Full documentation on [read the docs](https://q2-sidle.readthedocs.io/)

## Installation

q2-sidle requires a QIIME 2 enviroment. [Install QIIME2](https://docs.qiime2.org/2020.8/install/) according to the method that works best for your system.

Sidle depends on two conda libraries, `dask` and `regex`, as well as the [RESCRIPt]() qiime2 plugin. To install the plugin, start with adding the addtional conda libraries:

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
