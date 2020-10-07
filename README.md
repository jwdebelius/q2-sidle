![Build Passing](https://travis-ci.com/jwdebelius/q2-sidle.svg?branch=main)
[![Coverage Status](https://coveralls.io/repos/github/jwdebelius/q2-sidle/badge.svg)](https://coveralls.io/github/jwdebelius/q2-sidle)

# q2-sidle
A QIIME 2 implementation of SMURF. 

Full documentation on [read the docs](https://q2-sidle.readthedocs.io/)

## Installation

q2-sidle requires a QIIME 2 enviroment. [Install QIIME2](https://docs.qiime2.org/2020.8/install/) according to the method that works best for your system.
Additionally, you will need `dask` and `regex`. Install these using conda inside your QIIME 2 enviroment:

```bash
conda install dask regex
```

Finally, you will need the sidle plugin:

```bash
pip install git+https://github.com/jwdebelius/q2-sidle
```

## Tutorial

Please see our tutorial on [read the docs](https://q2-sidle.readthedocs.io/)
