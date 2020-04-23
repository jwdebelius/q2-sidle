.. SIdle: SMURF Implementation Done to acceLerate Effeciency documentation master file, created by
	sphinx-quickstart on Wed Feb 12 10:46:59 2020.
	You can adapt this file completely to your liking, but it should at least
	contain the root `toctree` directive.

Welcome to Sidle
================
SMURF Implementation Done to acceLerate Effeciency

Sidle is a python version of the Short MUliple Reads Framework (SMURF) algorithm originally developed by `Fuks et al`_ (2018). This allows the reconstruction of multiple short, fragmented amplicons against a known database to improve the resolution fo the reconstructed community over single amplicons. 

..with a novel tree build approach.
..It can also be used for meta analysis to provide improved asv-based resolution

.. toctree::
	:maxdepth: 2
	:caption: Contents:

	quickstart
	database_preperation
.. read_preperation
.. reconstruction
.. commands
	
.. parallel_processing
.. quick-start


Installation
------------

The current working version is a qiime2-based implementation. It requires that you have installed qiime2 according to the `installation instructions`_ and have activated the enviroment. Once inside the enviroment, sidle can be installed with the following commands.

.. code-block:: bash

   conda install dask regex
   conda install --channel conda-forge sparse
   pip install -e . --no-deps

.. You can test the installation by going to the `sidle/tests/` folder and running `nosetests`.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. Links
.. _Fuks et al: https://www.ncbi.nlm.nih.gov/pubmed/29373999
.. _installation instructions: https://docs.qiime2.org/2020.2/install/
.. .. _github: https://github.com/jwdebelius/sidle