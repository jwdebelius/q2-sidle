.. SIdle: SMURF Implementation Done to acceLerate Effeciency documentation master file, created by
	sphinx-quickstart on Wed Feb 12 10:46:59 2020.
	You can adapt this file completely to your liking, but it should at least
	contain the root `toctree` directive.

Welcome to Sidle
================
SMURF Implementation Done to acceLerate Effeciency

Sidle is a python version of the Short MUliple Reads Framework (SMURF) algorithm originally developed by `Fuks et al`_ (2018) with a novel tree-building solution. This allows the reconstruction of multiple short, fragmented amplicons against a known database to improve the resolution fo the reconstructed community over single amplicons. It can also be used for meta-analysis to provide improved asv-based resolution.

.. toctree::
	:maxdepth: 2
	:caption: Contents:

	install
	database_preparation
	read_preparation
	reconstruction
	parallel_processing


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. Links
.. _Fuks et al: https://www.ncbi.nlm.nih.gov/pubmed/29373999
