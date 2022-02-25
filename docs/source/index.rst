.. SIdle: SMURF Implementation Done to acceLerate Effeciency documentation master file, created by
	sphinx-quickstart on Wed Feb 12 10:46:59 2020.
	You can adapt this file completely to your liking, but it should at least
	contain the root `toctree` directive.

Welcome to Sidle
================
SMURF Implementation Done to acceLerate Effeciency

Sidle is a python version of the Short MUltiple Reads Framework (SMURF) algorithm originally developed by `Fuks et al`_ (2018) with a novel tree-building solution. This allows the reconstruction of multiple short, fragmented amplicons against a known database to improve the resolutionf of the reconstructed community over single amplicons.

.. toctree::
	:maxdepth: 2
	:caption: Contents:

	install
	database_preparation
	read_preparation
	reconstruction
	parallel_processing


Citations
=========

If you use Sidle, please cite:

1. Debelius, J.W.; Robeson, M.; Lhugerth, L.W.; Boulund, F.; Ye, W.; Engstrand, L. "A comparison of approaches to scaffolding multiple regions along the 16S rRNA gene for improved resolution." Preprint in BioRxiv. doi: 10.1101/2021.03.23.436606

1. Fuks, G.; Elgart, M.; Amir, A.; Zeisel, A.; Turnbaugh, P.J., Soen, Y.; and Shental, N. (2018). "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**: 17. doi: 10.1186/s40168-017-0396-x


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. Links
.. _Fuks et al: https://www.ncbi.nlm.nih.gov/pubmed/29373999
