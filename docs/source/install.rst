Installation
============

The current working version is a qiime2-based implementation. It requires that you have installed qiime2 according to the `installation instructions`_ and have activated the environment. In addition to vanilla qiime2, you will also need to install `RESCRIPt`_; these installation instructions will guide you through that process.

.. note::
	Sidle has been tested against qiime2-2020.11 and later; it may not function with earlier versions of qiime2.


One you have activated your qiime2 enviroment, sidle can be installed with the following commands.

.. code-block:: bash
	
	conda install dask
	conda install -c conda-forge -c bioconda -c qiime2 -c defaults xmltodict
	pip install git+https://github.com/bokulich-lab/RESCRIPt.git
	pip install git+https://github.com/jwdebelius/q2-sidle
	qiime dev refresh-cache

You can test the installation by running

.. code-block:: bash
	
	qiime 

You should see a print out which includes ``sidle``.

Now, you're all ready to start using Sidle.
	
.. Now, you're read to analyze your data. We recommend follow
.. 
.. starting with the :ref:`quickstart tutorial <quickstart>` to start doing regional alignment on a pre-prepared database.

.. _installation instructions: https://docs.qiime2.org/2021.4/install/
.. .. _github: https://github.com/jwdebelius/sidle
.. _RESCRIPt: https://github.com/bokulich-lab/RESCRIPt

