Installation
============

The current working version is a qiime2-based implementation. It requires that you have installed qiime2 according to the `installation instructions`_ and have activated the environment. 

Once inside the environment, sidle can be installed with the following commands.

.. code-block:: bash
	
	conda install dask regex
	pip install git+https://github.com/jwdebelius/q2-sidle
	qiime dev refresh-cache

You can test the installation by running

.. code-block:: bash
	
	qiime 

You should see a print out which includes ``sidle``.

Now, you're all ready to start using Sidle.
	
.. Now, you're read to analyze your data. We recommend followi
.. 
.. starting with the :ref:`quickstart tutorial <quickstart>` to start doing regional alignment on a pre-prepared database.

.. _installation instructions: https://docs.qiime2.org/2020.2/install/
.. .. _github: https://github.com/jwdebelius/sidle
	

