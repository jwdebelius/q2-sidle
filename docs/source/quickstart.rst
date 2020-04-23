Quick Start
===========

This is intended to quickly provide a workflow for reconstruction. We assume that you do not have time to read additional documentation and wish to have commands thrust before you like a five course meal. If you are intested in something more complicated, please refer to the task-specific descriptions. This tutorial assumes that you are working with one of the pre-formatted databases. If 


Reconstruct Your Samples
------------------------



The next step is to prepare your reads. Alignment requires having high quality reads 


.. Demultiplex a sample into regions
.. +++++++++++++++++++++++++++++++++


.. Denoise your samples
.. ++++++++++++++++++++


.. Reconstruct the table
.. ---------------------


.. Perform regional alignment
.. ++++++++++++++++++++++++++


.. Reconstruct the table
.. +++++++++++++++++++++


.. Build your tree
.. ---------------


.. Have fun!
.. ---------


Format your database
--------------------

If you do not already have a database for reconstruction, you will need to prepare your own. (We recommend checking the Resources page to see if you can find one that suits your needs already.) If you do, skip down to the alignment section.

Otherwise, 

You will need:

	* A list of the forward and reverse files used to amplify each region of interest
	* A reference database of your choice
	* Patience to run the extraction


For this example, you can get the reference database files with

.. code-block:: bash
	
	wget X
	tar -xzf sidle_db_tutorial.tgz
	cd sidle_db_tutorial


Pre-filtering
+++++++++++++

You will significantly improve performance and memory consumption during extraction if you pre-filter your daatabase. The original SMURF paper recommended removing any sequence with more than 3 degeneracies, here we'll try 5.

.. code-block:: bash

    qiime sidle filter-degenerate-sequences \
     --i-sequences 88_otus.qza \
     --p-max-degen 5 \
     --o-filtered-sequences 88_otus-filtered.qza \

If you want to learn more about pre-filtering, please refer to the `pre filtering page`_. 


Extract the database
++++++++++++++++++++

We'll perform database extraction using a single step with the ``extract-regional-database`` command. This command will take a forward primer, reverse primer, and trimming length and return an extracted region of kmers ready for alignment with ASVs. For this example, we'll use the 74F (``TGGCGGACGGGTGAGTAA``)-315R(``CTGCTGCCTCCCGTAGGA``) 16s rRNA primers from the original SMURF paper [1]_. Let's trim the region to 100nt (``--p-trim-length``). We'll also include a description of the region used (``--p-region``). This can be the primers, a description of the database, or some other string that makes you happy. If you're using strings that make you happy, keep in mind that the region name must be unique for each region used in reconstruction. To make sure that we have a unique region name, let's use *smurf 74F-315R 100nt*.

.. .. code-block:: bash

	sidle extract-regional-database \
	 --i-sequences 88_otus-filtered.qza \
	 --p-fwd-primer TGGCGGACGGGTGAGTAA \
	 --p-rev-primer CTGCTGCCTCCCGTAGGA \
	 --p-trim-length 100 \
	 --p-region "smurf 74F-315R 100nt"
	 --o-collapsed-kmers 88_outs_74_315_100_seqs.qza \
	 --o-kmer-map 88_outs_74_315_100_maps.qza

The command returns two files: the collapsed kmers, which are the sequences we'll use for alignment (``--o-collapsed-kmer``) and a file describing the relationship between the original database and the extracted kmers which is used for reconstruction (``--o-kmer-map``). 

Explore more about `extracting and collapsing databases`_ to learn how you can work with a region you've already extracted. 

