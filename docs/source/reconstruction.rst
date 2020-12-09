Sequence Reconstruction
=======================


The core of the SMURF algorithm is based on the kmer-based reconstruction of short regions into a full length framework. Within Sidle, there are two steps in database reconstruction. First, ASVs are aligned on a regional basis to generate the local kmer-based alignment. Then, the full collection of sequences is assembled into a reconstructed table of counts. For this example, we’ll work with a small, entirely artificial subset of samples that are designed to run quickly.

You can get the tutorial data by running 

.. code-block:: bash
	
    mkdir -p sidle_tutorial
    cd sidle_tutorial
    pwd
    
    wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/data.zip
    unzip data.zip
    wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/alignment.zip
    unzip alignment.zip
    wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/reconstruction.zip
    unzip reconstruction.zip

If you have not run the database tutorial, you will also want to get the database data.

.. code-block:: bash
	
	wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/database.zip
	unzip database.zip

Regional alignment
------------------

The first step in reconstruction is to perform per-region alignment between the sequences and the database. We’ll do this with the ``align-regional-kmers`` command. We set the reference database that we extracted previously as the (``--kmer-db-fp``). The ASV represent sequences are passed as the ``--rep-seq-fp``. Finally, we supply a regional definition. This should be the same as the region name that you gave when you extracted the kmers. In this case, the region name was "WonderWoman".

.. code-block:: bash
	
	qiime sidle align-regional-kmers \
     --i-kmers database/sidle-db-wonder-woman-100nt-kmers.qza \
     --i-rep-seq data/wonder-woman-100nt-rep-set.qza \
     --p-region WonderWoman \
     --o-regional-alignment alignment/wonder-woman-align-map.qza

This will output an alignment file and any ASV sequences which wouldn't be aligned to the database, for your own record keeping.

.. Note::

	If you are unsure of the region name or read length for your database kmers, you can always check the provenance by visiting `view.qiime2.org`_

Optionally, you can also modify parameters for the number of base pairs that differ between the reference and representative sequences (``--p-max-mismatch``); the original paper uses a mismatch of 2 with 130nt sequences. You may find that if you have longer kmers, you may want to increase this parameter accordingly. A lower (more stringent) value will increase the number of discarded sequences, a higher number may mean your matches are lower quality.

Using the same parameters, you will need to align the other two regions.

.. code-block:: bash
	
	qiime sidle align-regional-kmers \
	 --i-kmers database/sidle-db-batman-100nt-kmers.qza \
	 --i-rep-seq data/batman-100nt-rep-set.qza \
	 --p-region Batman \
	 --o-regional-alignment alignment/batman-align-map.qza

	qiime sidle align-regional-kmers \
	 --i-kmers alignment/green-lantern-kmer-db.qza \
	 --i-rep-seq table/green-lantern-rep-seq.qza \
	 --p-region GreenLantern \
	 --o-regional-alignment alignment/green-lantern-align-map.qza

Now, you have all three local alignments prepared, you’re ready to reconstruct your table.

Table Reconstruction
--------------------

The table is reconstructed in 3 parts. First, the regional fragments get re-assembled into complete database sequences. During this process, the algorithm will calculate how well a sequence can be reconstructed by determining the average number of kmers per region that could be mapped. This can be calculated using degenerate sequences in the reference as unique (``count-degenerates``) or treating all variations on the original database sequence as the same (``no-count-degenerates``).

Then, the relative abundance of the pooled counts gets computed through an optimization process. Optimized builds on the per-nucleotide error rate (``per-nucleotide-error``) and the number of nucleotides difference between the reference sequence and mapped ASV (calculated during alignment). Fuks et al [1]_ reported the per-nucleotide error did not affect their reconstruction substantially, therefore we have used their value by default. However, interested parties may wish to explore their supplementary materials. The ``min-relative-abundance`` sets a threshold below which we assume an assigned read is an error. This is somewhat dependent on the sequencing depth.

Finally, the relative abundance is used to reconstruct the table of counts by assigning counts from reads to the table. While this is somewhat artificial due to the differences in regional depth and coverage, the ``region-normalize`` will determine the observed depth. The default, ``average`` will calculate the depth as the average of the sequencing counts from the original tables and the final table will have a depth close to the regional input tables. ``weighted`` will calculate the depth as the sum of the sequencing counts from the original table, so the final table will have more counts than the original table.

We will perform reconstructing using the defaults.

.. code-block:: shell
	
    qiime sidle reconstruct-counts \
     --p-region WonderWoman \
      --i-kmer-map database/sidle-db-wonder-woman-100nt-map.qza \
      --i-regional-alignment alignment/wonder-woman-align-map.qza \
      --i-regional-table data/data/wonder-woman-100nt-table.qza \
     --p-region Batman \
      --i-kmer-map database/sidle-db-batman-100nt-map.qza \
      --i-regional-alignment alignment/batman-align-map.qza \
      --i-regional-table data/batman-100nt-table.qza \
     --p-region GreenLantern \
      --i-kmer-map database/sidle-db-green-lantern-100nt-map.qza \
      --i-regional-alignment alignment/green-lantern-align-map.qza \
	  --i-regional-table data/green-lantern-100nt-table.qza \
     --o-reconstructed-table reconstruction/league_table.qza \
     --o-reconstruction-summary reconstruction/league_summary.qza \
     --o-reconstruction-map reconstruction/league_map.qza

The command will produce a count table, a file containing details about the number of database kmers mapped to a region along with the ASV IDs, and a mapping that’s needed if you want to do taxonomic reconstruction.

Let’s take a look at the count table.

.. code-block:: shell
	
    qiime feature-table summarize \
     --i-table reconstruction/league_table.qza \
     --o-visualization reconstruction/league_table.qzv


You’ll notice that some of the feature IDs contain a ``|`` character, for example, ``1764594|195532|4471854``. This means the two databases sequences could not be resolved during the reconstruction, and so we assign the sequence to both regions. The more regions that are used in the reconstruction, the more likely you are to be able to accurately reconstruct the database sequences.

The second output is a summary. The summary can be used to evaluate the quality of the reconstruction; see the `original manuscript`_ [1]_ for more details. By default, the summary will consider degenerate kmers as unique sequences; you can change the behavior using the ``count-degenerates`` parameter; when False, kmers will only be counted if they belong to unique reference sequences. You can view the summary by tabulating the metadata.

.. code:: bash

    qiime metadata tabulate \
     --m-input-file reconstruction/league_summary.qza \
     --o-visualization reconstruction/league_summary.qzv


Let’s look at the information for the unresolve feature, ``1764594|195532|4471854``. How many regions is it found in?

Taxonomic Reconstruction
------------------------

Now you have the table reconstructed, you’re ready to reconstruct the taxonomy to match. Specifically, this process addresses cases where multiple database sequences cannot be untangled. The function takes the database map generated during reconstruction and the taxonomy associated with the database, and returns the reconstructed taxonomy.

There are three possible general cases for a set of shared sequences. First, they can share the full taxonomic string, second they may differ at some point, or third, they may be the same until one is missing an assignment. Let’s start with the simplest case. If we have two database sequences::

   1234    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum
   1235    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum

Then, when we reconstruct taxonomy, everything is the same and the final taxonomic label should be::

   1234 | 1235 k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum

There’s also the possibility that sequences differ at some higher level, for example::

   1236    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum
   1237    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Roseburia; s__

In that case, the algorithm would keep the taxonomic assignment associated with the most recent common ancestor::

   1236 | 1237 k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia | g__Roseburia; g__Blautia | g__Rosburia

If the ``--database`` parameter allows the user to select the type of database being used (``greengenes``, ``silva`` or ``none``). If the database is a defined database(``greengenes`` or ``silva``), some ad-hoc database cleaning will be performed automatically ✨. For example, if a defined string is::

   k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Entrobacteriales; f__Enterobacteriaceae; g__; s__

Then, the new, cleaned string will be::

    k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Entrobacteriales; f__Enterobacteriaceae; g__unsp. f. Enterobacteriaceae; s__unsp. f. Enterobacteriaceae

The ``--database`` parameter allows the user to select the type of database being used (``greengenes``, ``silva`` or ``none``). If the database is a defined database(``greengenes`` or ``silva``), some ad-hoc database cleaning will be performed and uncultured sequences will be handled, specifically with regard to the ``define-missing`` and ``ambiguity-handling`` parameters.

Our database is a subset of the greengenes database, so let’s specify that we used the greengenes database and inherit the missing strings.

.. code-block:: shell
    
    qiime sidle reconstruct-taxonomy \
     --i-reconstruction-map reconstruction/league_map.qza \
     --i-taxonomy database/sidle-db-taxonomy.qza \
     --p-database 'greengenes' \
     --p-define-missing 'inherit' \
     --o-reconstructed-taxonomy reconstruction/league_taxonomy.qza

You can check the taxonomic reconstruction by tabulating the taxonomy.

.. code-block:: shell

    qiime metadata tabulate \
     --m-input-file reconstruction/league_taxonomy.qza \
     --o-visualization reconstruction/league_taxonomy.qzv

What’s the taxonomy assignment for ``1764594|195532|4471854``?


Reconstructing the Phylogenetic Tree
------------------------------------

The last step in reconstruction is to reconstruct fragments for the phylogenetic tree. Unfortunately, if the reference sequences cannot be resolved, the phylogenetic tree cannot simply be inherited from the database. So, we need to reconstruct a new phylogenetic tree. We handle sequences in two ways.

1. Any database sequence which could full resolved can keep it’s position in the reference tree
2. Sequences which can’t be resolved need to be handled somehow.

We could randomly select a sequence to map the reconstructed region to. However, that might not work when there are several sequences that got combined. So, instead, if we can’t resolve the database sequence, we calculate a consensus sequence from the combined data, extract them over the regions we were able to map, and then those consensus sequences can be inserted into a phylogenetic reference backbone using SEPP or something similar.

.. Note::

	Successful reconstruction requires that the ids in the database you used as your reference for reconstruction and the database you’re using for alignment are the same. Make sure that you are using the same database release version and the same level of sequence identity

So, our first step is to reconstruct the consensus fragments from sequences that could not be resolved.

.. code-block:: shell

    qiime sidle reconstruct-fragment-rep-seqs \
    --p-region WonderWoman \
      --i-regional-alignment alignment/wonder-woman-align-map.qza \
     --p-region Batman \
      --i-regional-alignment alignment/batman-align-map.qza \
     --p-region Green-Lantern \
      --i-regional-alignment alignment/green-lantern-align-map.qza \
     --i-reconstruction-map reconstruction/league_map.qza \
     --i-reconstruction-summary reconstruction/league_summary.qza \
     --i-aligned-sequences database/sidle-db-aligned-sequences.qza \
     --o-representative-fragments reconstruction/league-rep-seq-fragments.qza

We can then insert the sequences into the reference tree. Let's first get the reference tree.

.. code-block:: shell

	wget \
	 -O "sepp-refs-gg-13-8.qza" \
	 "https://data.qiime2.org/2020.11/common/sepp-refs-gg-13-8.qza"

Then, we'll do the fragment insertion. 

.. code-block:: shell

    qiime fragment-insertion sepp \
     --i-representative-sequences reconstruction/league-rep-seq-fragments.qza \
     --i-reference-database ../../../medda-bench/simulations/refs/greengenes/sepp-refs-gg-13-8.qza \
     --o-tree reconstruction/league-tree.qza \
     --o-placements reconstruction/league-placements.qza

Now, you're ready to analyze your data.

Next Steps: Analysis!
---------------------

You now have a reconstructed table, and associated taxonomy. Go forth and enjoy your analysis. The `QIIME 2 tutorials`_ offer some good options of downstream diversity and statistical analyses, as does the `qiime2 library`_ The `qiime2R`_ package allows easy import of qiime2 Artifacts into R.

TL;DR Reconstruction
--------------------

Regional Alignment Commands
+++++++++++++++++++++++++++

* The region name for the alignment **must match** the region name used for building the kmer map
* Kmers and representative sequences must be the same length
* This step is performed on a per-region basis

**Syntax**

.. code-block:: bash
	
	qiime sidle align-regional-kmers \
	 --i-kmers [kmer sequences from extracted database] \
	 --i-rep-seq [ASV representative sequences] \
	 --p-region [Region name] \
	 --o-regional-alignment [regional alignment]

**Example**

.. code-block:: bash
	
	qiime sidle align-regional-kmers \
	 --i-kmers wonderwoman-kmer-db.qza \
	 --i-rep-seq wonderwoman-rep-seq.qza \
	 --p-region WonderWoman \
	 --o-regional-alignment wonderwoman-align-map.qza

Reconstructing the Table
++++++++++++++++++++++++

* Make sure your region names match between the alignment artifact, the database kmer map, and the ``region`` parameter.
* ``count-degenerates`` will control how the summary describes differences in the sequences
* ``region-normalize`` will affect how many counts are assigned in the final table

**Syntax**

For *n* regions

.. code-block:: bash

	qiime sidle reconstruct-counts \
     --p-region [region 1 name] \
      --i-kmer-map [region 1 kmer map] \
      --i-regional-alignment [region 1 alignment] \
      --i-regional-table [region 1 counts table] \
      ... \
   	  --p-region [region n name] \
      --i-kmer-map [region n kmer map] \
      --i-regional-alignment [region n alignment] \
      --i-regional-table [region n counts table] \
     --o-reconstructed-table [reconstructed table] \
     --o-reconstruction-summary [reconstruction summary] \
     --o-reconstruction-map [reconstructed database map]

**Example**

.. code-block:: bash

	qiime sidle reconstruct-counts \
     --p-region WonderWoman \
      --i-kmer-map database/sidle-db-wonder-woman-100nt-map.qza \
      --i-regional-alignment alignment/wonder-woman-align-map.qza \
      --i-regional-table data/data/wonder-woman-100nt-table.qza \
     --p-region Batman \
      --i-kmer-map database/sidle-db-batman-100nt-map.qza \
      --i-regional-alignment alignment/batman-align-map.qza \
      --i-regional-table data/batman-100nt-table.qza \
     --p-region GreenLantern \
      --i-kmer-map database/sidle-db-green-lantern-100nt-map.qza \
      --i-regional-alignment alignment/green-lantern-align-map.qza \
	  --i-regional-table data/green-lantern-100nt-table.qza \
     --o-reconstructed-table reconstruction/league_table.qza \
     --o-reconstruction-summary reconstruction/league_summary.qza \
     --o-reconstruction-map reconstruction/league_map.qza

Reconstructing taxonomy
+++++++++++++++++++++++

* A database specification is required 

**Syntax**

.. code-block:: bash

	qiime sidle reconstruct-taxonomy \
	 --i-reconstruction-map [reconstruction map] \
	 --i-taxonomy [taxonomy path] \
	 --p-database [database name] \
	 --o-reconstructed-taxonomy [reconstructed taxonomy]

**Example**

.. code-block:: bash

	qiime sidle reconstruct-taxonomy \
	 --i-reconstruction-map reconstruction/league_map.qza \
	 --i-taxonomy database/sidle-db-taxonomy.qza \
	 --p-database 'greengenes' \
	 --p-define-missing 'inherit' \
	 --o-reconstructed-taxonomy reconstruction/league_taxonomy.qza

Reconstructing the Tree
+++++++++++++++++++++++

* A phylogenetic tree can be reconstructed by first, estimating the consensus fragments for the original sequences and then inserting them into a tree.
* See the `q2-fragment-insertion`_ documentation for more information

**Fragment reconstruction syntax**

.. code-block:: shell
	
	qiime sidle reconstruct-fragment-rep-seqs \
	 --i-reconstruction-map [reconstruction map] \
	 --i-reconstruction-summary [reconstruction summary] \
	 --i-aligned-sequences [aligned sequences] \
	 --m-manifest-file [manifest] \
	 --o-representative-fragments [consensus fragments]

**Example reconstruction syntax**

.. code-block:: shell
	
	qiime sidle reconstruct-fragment-rep-seqs \
	 --i-reconstruction-map reconstruction/league_map.qza \
	 --i-reconstruction-summary reconstruction/league_summary.qza \
	 --i-aligned-sequences database/sidle-db-aligned-sequences.qza \
	 --m-manifest-file manifest.txt \
	 --o-representative-fragments reconstruction/league-rep-seq-fragments.qza

References
++++++++++

.. [1] Fuks, C; Elgart, M; Amir, A; et al (2018) "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**:17. doi: 10.1186/s40168-017-0396-x

.. links

.. _here: https://github.com/jwdebelius/q2-sidle/tree/main/docs/tutorial_data
.. _view.qiime2.org: https://view.qiime2.org
.. _qiime2 library: https://library.qiime2.org/
.. _qiime2R: https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121
.. _QIIME 2 tutorials: https://docs.qiime2.org/2020.6/tutorials/
.. _original manuscript: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0396-x
.. _q2-fragment-insertion: https://docs.qiime2.org/2020.8/plugins/available/fragment-insertion/



