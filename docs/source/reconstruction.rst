Sequence Reconstruction
=======================

The core of the SMURF algorithm is based on the kmer-based reconstruction of short regions into a full-length framework. Within Sidle, there are two steps in database reconstruction. First, ASVs are aligned on a regional basis to generate the local kmer-based alignment. Then, the full collection of sequences is assembled into a reconstructed table of counts. For this example, we’ll work with a small, entirely artifical subset of samples that are designed to run quickly.

If you've already done the database tutorial, make sure that you're in the ``sidle_tutorial`` directory.

.. code-block:: bash

	pwd

If you're new to the tutorial, you can make a new tutorila directory by running

.. code-block:: bash

	mkdir -p sidle_tutorial
	cd sidle-tutorial

Next, you will need to get the tutorial data. You will need to download three sets of files: `the sequencing data`_, `the alignments`_, and the `reconstruction files`_.

.. code-block:: bash

	wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/data.zip
	unzip data.zip
	wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/alignment.zip
	unzip alignment.zip
	wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/reconstruction.zip
	unzip reconstruction.zip


If you have not run the database tutorial, you will also want to get the
`prepared database`_.

.. code-block:: bash

	wget  https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/database.zip
	tar -xzf database.tgz

Regional alignment
------------------

The first step in reconstruction is to perform per-region alignment between the sequences and the database. We’ll do this with the ``align-regional-kmers`` command. We set the reference database that we extracted previously as the ``--kmer-db-fp``. The ASV represent sequences and are passed as the ``--rep-seq-fp``. Finally, we supply a regional definition. This should be the same as the region name that you gave when you extracted the kmers. In this case, the region name was “WonderWoman”.

Alignment is a pleasantly parallelizable problem, meaning that we can get improved performance by distributing the jobs. You can set this in almost any ``sidle`` command using the ``--p-n-workers`` parameter. (And you can learn more about parallel processing doc:`here <parallel_processing>`). For this example, we'll use 2 cores; if you have more avaliable you can or should use them.

.. code-block:: bash

	qiime sidle align-regional-kmers \
     --i-kmers database/sidle-db-wonder-woman-100nt-kmers.qza \
     --i-rep-seq data/wonder-woman-100nt-rep-set.qza \
     --p-region WonderWoman \
     --p-n-workers 2 \
     --o-regional-alignment alignment/wonder-woman-align-map.qza

This will output an alignment file and any ASV sequences which wouldn't be aligned to the database, for your own record keeping.

.. Note::

	If you are unsure of the region name or read length for your database kmers, you can always check the provenance by visiting `view.qiime2.org`_

Optionally, you can also modify parameters for the number of basepairs that differ between the reference and representative sequences (``--p-max-mismatch``); the original paper uses a mismatch of 2 with 130nt sequences.

You may find that if you have longer kmers, you might want to increase this parameter accordingly. A lower (more stringent) value will increase the number of discarded sequences, while a higher number may mean your matches are lower quality. You may find that if you have longer kmers, you may want to increase this parameter accordingly. A lower (more stringent) value will increase the number of discarded sequences, a higher number may mean your matches are lower quality.

Using the same parameters, you will need to align the other two regions.

.. code-block:: bash

	qiime sidle align-regional-kmers \
	 --i-kmers database/sidle-db-batman-100nt-kmers.qza \
	 --i-rep-seq data/batman-100nt-rep-set.qza \
	 --p-region Batman \
	 --p-n-workers 2 \
	 --o-regional-alignment alignment/batman-align-map.qza

	qiime sidle align-regional-kmers \
	 --i-kmers alignment/green-lantern-kmer-db.qza \
	 --i-rep-seq table/green-lantern-rep-seq.qza \
	 --p-region GreenLantern \
	 --p-n-workers 2 \
	 --o-regional-alignment alignment/green-lantern-align-map.qza

Now, you have all three local alignments prepared, you're ready to
reconstruct your table.

Table Reconstruction
--------------------

The table is reconstructed in three steps. First, the regional fragments get re-assembled into complete database sequences. Then, the relative abundance of the pooled counts gets computed through an optimization process. Finally, the relative abundance is used to reconstruct a table of counts.

The ``per-nucleotide-error`` is combined with the ``maximum-mismatch`` parameter from alignment to the probability that a sequence that differs from the reference. So, for instance, this algorithm allows a single ASV to be mapped to multiple sequences in the reference database. During reconstruction, the alignment mismatch, sequencing error, and relative abundance are combined to calculate the mapped abundance.

The ``min-abundance`` determines the relative abundance of a database sequence to be excluded during optimization. This is, to some degree, a function of the avaliable sequencing depth and the desired specificity of the fit.

Finally, let's plan on running the command in parallel, using the ``--p-n-workers`` flag; this is particularly useful in the per-sample reconstruction step. We'll use 2 workers in this tutorial, if you have more avaliable you may prefer that.

Now, let’s reconstruct the table, using the default settings.

.. code-block:: shell

    qiime sidle reconstruct-counts \
      --p-region WonderWoman \
        --i-kmer-map database/sidle-db-wonder-woman-100nt-map.qza \
        --i-regional-alignment alignment/wonder-woman-align-map.qza \
        --i-regional-table data/wonder-woman-100nt-table.qza \
      --p-region Batman \
        --i-kmer-map database/sidle-db-batman-100nt-map.qza \
        --i-regional-alignment alignment/batman-align-map.qza \
        --i-regional-table data/batman-100nt-table.qza \
      --p-region GreenLantern \
        --i-kmer-map database/sidle-db-green-lantern-100nt-map.qza \
        --i-regional-alignment alignment/green-lantern-align-map.qza \
        --i-regional-table data/green-lantern-100nt-table.qza \
      --p-n-workers 2 \
      --o-reconstructed-table reconstruction/league_table.qza \
      --o-reconstruction-summary reconstruction/league_summary.qza \
      --o-reconstruction-map reconstruction/league_map.qza

The command will produce a count table, a file containing details about the number of database kmers mapped to a region along with the ASV IDs, and a mapping that’s needed if you want to do taxonomic reconstruction.

Let’s take a look at the count table.

.. code-block:: shell

    qiime feature-table summarize \
     --i-table reconstruction/league_table.qza \
     --o-visualization reconstruction/league_table.qzv


You’ll notice that some of the feature IDs contain a ``|`` character, for example, ``1764594|195532|4471854``. This means the three databases sequences could not be resolved during the reconstruction, and so we assign the sequence to both regions. The more regions that are used in the reconstruction, the more likely you are to be able to accurately reconstuct the database sequences.

The second output is a summary. The summary can be used to evaluate the quality of the reconstruction; see the `original manuscript`_ [1]_ for more details. By default, the summary will consider degenerate kmers as unique sequences; you can change the behavior using the ``count-degenerates`` parameter; when False, kmers will only be counted if they belong to unique reference sequences. You can view the summary by tabulating the metadata.

.. code:: bash

    qiime metadata tabulate \
     --m-input-file reconstruction/league_summary.qza \
     --o-visualization reconstruction/league_summary.qzv


Let’s look at the information for the unresolved feature, ``133719|158591|190649``. How many regions has it found?

Taxonomic Reconstruction
------------------------

Now you have the table reconstructed, you’re ready to reconstruct the taxonomy to match. Specifcially, this process addresses cases where multiple database sequences cannot be untangled. The function takes the database map generated during reconstruction and the taxonomy associated with the database, and returns the reconstructed taxonomy.

There are three possible general cases for a set of shared sequences. First, they can share the full taxonomic string; second, they may differ at some point: or third, they may be same until one is missing an assignment. Let’s start with the simplest case. If we have two database sequences::

   1234    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum
   1235    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum

Then, when we reconstruct taxonomy, everything is the same and the final taxonomic label should be::

   1234 | 1235 k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum

There’s also thee possibility that sequences differ at some higher level, for example::

   1236    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum
   1237    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Roseburia; s__

In that case, the algorithm would keep the taxonomic assignment associated with the most recent common ancestor::

   1236 | 1237 k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia | g__Roseburia; g__Blautia | g__Rosburia

The ``--database`` parameter allows the user to select the type of database being used (``greengenes``, ``silva`` or ``none``). If the database is a defined database(``greengenes`` or ``silva``), some ad-hoc database cleaning will be performed automatically ✨, specifically with regard to the ``define-missing`` and ``ambiguity-handling`` parameters. For example, if a defined string is::

   k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Entrobacteriales; f__Enterobacteriaceae; g__; s__

Then, the new, cleaned string will be::

    k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Entrobacteriales; f__Enterobacteriaceae; g__unsp. f. Enterobacteriaceae; s__unsp. f. Enterobacteriaceae

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

The last step in reconstruction is to reconstruct fragments for the phylogenetic tree. Unfortunately, if the reference sequences cannot be resolved, the phylogenetic tree cannot simply be inherited from the database. So, we need to reconstruct a new phylognetic tree. We handle sequences in two ways.

1. Any database sequence which could full resolved can keep it’s position in the reference tree
2. Sequences which can’t be resolved need to handled somehow.

We could randomly select a sequence to map the reconstructed region to. However, that might not work when there are several sequences that got combined. So, instead, if we can’t resolve the database sequence, we calculate a concensus sequence from the combined data, extract them over the regions we were able to map, and then those consensus sequences can be inserted into a phylogenetic reference backbone using SEPP or something similar.

.. Note::

	Sucessful reconstruction requires that the ids in the database you used as your reference for reconstruction and the database you’re using for alignment are the same. Make sure that you are using the same database release version and the same level of sequence identity.

So, our first step is to reconstruct the consensus fragments from sequences that could not be resolved.

.. code-block:: shell

    qiime sidle reconstruct-fragment-rep-seqs \
     --p-region WonderWoman \
      --i-kmer-map database/sidle-db-wonder-woman-100nt-map.qza \
     --p-region Batman \
      --i-kmer-map database/sidle-db-batman-100nt-map.qza \
     --p-region GreenLantern \
      --i-kmer-map database/sidle-db-green-lantern-100nt-map.qza \
     --i-reconstruction-map reconstruction/league_map.qza \
     --i-reconstruction-summary reconstruction/league_summary.qza \
     --i-aligned-sequences database/sidle-db-aligned-sequences.qza \
     --o-representative-fragments reconstruction/league-rep-seq-fragments.qza

We can then insert the sequences into the reference tree. Let's first get the reference tree.

.. code-block:: shell

  wget \
    -O "sepp-refs-gg-13-8.qza" \
    "https://data.qiime2.org/2021.2/common/sepp-refs-gg-13-8.qza"

Then, we'll do the fragment insertion.

.. code-block:: shell

  qiime fragment-insertion sepp \
    --i-representative-sequences reconstruction/league-rep-seq-fragments.qza \
    --i-reference-database sepp-refs-gg-13-8.qza \
    --o-tree reconstruction/league-tree.qza \
    --o-placements reconstruction/league-placements.qza

Now, you're ready to analyze your data.

Next Steps: Analysis!
---------------------

You now have a reconstructed table, and associated taxonomy. Go forth and enjoy your analysis. The `QIIME 2 tutorials`_ offer some good options of downstream diversity and statistical analyses that can be done with this data.

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
    --i-rep-seq [ASV representative sequnces] \
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

* A phylogenetic tree can be reconstructed by first estimating the consensus fragments for the original sequences and then inserting them into a tree.
* See the `q2-fragment-insertion`_ documentation for more information

**Fragment reconstruction syntax**

.. code-block:: shell

  qiime sidle reconstruct-fragment-rep-seqs \
    --p-region [region 1 name] \
      --i-kmer-map [region 1 kmer map] \
      ... \
    --p-region [region n name] \
      --i-kmer-map [region n kmer map] \
     --i-reconstruction-map [reconstruction map] \
     --i-reconstruction-summary [reconstruction summary] \
     --i-aligned-sequences [aligned sequences] \
     --o-representative-fragments [reconstructed fragments]

**Example reconstruction syntax**

.. code-block:: shell

  qiime sidle reconstruct-fragment-rep-seqs \
   --p-region WonderWoman \
    --i-kmer-map database/sidle-db-wonder-woman-100nt-map.qza \
   --p-region Batman \
    --i-kmer-map database/sidle-db-batman-100nt-map.qza \
   --p-region GreenLantern \
    --i-kmer-map database/sidle-db-green-lantern-100nt-map.qza \
   --i-reconstruction-map reconstruction/league_map.qza \
   --i-reconstruction-summary reconstruction/league_summary.qza \
   --i-aligned-sequences database/sidle-db-aligned-sequences.qza \
   --o-representative-fragments reconstruction/league-rep-seq-fragments.qza

References
++++++++++

.. [1] Fuks, C; Elgart, M; Amir, A; et al (2018) "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**:17. doi: 10.1186/s40168-017-0396-x

.. links

.. _here: https://github.com/jwdebelius/q2-sidle/tree/main/docs/tutorial_data
.. _view.qiime2.org: https://view.qiime2.org
.. _absloute paths: https://www.linux.com/training-tutorials/absolute-path-vs-relative-path-linuxunix/
.. _original manuscript: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0396-x
.. _QIIME 2 tutorials: https://docs.qiime2.org/2021.2/tutorials/
.. _q2-fragment-insertion: https://docs.qiime2.org/2020.8/plugins/available/fragment-insertion/
.. _the sequencing data: https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/data.zip
.. _the alignments: https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/alignment.zip
.. _reconstruction files: https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/reconstruction.zip
.. _prepared database:  https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/database.zip
