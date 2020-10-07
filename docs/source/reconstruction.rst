Sequence Reconstruction
=======================

The core of the SMURF algorithm is based on the kmer-based reconstruction of short regions into a full length framework. Within Sidle, there are two steps in database reconstruction. First, ASVs are aligned on a regional basis to generate the local kmer-based alignment. Then, the full collection of sequences is assembled into a reconstructed table of counts. For this example, we’ll work with a small, entirely artifical subset of samples that are designed to run quickly.

You can get the tutorial data `here`_ or by running 

.. code-block:: bash
	
    mkdir -p sidle_tutorial/alignment
    cd sidle_tutorial
    wget https://github.com/jwdebelius/q2-sidle/blob/main/docs/tutorial_data/data.tgz
    wget https://github.com/jwdebelius/q2-sidle/blob/main/docs/tutorial_data/.tgz


If you have not run the database tutorial, you will also want to get the
database data.

.. code-block:: bash
	
	wget 

Regional alignment
------------------

The first step in reconstruction is to perform per-region alignment between the sequences and the database. We’ll do this with the ``align-regional-kmers`` command. We set the reference database that we extracted previously as the (``--kmer-db-fp``). The ASV represent sequences are passed as the ``--rep-seq-fp``. Finally, we supply a regional defination. This should be the same as the region name that you gave when you extracted the kmers. In this case, the region name was “WonderWoman”.

.. code-block:: bash
	
	qiime sidle align-regional-kmers \
     --i-kmers database/sidle-db-wonder-woman-100nt-kmers.qza \
     --i-rep-seq data/wonder-woman-100nt-rep-set.qza \
     --p-region WonderWoman \
     --o-regional-alignment alignment/wonder-woman-align-map.qza \
     --o-discarded-sequences alignment/wonder-woman-discarded-sequences.qza 

This will output an alignment file and any ASV sequences whcih wouldn't be aligned to the database, for your own record keeping.

.. Note::

	If you are unsure of the region name or read length for your database kmers, you can always check the provenance by visiting `view.qiime2.org`_

Optionally, you can also modify parameters for the number of basepairs that differ between the reference and representative sequences (``--p-max-mismatch``); the original paper uses a mismatch of 2 with 130nt sequences.

You may find that if you have longer kmers, you may want to increase this parameter accordingly. A lower (more stringent) value will increase the number of discarded sequences, a higher number may mean your matches are lower quality.

Using the same parameters, you will need need to align the other two regions.

.. code-block:: bash
	
	qiime sidle align-regional-kmers \
	 --i-kmers database/sidle-db-batman-100nt-kmers.qza \
	 --i-rep-seq data/batman-100nt-rep-set.qza \
	 --p-region Batman \
	 --o-regional-alignment alignment/batman-align-map.qza \
	 --o-discarded-sequences alignment/batman-discarded-sequences.qza \

	qiime sidle align-regional-kmers \
	 --i-kmers alignment/region3-kmer-db.qza \
	 --i-rep-seq table/region3-rep-seq.qza \
	 --p-region Green-Lantern \
	 --o-regional-alignment region3-align-map.qza \
	 --o-discarded-sequences region3-discarded-sequences.qza

Now, you have all three local alignments prepared, you’re ready to
reconstruct your table.

Table Reconstruction
--------------------

The table is reconstucted in 3 parts. First, the regional fragments get re-assembled into complete database sequences. Then, the relative abundance of the pooled counts gets computed through an optimization process. Finally, the relative abundance is used to reconstruct table of counts.

Manifest
++++++++

To do this, the function needs to bring together the regional pieces
(database map, alignment, and table) to compute a pooled,
region-normalized table. We map together regions using a manifest
format, a tab-seperated (``\t``) file.

Your manifest should contain five, case-sensitive columns:

-  **id** – a string naming the region **This must match the name given when you built the database and during kmer alignment**
-  **region-order** – a number that indicates the order in which the regions appear along the marker gene
-  **kmer-map** – the path to the mapping between original names in the database and the kmer sequences. This file is produced by the ``prepare-extracted-region`` command and should be semantic type ``FeatureData[KmerMap]``
-  **alignment-map** – the path to the mapping between the ASV name and database kmers. The file is produced by ``align-regional-kmers`` (we just produced 3!) and is of semantic type ``FeatureData[KmerAlignment]``
-  **frequency-table** – the path to the regional count table produced in denoising. This has a semantic type ``FeatureTable[Frequency]``

The manifest is the only file that is absloutely required to perform
reconstruction.

**Note**

   *The manifest format is specific to the current version. This will be
   deepreicated shortly. However, a version that at least works is
   probably better than a perfct implementation. So… stay tuned?*

Let's look at an example::
	
	id			region-order kmer-map										alignment-map							frequency-table
	WonderWoman		1		 database/sidle-db-wonder-woman-100nt-map.qza	alignment/wonder-woman-align-map.qza	data/wonder-woman-100nt-table.qza
	GreenLantern	3		 database/sidle-db-green-lantern-100nt-map.qza	alignment/green-lantern-align-map.qza	data/green-lantern-100nt-table.qza
	Batman			2		 database/sidle-db-batman-100nt-map.qza			alignment/batman-align-map.qza			data/batman-100nt-table.qza


Parameters
++++++++++

The ``max-mismatch`` and ``per-nucleotide-error`` are used to estimate the probability that a sequence that from the reference is actually a sequencing error or belongs to that sequence. The ``max-mismatch`` value used in reconstruction should match the alignment; by default this is 2 but you may choose to change it in alignmnent with your sequencing length. The authors of the method claim the error rate doesn’t matter; we refer interested reader to original paper’s supplemental material.

The ``min-abundance`` determines the relative abundance of a database sequence to be excluded during optimization.

Now, let’s reconstruct the table, using the default settings.

.. code-block:: shell
	
    qiime sidle reconstruct-counts \
     --m-manifest-file manifest.txt \
     --o-reconstructed-table reconstruction/league_table.qza \
     --o-reconstruction-summary reconstruction/league_summary.qza \
     --o-reconstruction-map reconstruction/league_map.qza

The command will produce a count table, a file containing details about the number of database kmers mapped to a region along with the ASV IDs, and a mapping that’s needed if you want to do taxonomic reconstruction.

Let’s take a look at the count table.

.. code-block:: shell
	
    qiime feature-table summarize \
     --i-table reconstruction/league_table.qza \
     --o-visualization reconstruction/league_table.qzv


You’ll notice that some of the feature IDs contain a ``|`` character, for example, ``1764594|195532|4471854``. This means the two databases sequences could not be resolved during the reconstruction, and so we assign the sequence to both regions. The more regions that are used in the reconstruction, the more likely you are to be able to accurately reconstuct the database sequences.

The second output is a summary. The summary can be used to evaluate the quality of the reconstruction; see the `original manuscript`_ [1]_ for more details. By default, the summary will consider degenerate kmers as unique sequences; you can change the behavior using the ``count-degenerates`` parameter; when False, kmers will only be counted if they belpng to unique reference sequences. You can view the summary by tabulating the metadata.

.. code:: bash

    qiime metadata tabulate \
     --m-input-file reconstruction/league_summary.qza \
     --o-visualization reconstruction/league_summary.qzv


Let’s look at the information for the unresolve feature, ``1764594|195532|4471854``. How many regions is it found it?

Taxonomic Reconstruction
------------------------

Now you have the table reconstructed, you’re ready to reconstruct the taxonomy to match. Specifcially, this process addresses cases where multiple database sequences cannot be untangled. The function takes the database map generated during reconstruction and the taxonomy associated with the database, and returns the reconstructed taxonomy.

There are three possible general cases for a set of shared sequences. First, they can share the full taxonomic string, second they may differ at some point, or third, they may be same until one is missing an assignment. Let’s start with the simpliest case. If we have two database sequences::

   1234    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum
   1235    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum

Then, when we reconstruct taxonomy, everything is the same and the final taxonomic label should be::

   1234 | 1235 k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum

There’s also thee possibility that sequences differ at some higher level, for example::

   1236    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia; s__obeum
   1237    k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Roseburia; s__

In that case, the algorithm would keep the taxonomic assignment associated with the most recent common ancestor::

   1236 | 1237 k__Bacteria; p__Firmictues; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Blautia | g__Roseburia; g__Blautia | g__Rosburia

If the ``--database`` parameter allows the user to select the type of database being used (``greengenes``, ``silva`` or ``none``). If the database is a defined datavase(``greengenes`` or ``silva``), some ad-hoc database cleaning will be performed automatically ✨. For example, if a defined string is::

   k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Entrobacteriales; f__Enterobacteriaceae; g__; s__

Then, the new, cleaned string will be::

    k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Entrobacteriales; f__Enterobacteriaceae; g__unsp. f. Enterobacteriaceae; s__unsp. f. Enterobacteriaceae

The ``--database`` parameter allows the user to select the type of database being used (``greengenes``, ``silva`` or ``none``). If the database is a defined datavase(``greengenes`` or ``silva``), some ad-hoc database cleaning will be performed and uncultured sequences will be handled, specifically with regard to the ``define-missing`` and ``ambiguity-handling`` parameters.

Our database is a subset of the greengenes database, so let’s specify that we used the greengenes database and inheriet the missing strings.

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

The last step in reconstruction is to reconstruct fragments for the phylogenetic tree. Unfortunately, if the reference sequences cannot be resolved, the phylogenetic tree cannot simply be inherieted from the database. So, we need to reconstruct a new phylognetic tree. We handle sequences in two ways.

1. Any database sequence which could full resolved can keep it’s position in the reference tree
2. Sequences which can’t be resolved need to handled somehow.

We could randomly select a sequence to map the reconstructed region to. However, that might not work when there are several sequences that got combine. So, instead, if we can’t resolve the database sequence, we calculate a concensus sequence from the combined data, extract them over the regions we were able to map, and then those concensus sequences can be inserted into a phylogenetic reference backbone using SEPP or something similar.

.. Note::

	Sucessful reconstruction requires that the ids in the database you used as your reference for reconstruction and the database you’re using for alignment are the same. Make sure that you are using the same database release version and the same level of sequence identity

So, our first step is to reconstruct the concensus fragments from sequences that could not be resolved.

.. code-block:: shell

    qiime sidle reconstruct-fragment-rep-seqs \
     --i-reconstruction-map reconstruction/league_map.qza \
     --i-reconstruction-summary reconstruction/league_summary.qza \
     --i-aligned-sequences database/sidle-db-aligned-sequences.qza \
     --m-manifest-file manifest.txt \
     --o-representative-fragments reconstruction/league-rep-seq-fragments.qza

We can then insert the sequences into the reference tree.

.. code-block:: shell

    qiime fragment-insertion sepp \
     --i-representative-sequences reconstruction/league-rep-seq-fragments.qza \
     --i-reference-database ../../../medda-bench/simulations/refs/greengenes/sepp-refs-gg-13-8.qza \
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

* Make sure your :ref:`input manifest <Table Reconstruction>` conforms to the guidelines 
* Your region names must  match between the alignment, kmer, and manifest
* ``count-degenerates`` will control how the summary describes differences in the sequences
* ``max-mismatch`` helps determine the probability sequences should be retained. This should match what was passed to the alignment.
* **NOTE**: THIS WILL CHANGE IN THE NEAR FUTURE. DON'T LET PERFECT BE THE ENEMY OF GOOD ENOUGH

**Syntax**

.. code-block:: bash

	qiime sidle reconstruct-counts \
	 --m-manifest-file [manifest file] \
	 --o-reconstructed-table [reconstructed table] \
	 --o-reconstruction-summary [reconstruction summary] \
	 --o-reconstruction-map [reconstruction map]

**Example**

.. code-block:: bash

	qiime sidle reconstruct-counts \
	 --m-manifest-file region-manifest.tsv \
	 --o-reconstructed-table league_table.qza \
	 --o-reconstruction-summary league_summary.qza \
	 --o-reconstruction-map league_map.qza

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

* A phylogenetic tree can be reconstructed by first, estimating the concensus fragments for the original sequences and then inserting them into a tree.
* See the `q2-fragment-insertion`_ documentation for more inforation

**Fragment reconstruction syntax**

..code-block:: shell
	
	qiime sidle reconstruct-fragment-rep-seqs \
	 --i-reconstruction-map [reconstruction map] \
	 --i-reconstruction-summary [reconstruction summary] \
	 --i-aligned-sequences [aligned sequences] \
	 --m-manifest-file [manifest] \
	 --o-representative-fragments [concensus fragments]

**Example reconstruction syntax**

..code-block:: shell
	
	qiime sidle reconstruct-fragment-rep-seqs \
	 --i-reconstruction-map reconstruction/league_map.qza \
	 --i-reconstruction-summary reconstruction/league_summary.qza \
	 --i-aligned-sequences database/sidle-db-aligned-sequences.qza \
	 --m-manifest-file manifest.txt \
	 --o-representative-fragments reconstruction/league-rep-seq-fragments.qza

References
++++++++++

.. [1] Fuks, C; Elgart, M; Amir, A; et al (2018) "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**:17. doi: 10.1186/s40168-017-0396-x
.. .. .. [2] McDonald, D; Price, NM; Goodrich, J, et al (2012). "An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea." *ISME J*. **6**: 610. doi: 10.1038/ismej.2011.139
.. .. .. [3] Quast, C.; Pruesse, E; Yilmaz, P; et al. (2013) "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." *Nucleic Acids Research*. **41**:D560. doi: 10.1093/nar/gks1219
.. .. .. [4] Martin, M. (2011). "Cutadapt removes adapter sequences from high-throughput sequencing reads". *EMBnet.journal* **17**:10. doi: https://doi.org/10.14806/ej.17.1.200
.. .. .. [5] Callahan, B; McMurdie, P; Rosen, M; et al (2016) "Dada2: High resolution sample inference from Illumina amplicon dada." *Nature Methods*. **13**: 581. doi: https://doi.org/10.1038/nmeth.3869
.. .. .. [6] Amir, A; McDonald, D; Navas-Molina, JA et al. (2017) "Deblur Rapidly Resolves Single-Nucleotide Community Sequence Patterns". *mSystems*. **2**:e00191 doi: 10.1128/mSystems.00191-16
.. .. .. [7] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) "VSEARCH: a versatile open source tool for metagenomics." *PeerJ* 4:e2584 doi: 10.7717/peerj.2584

.. links

.. _here: https://github.com/jwdebelius/q2-sidle/tree/main/docs/tutorial_data
.. _view.qiime2.org: https://view.qiime2.org
.. _absloute paths: https://www.linux.com/training-tutorials/absolute-path-vs-relative-path-linuxunix/
.. _original manuscript: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0396-x
.. _QIIME 2 tutorials: https://docs.qiime2.org/2020.6/tutorials/
.. _q2-fragment-insertion: https://docs.qiime2.org/2020.8/plugins/available/fragment-insertion/

