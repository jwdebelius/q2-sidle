SMURF Sample Tutorial
=====================

This is intended to replicate the tutorial from the original SMURF paper using the same tutorial data and a prepared database similar to what was provided for that study.

To start, make sure you have the latest version `qiime installed`_; then install q2-sidle according in the :doc:`install instructions <install>`.


Accessing the data
------------------

For this tutorial, we'll be using data derived from the example presented in the `original smurf repository`_. The data comes from a single sample which was amplified using the `smurf primers` [1]_. We'll focus on the second, third, and fourth regions to simplify reconstruction, although the example can easily be extended to the full sequence. The primers for these three regions are extracted from Table 1 of the original paper:

+--------+-----------+------------------------+-----------+------------------------+-------------+
| Region | Forward   | Forward Primer         | Reverse   | Reverse Primer         | Approximate |
|        | position* | (5'-seq-3')            | Position* | (5'-seq-3')            | length      |
+========+===========+========================+===========+========================+=============+
| 1      | 74        | ``TGGCGGACGGGTGAGTAA`` | 315       | ``CTGCTGCCTCCCGTAGGA`` | 240nt       |
+--------+-----------+------------------------+-----------+------------------------+-------------+
| 2      | 316       | ``TCCTACGGGAGGCAGCAG`` | 484       | ``TATTACCGCGGCTGCTGG`` | 160nt       |
+--------+-----------+------------------------+-----------+------------------------+-------------+
| 3      | 486       | ``CAGCAGCCGCGGTAATAC`` | 650       | ``CGCATTTCACCGCTACAC`` | 160nt       |
+--------+-----------+------------------------+-----------+------------------------+-------------+
| 4      | 752       | ``AGGATTAGATACCCTGGT`` | 911       | ``GAATTAAACCACATGCTC`` | 160nt       |
+--------+-----------+------------------------+-----------+------------------------+-------------+
| 5      | 901       | ``GCACAAGCGGTGGAGCAT`` | 1057      | ``CGCTCGTTGCGGGACTTA`` | 160nt       |
+--------+-----------+------------------------+-----------+------------------------+-------------+
| 6      | 1143      | ``AGGAAGGTGGGGATGACG`` | 1336      | ``CCCGGGAACGTATTCACC`` | 190nt       |
+--------+-----------+------------------------+-----------+------------------------+-------------+

\* Forward and reverse position were calcualted as the average position in the Greengenes 13_5 database sequences.

The data was demultiplexed by sample, and has been loaded into QIIME 2 as an artifact. Because the primer pairs are not fully able to cover the regions in all cases, the reads were imported as paired (for regions 2, 3, 4, and 5) and as forward and reverse reads (for regions 1 and 6). For more about this, please see the qiime 2 documentation for `importing sequence data`_. 

So, we'll start by downloading the sequencing artifact and navigating into that directory.

.. code-block:: bash
	
	wget []
	tar -xzf sidle-tutorial
	cd sidle-tutorial

.. If you check the directory contents, you should find 3 files. 

.. * ``all_regions_fwd.qza``
.. * ``all_regions_rev.qza``
.. * ``manifest.tsv``

.. Additionally, you will need to download the reference database, extracted from greengenes 13_8 at 97% identity. 

.. Through out the course of the tutorial, you will have the option to download the pre-computed additional files. 

.. Preparing the Reads
.. -------------------

.. Although the original SMURF paper relied on a quality filtering protocol, we recommend the use of a well-established existing denoising algorith. For this tutorial, we'll use deblur [2]_; you can learn more about multiple options on the ref:`read preperation <read_preperation>` page.


.. Demultiplexing into regions
.. +++++++++++++++++++++++++++

.. We'll use the cutadapt [3]_ plugin in QIIME 2 to demultiplex into regions and remove the primers in the same step. As an example, we'll work with the first two regions of the example file. Let's start by extracting the first region. For this, we'll use the paired end sequences. We can both remove the sequence primer and get *only* sequences with that primer using the ``--p-discard-untrimmed`` flag. This will help us demultiplex into regions of interest. 

.. The only other parameter we'll change here is to change the ``--p-error-rate``, which we'll increase to 0.15. The original paper allowed for 2 errors in the primer during demultiplexing; on an 18nt primer, and will help retain sequences in reverse reads. You can learn more about the parameters for`q2 cutadapt`_; and `cutadapt`_ by looking at their respective documents.

.. We'll start by extracting the first region as paired end sequeces using the forward (``TGGCGGACGGGTGAGTAA``) and reverse (``CTGCTGCCTCCCGTAGGA``) primers.

.. .. code-block:: bash
	
.. 	qiime cutadapt trim-paired \
.. 	 --i-demultiplexed-sequences all_regions_pair.qza \
.. 	 --p-front-f TGGCGGACGGGTGAGTAA \
.. 	 --p-front-r CTGCTGCCTCCCGTAGGA \
.. 	 --p-discard-untrimmed \
.. 	 --p-error-rate 0.15 \
.. 	 --o-trimmed-sequences trimmed-regions/74_315_pair_demux.qza

.. You can check the number of sequences that have been retained by summarizing the trimmed sequences.

.. .. code-block:: bash

.. 	qiime demux summarize \
.. 	 --i-data trimmed-regions/74_315_pair_demux.qza \
.. 	 --o-visualization trimmed-regions/74_315_pair_demux.qza.qzv

.. If open `qiime2 view`_ (view.qiime2.org) and drop the ``e74_315_pair_demux.qza.qzv`` into the browers window, you should see about 20,300 sequences.

.. For the second region, let's try importing the forward and reverse sequences seperately. This may be appropriate if your reads do not overlap (i.e. 2 x 150 sequences of the V4 16s hypervariable region) or if you only have forward reads, for example, those generated by Ion Torrent. Here, we'll use the ``all_regions_fwd.qza``. The forward primer for this region is ``TCCTACGGGAGGCAGCAG``

.. .. code-block:: bash
	
.. 	qiime cutadapt trim-single \
.. 	 --i-demultiplexed-sequences all_regions_fwd.qza \
.. 	 --p-front TCCTACGGGAGGCAGCAG \
.. 	 --p-discard-untrimmed \
.. 	 --p-error-rate 0.15 \
.. 	 --o-trimmed-sequences trimmed-regions/316_484_fwd_demux.qza

.. You can summarized the forward sequences the same way you did the paired end sequences. If you do, how many sequences do you find?

.. We can do the same thing with the reverse sequences.

.. .. code-block:: bash
	
.. 	qiime cutadapt trim-single \
.. 	 --i-demultiplexed-sequences all_regions_rev.qza \
.. 	 --p-front TATTACCGCGGCTGCTGG \
.. 	 --p-discard-untrimmed \
.. 	 --p-error-rate 0.15 \
.. 	 --o-trimmed-sequences trimmed-regions/316_484_rev_demux.qza

.. You can find the remained of the trimmed regions in the ``trimmed-regions`` folder.

.. Denoising the Sequences
.. +++++++++++++++++++++++

.. The next step is to denoise the sequences. For this example, we'll use Deblur. Deblur has been selected because it was suggested by the original authors of the algorithm as an alternative approach to quality filtering, and because it generates ASVs of a standard length. Since kmer-based alignment relies on this standard length, using deblur makes it easier.

.. For the first region, we'll follow the instructions from the qiime2 `alternative methods of read joining`_ tutorial to prepare the reads. This is the only parameter we'll set. More detailed descriptions of additional parameters can be found in the `deblur documentation`_, `the original manuscript`_ [2]_, and the `q2-deblur`_ documentation; we recommend using the defaults.

.. .. code-block:: bash

.. 	mkdir -p joined-seqs/
.. 	qiime vsearch join-paired \
.. 	 --i-demultiplexed-seqs trimmed-regions/74_315_pair_demux.qza \
..  	 --o-joined-sequences joined-seqs/74_315_pair_joined.qza

..  	 mkdir -p quality-filtered/
..  	 qiime quality-filter q-score-joined \
..  	 --i-demux joined-seqs/74_315_pair_joined.qza \
..  	 --o-filtered-sequences quality-filtered/74_315_pair_qc.qza \
..  	 --o-filter-stats quality-filtered/74_315_pair_stats.qza

.. This database region has been pre-cut to 150nt, so we'll use that as the trim parameter for deblur (``--p-trim-length``). 

.. .. code-block:: bash
	
.. 	mkdir -p denoised
.. 	qiime deblur denoise-16S \
.. 	 --i-demultiplexed-seqs quality-filtered/74_315_pair_qc.qza \
.. 	 --p-trim-length 100 \
.. 	 --o-table denoised/74_315_pair_table.qza \
.. 	 --o-representative-sequences denoised/74_315_rep-seqs.qza \
.. 	 --o-stats denoised/74_315_stats.qza

.. The forward reads only need to be quality filtered then denoised. (You can see an example in the `qiime2 moving pictures tutorial`_). For single end regions, we're using a 100nt region, so we'll trim to that length.

.. .. code-block:: bash
	
.. 	qiime quality-filter q-score \
.. 	 --i-demux trimmed-regions/316_484_fwd_demux.qza \
.. 	 --o-filtered-sequences quality-filtered/316_484_fwd_qc.qza \
..  	 --o-filter-stats quality-filtered/316_484_fwd_stats.qza

..  	qiime deblur denoise-16S
..  	 --i-demultiplexed-seqs quality-filtered/316_484_fwd_stats.qza \
.. 	 --p-trim-length 100 \
.. 	 --o-table denoised/7316_484_fwd_table.qza \
.. 	 --o-representative-sequences denoised/316_484_fwd_rep-seqs.qza \
.. 	 --o-stats denoised/316_484_fwd_stats.qza

.. The remaining regions have been denoised and can be found in the ``denoised ``

.. Regional Alignment
.. ------------------

.. Now that we have denoised regions and a database, we can start to perform the regional alignments that are the basis of SMURF reconstruction. We'll do this with the ``align-regional-kmers`` command. 

.. If you haven't downloaded the database, yet, you should do that.

.. .. code-block:: bash
	
.. 	wget X
.. 	untar -czf database.tgz

.. There are two key properties of the database that need to be matched for reconstruction to work. First, the trim length must be the same to allow alignment. Second, the regions must have the same name for internal book keeping. In this case, we'll align with a 100 nt amplicon from the 74-316 region; the region label in the database is ``74-315-100nt-f``, so we need to use the same label here.

.. So, we'll simply align the sequences for each region. As an example, we'll align only the first region; all the regions are aligned using the same command.

.. .. code-block:: bash
	
.. 	mkdir -p alignment
	
.. 	qiime sidle align-regional-kmers \
.. 	 --i-kmers database/gg_97-74-315-100-kmers-seqs.qza \
.. 	 --i-rep-seq denoised/74_315_rep-seqs.qza \
.. 	 --p-region '74-315-100nt-f' \
.. 	 --o-regional-alignment alignment/74_315_aligned-seqs.qza \
.. 	 --o-discarded-sequences alignment/74-315_discard.qza

.. This command may take a while to run, 

.. .. Note::
.. 	If you're not sure about the trim length or region name of an extracted database, you can always drop the artifact into `qiime2 view`_, click on the "provenance" tab and look at the database preperation tab.


.. .. Reconstruction
.. .. --------------

.. .. Build Your Tree
.. .. ---------------

.. .. Analyze your data!
.. .. ------------------

.. .. Getting the Database
.. .. --------------------

.. .. For this tutorial, we'll use a version of the greengenes 13_8 database [2]_ that has been filtered to remove anything with an undefined kingdom or phylum or which has more than 5 degenerate nucleotides along the full length of the sequence. If you need to learn more this, we recommend the page on :ref: `database preperation`. You can see a list of databases associated with the package on `github page`_. 

.. .. Let's download a database

.. .. .. code-block:: shell

.. .. 	wget []

.. TL;DR : Commands Only
.. ---------------------

.. Data Prepereation
.. +++++++++++++++++

.. Importing Data
.. ^^^^^^^^^^^^^^

.. * Follow the `qiime import instructions`_ to  import the data and demultiplex into samples as appropriate. (Or see :ref:seq_demultiplexing) for more information.
.. * You may want to import  forward, reverse, and paired regions if you have regions that do not overlap.

.. Single End Reads
.. ++++++++++++++++

.. Seperating out a single end reads by region
.. """""""""""""""""""""""""""""""""""""""""""

.. * Foward and reverse reads can be treated in the same way - just keep track of which is which

.. **Syntax**

.. .. code-block:: bash
	
.. 	qiime cutadapt trim-single \
.. 	 --i-demultiplexed-sequences [all single end regions].qza \
.. 	 --p-front [primer for that region] \
.. 	 --p-discard-untrimmed \
.. 	 --p-error-rate [error rate] \
.. 	 --o-trimmed-sequences [regional sequences].qza

.. **Example**

.. .. code-block:: bash
	
.. 	qiime cutadapt trim-single \
.. 	 --i-demultiplexed-sequences all_regions_fwd.qza \
.. 	 --p-front TCCTACGGGAGGCAGCAG \
.. 	 --p-discard-untrimmed \
.. 	 --p-error-rate 0.15 \
.. 	 --o-trimmed-sequences trimmed-regions/316_484_fwd_demux.qza


.. Denoising single end reads
.. """"""""""""""""""""""""""

.. * See the `moving pictures tutorial`_ for all the gory details
.. * Make sure your trim length matches the database trim length for your regions

.. **Example**

.. .. code-block:: bash
	
.. 	qiime quality-filter q-score \
.. 	 --i-demux trimmed-regions/316_484_fwd_demux.qza \
.. 	 --o-filtered-sequences quality-filtered/316_484_fwd_qc.qza \
..  	 --o-filter-stats quality-filtered/316_484_fwd_stats.qza

..  	qiime deblur denoise-16S
..  	 --i-demultiplexed-seqs quality-filtered/316_484_fwd_stats.qza \
.. 	 --p-trim-length 100 \
.. 	 --o-table denoised/7316_484_fwd_table.qza \
.. 	 --o-representative-sequences denoised/316_484_fwd_rep-seqs.qza \
.. 	 --o-stats denoised/316_484_fwd_stats.qza


.. Paired end reads
.. ++++++++++++++++


.. Seperating out a paired end reads by region
.. """""""""""""""""""""""""""""""""""""""""""

.. **Syntax**

.. .. code-block::shell

.. 	qiime cutadapt trim-paired \
.. 	 --i-demultiplexed-sequences [reads].qza \
.. 	 --p-front-f [forward  primer] \
.. 	 --p-front-r [reverse primer] \
.. 	 --p-discard-untrimmed \
.. 	 --p-error-rate 0.15 \
.. 	 --o-trimmed-sequences [regional trimmed sequences].qza

.. **Example**

.. .. code-block::shell
	
.. 	qiime cutadapt trim-paired \
.. 	 --i-demultiplexed-sequences all_regions_pair.qza \
.. 	 --p-front-f TGGCGGACGGGTGAGTAA \
.. 	 --p-front-r CTGCTGCCTCCCGTAGGA \
.. 	 --p-discard-untrimmed \
.. 	 --p-error-rate 0.15 \
.. 	 --o-trimmed-sequences 74_315_pair_demux.qza


.. Denoising paired end reads
.. """"""""""""""""""""""""""

.. * See the `Atacama soil microbiome`_ and `alternative methods of read joining`_ tutorials for more information
.. * Make sure your trim length matches the database trim length for your regions

.. **Example**

.. .. code-block:: bash
	
.. 	qiime quality-filter q-score \
.. 	 --i-demux trimmed-regions/316_484_fwd_demux.qza \
.. 	 --o-filtered-sequences quality-filtered/316_484_fwd_qc.qza \
..  	 --o-filter-stats quality-filtered/316_484_fwd_stats.qza

..  	qiime deblur denoise-16S
..  	 --i-demultiplexed-seqs quality-filtered/316_484_fwd_stats.qza \
.. 	 --p-trim-length 100 \
.. 	 --o-table denoised/7316_484_fwd_table.qza \
.. 	 --o-representative-sequences denoised/316_484_fwd_rep-seqs.qza \
.. 	 --o-stats denoised/316_484_fwd_stats.qza


.. .. 

.. .. websites
.. .. _qiime installed: https://docs.qiime2.org/2020.6/install/
.. .. _original smurf repository: https://github.com/NoamShental/SMURF
.. .. _importing sequence data: https://docs.qiime2.org/2020.2/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq
.. .. _q2 cutadapt: https://docs.qiime2.org/2020.2/plugins/available/cutadapt/ 
.. .. _cutadapt: https://cutadapt.readthedocs.io/en/stable/
.. .. _qiime2 view: https://view.qiime2.org/
.. .. _alternative methods of read joining: https://docs.qiime2.org/2020.6/tutorials/read-joining/
.. .. _deblur documentation: https://github.com/biocore/deblur
.. .. _the original manuscript: https://msystems.asm.org/content/2/2/e00191-16
.. .. _q2-deblur: https://docs.qiime2.org/2020.6/plugins/available/deblur/
.. .. _qiime2 moving pictures tutorial: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/#option-2-deblur
.. .. _qiime import instructions: https://docs.qiime2.org/2020.2/install/
.. .. _Atacama soil microbiome: https://docs.qiime2.org/2020.6/tutorials/atacama-soils/

.. .. references/footnotes
.. .. [1] Fuks, C; Elgart, M; Amir, A; et al (2018) "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**:17. doi: 10.1186/s40168-017-0396-x
.. .. [2] Amir, A; McDonald, D; Navas-Molina, JA et al. (2017) "Deblur Rapidly Resolves Single-Nucleotide Community Sequence Patterns". *mSystems*. **2**:e00191 doi: 10.1128/mSystems.00191-16
.. .. [3]  Martin, M. (2011). "Cutadapt removes adapter sequences from high-throughput sequencing reads". *EMBnet.journal* **17**:10. doi: https://doi.org/10.14806/ej.17.1.200
.. .. .. [4] Callahan, B; McMurdie, P; Rosen, M; et al (2016) "Dada2: High resolution sample inference from Illumina amplicon dada." *Nature Methods*. **13**: 581. doi: https://doi.org/10.1038/nmeth.3869
.. .. .. [5] 

.. .. .. [2] McDonald, D; Price, NM; Goodrich, J, et al (2012). "An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea." *ISME J*. **6**: 610. doi: 10.1038/ismej.2011.139