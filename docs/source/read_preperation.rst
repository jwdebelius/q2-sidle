
Read Preparation
================

Although the original SMURF paper relied on a quality filtering protocol, we have elected to recommend the use of existing denoising algorithms. In Sidle, we recommend using existing tools to perform pre-processing. Here, we provide example code for processing samples, however, this can be accomplished any number of ways.

.. note::

    The following steps represent a **suggested** pipeline for preparing reads to be used with sidle. Alternatively pipelines are possible and may be more useful for your specific circumstances.

As an initial example of read preparation, your pipeline may include the following steps. To make navigation easy, what stage are you at?

* I don't know. Someone handed me files and told me to analyze them!
	* Contact your sequencing provider or the person the files came from if they were sequenced locally
	* Check the repository where they came from. EBI and SRA typically provide files `demultiplexed by sample`_; Qiita provides `deblured tables`_
* I have multiplexed sequences (multiple samples in the same fastq file)
	* You need to `demultiplex by sample and region`_
	* Then, you want to `denoise the sequences`_
* I have two fastq files per sample
	* Your sequences are already demultiplexed by samples; you need to `demultiplex by region`_
	* Then, you want to `denoise the sequences`_
* I have two fastq files per sample and per region
	* Your sequences are already demultiplexed. Make you check you've `demultiplexed correctly`_
	*  Then, you need to `denoise the sequences`_
* I have an ASV table for each region and a corresponding representative sequence file
	* It was `denoised with dada2`_
	* It was `denoised with deblur`_


Demultiplex your reads
----------------------

Fully Multiplexed sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _demux_sample_and_region: 

The first step of sample processing is demultiplexing your sequences into sample x region pairs. You may have fully multiplexed sequences, in which case, you will need to demultiplex into both reads and regions. If this is the case, you will likely have two or three fastq files, which likely represent forward, reverse, and index reads. Please refer to the QIIME 2 documentation for `demultiplexing EMP sequences`_ and `demultiplexing with cutadapt`_.

Now, you're ready to `denoise the regions`_.

.. _demuxed_by_sample:

Sequences barcoded by sample (mixed regions per sample) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If your samples are demultiplexed to include a single set of files (probably forward and reverse), you will need to demultiplex the files into regions. You can use the cutadapt [1]_ plugin in QIIME 2 to demultiplex into regions and remove the primers in the same step. If you have paired-end reads, you will want to use the ``trim-paired`` command; if you have Ion-Torrent or are using a subset of reads, you might choose to use ``trim-single``. (See the `q2 cutadapt`_ and `cutadapt`_ docs for more details.) Using the ``--p-discard-untrimmed`` flag will remove any sequence which does not have the primer region,  allowing you to find our sequences of interest.

.. code-block::shell

	qiime cutadapt trim-paired \
	 --i-demultiplexed-sequences example_seqs.qza \
	 --p-front-f TCCTACGGGAGGCAGCAG \
	 --p-adapter-r TATTACCGCGGCTGCTGG \
	 --p-error-rate 0.1 \
	 --p-indels \
	 --p-discard-untrimmed \
	 --o-trimmed-sequences example_seqs_316_464.qza

This will give you a table per region. Next, continue on to `denoise the regions`_.

.. _demuxed_by_sample_and_sample_region:

Sequences barcoded by sample and region
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are lucky, or planned ahead, you may have added a unique barcode to each region of a sample. So, if, for example, you sequenced 5 samples with 6 regions, after demultiplexing, you should have 30 demultiplexed files and a map which describes the relationship between the name and region. 

You'll need that later. We suggest a file with at least three columns:  ``sample-id``, ``original-sample-id``, and ``regional-id``. (Note that ``sample-id`` is a QIIME-required header column.) You could potentially combine this with a `manifest`_ or `barcode file`_  that you used to import and/or demultiplex your data. (If you don't if your data is demultiplexed already, check in with your sequencing center! Make friends with your sequencing center. Invite them to drink a hot beverage with you and talk cool research.)

Depending on your region lengths and sequencing length, you may need to map the forward and reverse reads separately. So, you may need final sets of demultiplexed reads which represent the

* paired end reads for regions which can be merged
* forward reads only for regions which cannot be merged
* reverse reads only for regions which cannot be merged

.. note::

	If you have a region which does not overlap, the forward and reverse reads should be processed as separate kmers and then reconstructed, rather than concatenated during denoising.

Next, make sure that you trim your primers. If you plan to denoise all the data together, this should be done sequentially for each primer pair that you plan to use for each block of sequences. (i.e. you should trim the paired reads, forward reads, and reverse reads separately.)

.. code-block::shell

	qiime cutadapt trim-paired \
	 --i-demultiplexed-sequences example_seqs.qza \
	 --p-front-f TCCTACGGGAGGCAGCAG \
	 --p-adapter-r TATTACCGCGGCTGCTGG \
	 --p-error-rate 0.1 \
	 --p-indels \
	 --o-trimmed-sequences example_seqs_trimmed1.qza

This will give you a table where each sample is named whatever you've linked to your barcode. From here, you can either `filter your sequences`_ before denoising, or proceed combining all regions, and then filter the table later. If you denoise with Dada2, you may find better performance if you leave the sequences together; this will not affect denoising with deblur.

.. _denoising:

Denoise reads with your favorite algorithm
------------------------------------------

Reads should be denoised, possibly merged, and *trimmed* to a standard length for kmer-based alignment. Depending on the algorithmic approach, `Dada2`_ [2]_  and `Deblur`_ [3]_ might both be good approaches. (If you're interested in an independent comparison of the two methods outside reconstruction, `Nearing et al`_ provides an independent benchmark [4]_).

However, there are some limitations. Ion Torrent and 454 pyrosequencing results should be `denoised with dada2`_. Dada2 tends to retain more high quality sequences than deblur, but may inflate the number of features. Because of the algorithm, it also has a longer run time. 

Illumina data which has already been joined or quality filtered should be `denoised with deblur`_. It's a faster algorithm and highly parallelizable but it's also more conservative. 


.. _denoise_dada2:

DADA2
^^^^^

There are several helpful tutorials on the QIIME 2 website that describe running `dada2 on forward reads`_ and `dada2 on paired reads`_. Minimal pre-processing should be applied before DADA2: simply demultiplex your data and pass it into the command. 

Once DADA2 has been run, you will need to trim the reads to a consistent length. This can be done using the qiime dada2 parameters during denoising, or with the ``trim-dada2-posthoc`` method in q2-sidle. 

As an example of the command, we can download the feature table and representative sequences from the qiime2 `Moving Pictures Tutorial`_ and then practice. 

.. code-block:: bash

	wget https://docs.qiime2.org/2020.6/data/tutorials/moving-pictures/table-dada2.qza .
	wget https://docs.qiime2.org/2020.6/data/tutorials/moving-pictures/rep-seqs-dada2.qza .

If you look at the sequence summary (`viewable here`_), you'll find the sequences have already been trimmed to 120nt. However, for the alignment we plan to do, it may be useful to trim them to 100nt.

.. code-block:: bash

	qiime sidle trim-dada2-posthoc \
	 --i-table table-dada2.qza \
	 --i-representative-sequences rep-seqs-dada2.qza \
	 --p-trim-length 100 \
	 --o-trimmed-table table-dada2-100nt.qza \
	 --o-trimmed-representative-sequences rep-seq-dada2-100nt.qza

You can check the length by tabulating the sequences.

.. code-block:: bash

	qiime feature-table tabulate-seqs \
	 --i-data rep-seq-dada2-100nt.qza \
	 --o-visualization rep-seq-dada2-100nt.qzv

You should find the sequences all trimmed to 100nt, and ready for alignment.

.. _denoise_deblur:

Deblur
^^^^^^

If you have sequenced using Illumina, Deblur may be easier to use and is recommended by the authors/original developers of SMURF. You can find a tutorial for deblurring `single end reads`_  or `paired end reads`_ on the QIIME webpage. Simply set your Deblur trim length to the final kmer length you'll use and proceed. 


Next Step: Reconstruction!
--------------------------

Now, you're ready to proceed to reconstruction!

TL;DR: Read Preparation
-----------------------

A quick flowchart for figuring out how to demultiplex and pre-process your reads.



Demultiplexing
^^^^^^^^^^^^^^
 
* You need to determine if your reads have already been multiplexed and how, and import/demultiplex accordingly
	* `EMP Demultiplexing`_
	* `Cutadapt Demultiplexing`_
	* `Import already demultiplexed reads into QIIME 2`_

* Samples with mixed regions can be extracted using cutadapt to trim primers and discard untrimmed reads

**Paired End Command**

.. code-block::bash
	
	qiime cutadapt trim-paired \
	 --i-demultiplexed-sequences example_seqs.qza \
	 --p-front-f TCCTACGGGAGGCAGCAG \
	 --p-adapter-r TATTACCGCGGCTGCTGG \
	 --p-error-rate 0.1 \
	 --p-indels \
	 --o-trimmed-sequences example_seqs_trimmed1.qza

**Single End Command**

.. code-block:: bash
	
	qiime cutadapt trim-single \
	 --i-demultiplexed-sequences all_regions_fwd.qza \
	 --p-front TCCTACGGGAGGCAGCAG \
	 --p-discard-untrimmed \
	 --p-error-rate 0.15 \
	 --o-trimmed-sequences trimmed-regions/316_484_fwd_demux.qza

Denoising and Table Preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* See the relevant QIIME 2 tutorials:
	* Dada2
		* Single end reads: `Moving Pictures Option 1`_
		* Paired end reads: `Atacama Soils`_ Tutorial
	* Deblur
		* Single end reads: `Moving Pictures Option 2`_
		* Paired end reads: `Alternative Methods of Read Joining`_ Tutorial

* Make sure to trim your sequences to the same length that was used for your database. You can do this with  the command,

.. code-block:: bash

	qiime sidle trim-dada2-posthoc \
	 --i-table table-dada2.qza \
	 --i-representative-sequences rep-seqs-dada2.qza \
	 --p-trim-length 100 \
	 --o-trimmed-table table-dada2-100nt.qza \
	 --o-trimmed-representative-sequences rep-seq-dada2-100nt.qza

References
----------

.. [1] Martin, M. (2011). "Cutadapt removes adapter sequences from high-throughput sequencing reads". *EMBnet.journal* **17**:10. doi: https://doi.org/10.14806/ej.17.1.200
.. [2] Callahan, B; McMurdie, P; Rosen, M; et al (2016) "Dada2: High resolution sample inference from Illumina amplicon dada." *Nature Methods*. **13**: 581. doi: https://doi.org/10.1038/nmeth.3869
.. [3] Amir, A; McDonald, D; Navas-Molina, JA et al. (2017) "Deblur Rapidly Resolves Single-Nucleotide Community Sequence Patterns". *mSystems*. **2**:e00191 doi: 10.1128/mSystems.00191-16
.. [4] Nearing, J.T.; Douglas, G.M.; Comeau, A.M.; Langille, M.G.I. (2018) "Denoising the Denoisers: an independent evaluation of microbiome sequencing error-correction approaches." *Peer J*. **6**: e5364 doi: 10.7717/peerj.5364
.. .. [4] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) "VSEARCH: a versatile open source tool for metagenomics." *PeerJ* 4:e2584 doi: 10.7717/peerj.2584

.. links
.. _demultiplexed by sample: _demuxed_by_sample
.. _deblured tables: _denoise_deblur
.. _demultiplex by sample and region: demux_sample_and_region
.. _denoise the sequences: _denoising
.. _demultiplex by region: _demuxed_by_sample
.. _demultiplexed correctly: _demuxed_by_sample_and_sample_region
.. _denoised with dada2: _denoise_dada2
.. _denoised with deblur: _denoise_deblur
.. _demultiplexing EMP sequences: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/#demultiplexing-sequences
.. _EMP Demultiplexing: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/#demultiplexing-sequences
.. _demultiplexing with cutadapt: https://forum.qiime2.org/t/demultiplexing-and-trimming-adapters-from-reads-with-q2-cutadapt/2313
.. _Cutadapt Demultiplexing: https://forum.qiime2.org/t/demultiplexing-and-trimming-adapters-from-reads-with-q2-cutadapt/2313
.. _Import already demultiplexed reads into QIIME 2: https://docs.qiime2.org/2020.2/tutorials/importing/#fastq-manifest-formats
.. _manifest: https://docs.qiime2.org/2020.2/tutorials/importing/#fastq-manifest-formats
.. _barcode file: https://forum.qiime2.org/t/demultiplexing-and-trimming-adapters-from-reads-with-q2-cutadapt/2313
.. _q2 cutadapt: https://docs.qiime2.org/2020.6/plugins/available/cutadapt/
.. _cutadapt: https://cutadapt.readthedocs.io/en/stable/
.. _filter your sequences: https://docs.qiime2.org/2020.2/plugins/available/demux/filter-samples/
.. _denoise the regions: _denoising
.. _Dada2: https://docs.qiime2.org/2020.6/plugins/available/dada2/
.. _Deblur: https://docs.qiime2.org/2020.6/plugins/available/deblur/
.. _Nearing et al: https://peerj.com/articles/5364/
.. _dada2 on forward reads: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/#option-1-dada2
.. _dada2 on paired reads: https://docs.qiime2.org/2020.6/tutorials/atacama-soils/#paired-end-read-analysis-commands
.. _Moving Pictures Tutorial: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/
.. _viewable here: https://view.qiime2.org/?src=https%3A%2F%2Fdocs.qiime2.org%2F2020.6%2Fdata%2Ftutorials%2Fmoving-pictures%2Frep-seqs.qzv
.. _single end reads: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/#option-2-deblur
.. _paired end reads: https://docs.qiime2.org/2020.6/tutorials/read-joining/
.. _moving pictures option 1: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/#option-1-dada2
.. _moving pictures option 2: https://docs.qiime2.org/2020.6/tutorials/moving-pictures/#option-2-deblur
.. _Atacama Soils: https://docs.qiime2.org/2020.6/tutorials/atacama-soils/#paired-end-read-analysis-commands
.. _alternative methods of read joining: https://docs.qiime2.org/2020.6/tutorials/read-joining/


