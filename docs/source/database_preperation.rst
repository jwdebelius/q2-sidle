Database Preparation
====================

The first step in running Sidle (and in the SMURF algorithm) is the preparation of a reference database for alignment. This only needs to be done once for each database, primer pair and kmer length, since the reference files can  be reused. 

This page will walk through the process of preparing a database and extracting a primer region.

Filtering the Database
----------------------

Database preparation can optionally begin by filtering the database to remove sequences with too many degenerate nucleotides or with taxonomic assignments that will not be used. Because the SMURF algorithm depends on alignment against the database, 


, which slow extraction, alignment, and reconstruction. The authors of SMURF [1]_ recommend filtering the database to remove sequences with more than 3 degenerate nucleotides. This represents about X% of the greengenes 13_8 database at 99% [2]_ specificity about Y% of the Silva 132 database [3]_.

As an example, we will work with the greengenes 13_8 database with 88% identity. This is a relatively small database to allow rapid running. You can download a tarball with the necessary example files and move into that folder.

We'll start by filtering to remove anything with more than 3 degenerate nucleotides using the ``--p-max-degen`` parameter.

.. code-block:: shell

    qiime sidle filter-degenerate-sequences \
     --i-sequences 88_otus.qza \
     --p-max-degen 3 \
     --o-filtered-sequences 88_otus-degenerate-filtered.qza

Some users may also want to filter out sequences which are undefined below a certain level, for example, a sequence which does not have a phylum-level assignment. Database filtering is 

.. code-block:: shell

    qiime sidle filter-unannotated-sequences \
     --i-sequences 88_otus-degenerate-filtered.qza \
     --i-taxonomy 88_otus-taxonomy.qza \
     --p-database greengenes \
     --p-empty-level 1 \
     --o-filtered-sequences 88_otus-degenerate-phylum-filtered.qza \



or those that map to mitochondria. This is also designed to decrease memory requirements and remove sequences which will likely be discarded in the final analysis.

.. code-block:: shell

    qiime taxa filter-seqs \
     --i-sequences 88_otus-degenerate-filtered.qza \
     --i-taxonomy 88_otus-taxonomy.qza \
     .. --p-

This example was adapted from the  `QIIME 2 filtering tutorial`_; more examples can be found on their webpage.

.. Note::
    
    The taxonomic filtering should be considered carefully and pre-filtering should be very permissive. Many common databases lack clear taxonomic resolution at lower taxonomic levels and these sequences still provide meaningful information in reconstruction. 


.. There is an option to do minimal curation and automatically remove any sequence which is undefined at a given taxonomic level, for example, any sequence which does not have a defined phylum assignment. This is designed to decrease memory requirements and remove sequences which may be discarded in final analysis anyway. 

.. This can be specified by passing in the path for the taxonomy file (``--reference-taxa-fp``). Currently, the only two databases (``--database``) which are supported are Greengenes and Silva. If you wish to use a custom database, you can format your taxonomy to match the style of the databases by either leaving missing taxonomic strings blank (i.e. ``p__``) in the greengenes style or the Silva convention by including the word "uncultured" or "uncultivated" in the organism name.

.. The lowest taxonomic level to which a sequence must be defined (``--lowest-empty-level``) is then used to filter the sequences. For example, if ``--lowest-taxonomic-level`` is "class", then in a database with greengenes-style taxonomy the sequence linked with this taxonomy would be retained::

..     k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Entrobacteriales;f__;g__;s__

.. but, a sequence without a class would not:::

..     k__Bacteria;p__Proteobacteria;c__;o__;f__;g__;s__

.. Let's try filtering to remove anything with taxonomy undefined at phylum level.

.. .. code-block:: bash

..     sidle filter-reference-database \
..      --reference-seq-fp 88_otus.fasta \
..      --reference-taxa-fp 88_otus_taxonomy.txt \
..      --filt-seq-fp 88_otus-filtered-phylum.fasta \
..      --database greengenes
..      --lowest-empty-level phylum


.. Prepare a regional database
.. ---------------------------

.. The next step is to extract a region of the database. Alignment with the SMURF algorithm relies on extracting the exact kmer to be aligned with your ASVs, so the primer pair and read length must match exactly. Unlike other techniques, there is, unfortunately, no “good enough '' regional approach. 

.. To maximize memory efficiency, the database is also prepared by expanding degenerate nucleotides and collapsing duplicated kmers into a single sequence. This can be accomplished in a single step using SIdle or in two steps combining SIdle with your favorite *in Silico* PCR approach using a two-step preparation. 


.. Mass Database Preperation
.. +++++++++++++++++++++++++

.. The most effecient way to process a database in Sidle is using the `extract` command. This takes a full length 16s database and a list of primer pairs and extracts and tidies the corresponding primer pairs and an output directory to the extracted regions and regional maps. 

.. This is accomplished by preparing a database map file, which lists the region identifier (``region-id``), the forward primer sequence (``forward-primer``) and reverse primer sequence (``reverse-primer``) in a tab-seperated format. For example, if we were using the first three SMURF primers [1]_, our database map would be ::

..     region-id   forward-primer  reverse-primer
..     74-315  TGGCGGACGGGTGAGTAA  CTGCTGCCTCCCGTAGGA
..     316-484 TCCTACGGGAGGCAGCAG  TATTACCGCGGCTGCTGG
..     486-650 CAGCAGCCGCGGTAATAC  CGCATTTCACCGCTACAC

.. In this example, this file is the `primer_map.tsv`. 
 
.. It should alos be noted that this type of extraction apply consistent parameters across all regions. For example, the sample ``--trim-length`` is applied to all regions, so if you want a unique trim length in each region it is better to consider per region extractions (below).

.. As an example, we'll extract 100nt sequences using the first 3 smurf primers.

.. .. code-block:: bash

..     sidle extract \
..      --full-seq-fp 88_otus-filtered-seqs.fasta \
..      --primer-fp primer_map.tsv \
..      --outdir extracted-100 \
..      --trim-length 100

.. The function outputs to a directory (``--outdir``) which will contain all the per-region extracted files where the database maps are named `[region]-map.tsv` and `[region]-seqs.tsv`. In this example, the output files should be

.. .. code-block:: bash

..     ls extracted
..     75-315-map.tsv   74-315-seqs.fasta   316-484-map.tsv
..     316-484-seqs.fasta   486-650-map.tsv 484-650-seqs.fasta


.. Single Region Extraction
.. ++++++++++++++++++++++++

.. A regional database can be extracted and processed in a single step using ``extract-regional-database``. As an example, we’ll look at extracting a region between 74f and 315R using he first primer pair from the SMURF [1]_ paper (5’-TGGCGGACGGGTGAGTAA-3’) and (5’-CTGCTGCCTCCCGTAGGA-3’). For kmer-based alignment we’ll need a set trim length (``--trim-length``), let’s use 100nt. We’ll also trim the primers (``--trim-primers``) before we trim the sequences.

.. Once the sequences are extracted, we can perform database preperation for reconstruction by piping it to a function which will expand degenerate sequences and then collapse duplicated sequences into a single read for alignment. This can be specified using the ``--db-map`` output parameter. We can optionally specify a regional identifier in the database using the ``--region`` parameter. If nothing is specified, then region will be specified as the forward and reverse primer.

.. .. code-block:: bash

..     sidle extract-regional-database \
..      --full-seq-fp 88_otus-filtered-seqs.fasta \
..      --region-fp 88-otus-74-315-tidy-100nt.fasta \
..      --db-map-fp 88-otus-74-315-map.tsv \
..      --region '74-315' \
..      --fwd-primer TGGCGGACGGGTGAGTAA \
..      --rev-primer CTGCTGCCTCCCGTAGGA \   
..      --trim-length 100


.. Two step database preparation
.. +++++++++++++++++++++++++++++

.. SIdle can also be used in conjunction with other regions that have been extracted using QIIME, ipcress, or other techniques. As an example, this might be a region you've already extracted to train a naive bayesian classifier. All SIdle needs to do with these regions is to expand the degenerate sequences, collapse the duplicates, and map the regions. This can be done using the `collapse-regional-database` command.

.. Here, we need the extracted sequences (``--region-fp``), which were extracted using the SMURF 316F-484R primer pair. The command will output the collapsed sequences (``--kmer-db-fp``) and mapping between the kmer name and the original sequence name (``--db-map-fp``). In this case, we need to specify a regional name (``--region``) because it can't be intuited from the data.

.. The SMURF algorithm relies on sequence alignment using sequences that are a consistent length. So, if the extraction technique didn't lead to a consistent sequence length, that can be specified here with the ``--trim-length`` parameter. Let's trim these to 100nt.

.. .. code-block:: bash

..     sidle collapse-regional-database \
..      --region-fp 88-otus-trimmed-316-484-100nt.fasta \
..      --kmer-db-fp 88-otus-316-484-tidy-100nt.fasta \
..      --kmer-map-fp 88-otus-316-484-map.tsv \
..      --trim-length 100 \
..      --region '316-484'

.. Using either the one or two step database extraction, you should prepare database regions for all the primer pairs and read lengths you wish to use. However, once this primer x read length x database combination has been set up, the regions can be re-used for as many experiments as you wish. 


.. citations

.. [1] Fuks, C; Elgart, M; Amir, A; et al (2018) "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**:17. doi: 10.1186/s40168-017-0396-x
.. [2] McDonald, D; Price, NM; Goodrich, J, et al (2012). "An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea." *ISME J*. **6**: 610. doi: 10.1038/ismej.2011.139
.. [3] Quast, C.; Pruesse, E; Yilmaz, P; et al. (2013) "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." *Nucleic Acids Research*. **41**:D560. doi: 10.1093/nar/gks1219

.. websites

.. _QIIME 2 filtering tutorial: https://docs.qiime2.org/2020.2/tutorials/filtering/#taxonomy-based-filtering-of-tables-and-sequences

