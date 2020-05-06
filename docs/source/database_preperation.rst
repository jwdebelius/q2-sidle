Database Preparation
====================

If you do not already have a database for reconstruction, you will need to prepare your own. (We recommend checking the Resources page to see if you can find one that suits your needs already.) Database preperation will only need to be done once his only each database, primer pair and kmer length, since the reference files can  be reused. For database preperation, you will need


You will need:

    * A working `installation`_
    * A list of the forward and reverse files used to amplify each region of interest
    * A reference database of your choice
    * Patience to run the extraction


You can start the tutorial by downloading the database sequences and taxonomy. These have already been imported into qiime2 and as Artifacts.
    
.. code-block:: shell
    mkdir sidle-database-tutorial
    cd sidle-database-tutorial
    wget https://docs.qiime2.org/2020.2/data/tutorials/feature-classifier/85_otus.qza
    wget https://docs.qiime2.org/2020.2/data/tutorials/feature-classifier/ref-taxonomy.qza


Filtering the Database
----------------------

Database preparation can optionally begin by filtering the database to remove sequences with too many degenerate nucleotides or with taxonomic assignments that will not be used. Degenerate filtering limits memory consumption through out. The authors of SMURF [1]_ recommend filtering the database to remove sequences with more than 3 degenerate nucleotides. This represents about X% of the greengenes 13_8 database at 99% [2]_ specificity about Y% of the Silva 132 database [3]_.

We'll start by filtering to remove anything with more than 3 degenerate nucleotides using the ``--p-max-degen`` parameter.

.. code-block:: shell

    qiime sidle filter-degenerate-sequences \
     --i-sequences 85_otus.qza \
     --p-max-degen 3 \
     --o-filtered-sequences 85_otus-filtered.qza

If you check the dataset, you'll find the original database contained 5088 sequences, and the new database contains 4790.

Some users may also want to filter out sequences which have undefined taxonomy below a certain level, for example, a sequence without a phylum level assignment. For greengenes, this can be accomplished using the existing qiime command, ``qiime taxa fitler-sequences`` and using the fact that greengenes uses an empty string for missing assignments.

.. code-block:: shell
    
    qiime taxa filter-seqs \
     --i-sequences 85_otus-filtered.qza \
     --i-taxonomy ref-taxonomy.qza \
     --p-exclude "p__;,k__;" \
     --p-mode contains \
     --o-filtered-sequences 85-otus-filtered-defined-phylum.qza

For databases with more complicated strings that include taxonomy, it will be necessary to the level designation to avoid removing taxa which may be undefined at lower levels. Once you have finished pre-filtering, you are ready to start extracting regions. 

.. Note::
    
    The taxonomic filtering should be considered carefully and pre-filtering should be very permissive. Many common databases lack clear taxonomic resolution at lower taxonomic levels and these sequences still provide meaningful information in reconstruction. 


Prepare a regional database
---------------------------

The next step is to extract a region of the database. Alignment with the SMURF algorithm relies on extracting the exact kmer to be aligned with your ASVs, so the primer pair and read length must match exactly. Unlike other techniques, there is, unfortunately, no "good enough" approach. To maximize memory efficiency, the database is also prepared by expanding degenerate nucleotides and collapsing duplicated kmers into a single sequence.

Single Region Extraction
++++++++++++++++++++++++

A regional database can be extracted and processed in a single step using ``extract-regional-database``. As an example, we’ll look at extracting a region between 74f and 315R using he first primer pair from the SMURF paper (5’-TGGCGGACGGGTGAGTAA-3’) and (5’-CTGCTGCCTCCCGTAGGA-3’). For kmer-based alignment we’ll need a set trim length (``--trim-length``), let’s use 100nt. We’ll also trim the primers (``--trim-primers``) before we trim the sequences. The command will also collapse the sequences in preperation for use with sidle. We can optionally specify a regional identifier in the database using the ``--region`` parameter. If nothing is specified, then region will be specified as the forward and reverse primer.

.. code-block:: bash

    qiime sidle extract-regional-database \
     --i-sequences 85-otus-filtered-defined-phylum.qza \
     --p-region "74-315" \
     --p-fwd-primer TGGCGGACGGGTGAGTAA \
     --p-rev-primer CTGCTGCCTCCCGTAGGA \
     --p-trim-length 100 \
     --o-collapsed-kmers 85_otus-74-315-100nt-kmers.qza \
     --o-kmer-map 85_otus-74-315-100nt-map.qza

The command will output the sequences (``--o-collapsed-kmers``) with degenerate sequences expanded and duplicated sequences removed and a mapping between the original sequence name and the kmer-name (``--o-kmer-map``).

It is finally worth noting that sequences prepared this way cannot be re-used for other purposese (for example, training a classifier).

Two step database preparation
+++++++++++++++++++++++++++++

You can also prepare a database from regional sequences you have already extract. As an example, this might be a region you've already extracted to train a naive bayesian classifier. All q2-sidle needs to do with these regions is to expand the degenerate sequences, collapse the duplicates, and map the regions. 

To do this, we'll start by using q2-feature-classifier to extract sequences from the second region in the SMURF paper.

.. code-block:: bash

    qiime feature-classifier extract-reads \
     --i-sequences 85-otus-filtered-defined-phylum.qza \
     --p-f-primer TCCTACGGGAGGCAGCAG \
     --p-r-primer TATTACCGCGGCTGCTGG \
     --p-min-length 50 \
     --p-max-length 200 \
     --o-reads 85-otus-316-484-extracted.qza


Once we have these sequences extracted, then we'll use the ``prepare-extracted-region`` to expand the degenerates and build the database map.

.. code-block::bash
    
    qiime sidle prepare-extracted-region \
     --i-sequences 85-otus-316-484-extracted.qza \
     --p-region "316-484" \
     --p-trim-length 100 \
     --o-collapsed-kmers 85_otus-316-484-100nt-kmers.qza \
     --o-kmer-map 85_otus-316-484-100nt-map.qza

Now, you have a regional kmer database that is ready for reconstruction.

.. citations

.. [1] Fuks, C; Elgart, M; Amir, A; et al (2018) "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**:17. doi: 10.1186/s40168-017-0396-x
.. [2] McDonald, D; Price, NM; Goodrich, J, et al (2012). "An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea." *ISME J*. **6**: 610. doi: 10.1038/ismej.2011.139
.. [3] Quast, C.; Pruesse, E; Yilmaz, P; et al. (2013) "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." *Nucleic Acids Research*. **41**:D560. doi: 10.1093/nar/gks1219



.. .. Mass Database Preperation
.. .. +++++++++++++++++++++++++

.. .. The most effecient way to process a database in Sidle is using the `extract` command. This takes a full length 16s database and a list of primer pairs and extracts and tidies the corresponding primer pairs and an output directory to the extracted regions and regional maps. 

.. .. This is accomplished by preparing a database map file, which lists the region identifier (``region-id``), the forward primer sequence (``forward-primer``) and reverse primer sequence (``reverse-primer``) in a tab-seperated format. For example, if we were using the first three SMURF primers [1]_, our database map would be ::

.. ..     region-id   forward-primer  reverse-primer
.. ..     74-315  TGGCGGACGGGTGAGTAA  CTGCTGCCTCCCGTAGGA
.. ..     316-484 TCCTACGGGAGGCAGCAG  TATTACCGCGGCTGCTGG
.. ..     486-650 CAGCAGCCGCGGTAATAC  CGCATTTCACCGCTACAC

.. .. In this example, this file is the `primer_map.tsv`. 
 
.. .. It should alos be noted that this type of extraction apply consistent parameters across all regions. For example, the sample ``--trim-length`` is applied to all regions, so if you want a unique trim length in each region it is better to consider per region extractions (below).

.. .. As an example, we'll extract 100nt sequences using the first 3 smurf primers.

.. .. .. code-block:: bash

.. ..     sidle extract \
.. ..      --full-seq-fp 88_otus-filtered-seqs.fasta \
.. ..      --primer-fp primer_map.tsv \
.. ..      --outdir extracted-100 \
.. ..      --trim-length 100

.. .. The function outputs to a directory (``--outdir``) which will contain all the per-region extracted files where the database maps are named `[region]-map.tsv` and `[region]-seqs.tsv`. In this example, the output files should be

.. .. .. code-block:: bash

.. ..     ls extracted
.. ..     75-315-map.tsv   74-315-seqs.fasta   316-484-map.tsv
.. ..     316-484-seqs.fasta   486-650-map.tsv 484-650-seqs.fasta
