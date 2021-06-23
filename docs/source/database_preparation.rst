Database Preparation
====================

If you do not already have a database for reconstruction, you will need to prepare your own. (We recommend checking the Resources page to see if you can find one that suits your needs already.) Database preparation will only need to be done once for each database, primer pair and kmer length, since the reference files can  be reused.


You will need:

* A working sidle :doc:`installation <install>`
* A list of the forward and reverse files used to amplify each region of interest
* A reference database of your choice
* Patience to run the extraction


You can start the tutorial by `downloading the database sequences and taxonomy`_. These have already been imported into qiime2 and as Artifacts.

.. code-block:: shell

    mkdir -p sidle_tutorial
    cd sidle_tutorial
    wget https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/database.zip
    unzip database.zip
    cd database

After you run this command, you will find two input files in your folder: ``sidle-db-full-sequences.qza``, which is the full length database sequences and ``sidle-db-taxonomy.qza``, the taxonomy for the sequences.

.. code:: bash

    qiime feature-table tabulate-seqs \
     --i-data sidle-db-full-sequences.qza \
     --o-visualization sidle-db-full-sequences.qzv

You can view the artifact using `qiime2 view`_. You should find that there are 5649 sequences when you check  the data.

Filtering the Database
----------------------

Database preparation can optionally begin by filtering the database to remove sequences with too many degenerate nucleotides or with taxonomic assignments that will not be used.

Degenerate filtering limits memory consumption throughout. The authors of SMURF [1]_ recommend filtering the database to remove sequences with more than 3 degenerate nucleotides. This represents about X% of the greengenes 13_8 database at 99% [2]_ specificity.  a(The the RESCRIPt formatted Silva 138 database is filtered to exclude sequences with more than 5 degenerates [3]_, [4]_). Increasing the number of allowed degenerates (the ``--p-max-degen`` parameter) will allow more sequences through the filter, and may mean more matches in downstream alignment. However, this comes at a substantial increase in the run time and memory needed, since degenerate sequences have to be expanded, meaning more alignments are required.

For this tutorial, we'll start by filtering to remove anything with more than 3 degenerate nucleotides, since this was the recommended threshold in the original algorithm.

.. code:: bash

    qiime sidle filter-degenerate-sequences \
     --i-sequences sidle-db-full-sequences.qza \
     --p-max-degen 3 \
     --o-filtered-sequences sidle-db-full-degen-filtered-sequences.qza


Try summarizing your database again. There should be 5440 sequences remaining. How many do you have?

Some users may also want to filter out sequences which may not be relevant to their analysis, for example, mitochondria or chloroplasts or sequences which are undefined at a high taxonomic level. (Phylum or class, for example.) You can learn more about `filtering by taxonomy`_ in the QIIME2 tutorial, but as a brief example, we'll show filtering a greengenes database for features missing a phylum (**p__;**) or kingdom(**k__;**) designation.

.. code:: bash

    qiime taxa filter-seqs \
     --i-sequences sidle-db-full-degen-filtered-sequences.qza \
     --i-taxonomy sidle-db-taxonomy.qza \
     --p-exclude "p__;,k__;" \
     --p-mode contains \
     --o-filtered-sequences sidle-db-full-degen-filtered-phylum-def-sequences.qza

If you summarize the sequences, you should find that you now have 5408 sequences remaining.

For databases with more complicated strings that include taxonomy, it will be necessary to include the level designation to avoid removing taxa which may be undefined at lower levels.

.. Note::

    The taxonomic filtering should be considered carefully and pre-filtering should be very permissive. Many common databases lack clear taxonomic resolution at lower taxonomic levels (family, genus, species) and these sequences still provide meaningful information in reconstruction.

Once you have finished pre-filtering, you are ready to start extracting regions.


Prepare a regional database for each primer set
-----------------------------------------------

The next step is to extract a region of the database. Alignment with the SMURF algorithm relies on extracting the exact kmer to be aligned with your ASVs, so the primer pair and read length must match exactly. Unlike other techniques, there is, unfortunately, no "good enough" approach. To maximize memory efficiency, the database is also prepared by expanding degenerate nucleotides and collapsing duplicated kmers into a single sequence.

First, the region is extracted from the pre-filtered database using the ``extract-reads`` function from the `feature classifier`_ plugin. As an example, we’ll look at extracting a region between 316F and 484R using the second primer pair from the SMURF paper (5’-``TCCTACGGGAGGCAGCAG``-3’) and (5’-``TATTACCGCGGCTGCTGG``-3’).

.. code:: bash

    qiime feature-classifier extract-reads \
     --i-sequences sidle-db-full-degen-filtered-phylum-def-sequences.qza \
     --p-f-primer TCCTACGGGAGGCAGCAG \
     --p-r-primer TATTACCGCGGCTGCTGG \
     --o-reads sidle-db-filt-316F-484R.qza

For this example, we used the default settings, although these are slightly different from the original SMURF algorithm: In QIIME, the primers are extracted if they have at least an 80% match with the sequence by default; the Matlab implementation of SMURF used a maximum difference of 2 nucleotides [1]_. If you wish to use a limit closer to the original algorithm, this can be changed using the ``--p-identity`` parameter; however, for the sake of this tutorial, we'll use the defaults.

Once the reads have been extracted, they need to be prepared to be used in alignment. This step will expand any degenerate reads that have been extracted, collapse duplicate reads, and trim them to a consistent length. For the full pipeline to work correctly, the primers need to be specified in this step, so once again, you'll need  to pass your primers. You'll also need to specify a trim length; let's use 100nt. Finally, we need to specify a regional identifier in the database using the ``--region`` parameter. This should be the same regional parameter that you use during alignment. We'll call it "WonderWoman" because (a) Diana Prince is amazing and (b) the regional name doesn't matter.

.. code:: bash

    qiime sidle prepare-extracted-region \
     --i-sequences sidle-db-filt-316F-484R.qza \
     --p-region "WonderWoman" \
     --p-fwd-primer TCCTACGGGAGGCAGCAG \
     --p-rev-primer TATTACCGCGGCTGCTGG \
     --p-trim-length 100 \
     --o-collapsed-kmers sidle-db-wonder-woman-100nt-kmers.qza \
     --o-kmer-map sidle-db-wonder-woman-100nt-map.qza

The command will output the sequences (``--o-collapsed-kmers``) with degenerate sequences expanded and duplicated sequences removed and a mapping between the original sequence name and the kmer name (``--o-kmer-map``). You can use qiime to visualize your kmer map, which gives you the relationship between  the original database sequence name (**db-seq**), an expanded name which accounts for degenerates (**seq-name**), the collapsed regional identifier (**kmer**), the primers (**fwd-primer** and **rev-primer**), the region identifier (**region**), and the  sequence length  (**trim-length**).

.. code:: bash

    qiime metadata tabulate \
     --m-input-file sidle-db-wonder-woman-100nt-map.qza \
     --o-visualization sidle-db-wonder-woman-100nt-map.qzv


In some cases, the reference region and sequence length may not be long enough to cover the full amplicon. If that's the case, you can extract the read starting from the reverse primer by setting the trim length to a negative value. You can even reverse complement the resultant amplicons using the ``--reverse_complement_result`` flag. Let's do an example using the same primer-pair region as before, but call the region "Batman". *Note we've swapped the forward and reverse primer sequences.*

.. code:: bash

    qiime sidle prepare-extracted-region \
     --i-sequences sidle-db-filt-316F-484R.qza \
     --p-region "Batman" \
     --p-fwd-primer TATTACCGCGGCTGCTGG \
     --p-rev-primer TCCTACGGGAGGCAGCAG \
     --p-trim-length -100 \
     --p-reverse-complement-result \
     --o-collapsed-kmers sidle-db-batman-100nt-kmers.qza \
     --o-kmer-map sidle-db-batman-100nt-map.qza

As an exercise, try using the 486-650 primers (3-``CAGCAGCCGCGGTAATAC``-5 forward; 3-``CGCATTTCACCGCTACAC``-5 reverse) to extract and prepare a 100nt region called "GreenLantern" as we outlined above. Use the same naming convention as the other two extracted regions (``sidle-db-green-lantern-100nt-kmers.qza``).

Recap
-----
We've essentially constructed three amplicon regional databases. *Again, Batman is simply a reverse-compliment extraction example of the same region we extracted for WonderWoman. This mimics the case in which we have paired-end reads that we could not merge, so we treat them as separate region data as mentioned above.*

+--------------------+----------------------------+--------------------+--------------------+
| Hero (Region Name) | Region                     |  Forward Primer    | Reverse Primer     |
+====================+============================+====================+====================+
| WonderWoman        | 316F-484R                  | TCCTACGGGAGGCAGCAG | TATTACCGCGGCTGCTGG |
+--------------------+----------------------------+--------------------+--------------------+
| Batman             | 316F-484R (rev compliment) | TATTACCGCGGCTGCTGG | TCCTACGGGAGGCAGCAG |
+--------------------+----------------------------+--------------------+--------------------+
| GreenLantern       | 486F-650R                  | CAGCAGCCGCGGTAATAC | CGCATTTCACCGCTACAC |
+--------------------+----------------------------+--------------------+--------------------+

Now, you have a database that's ready to use for alignment and reconstruction.

TL;DR: Database Preparation
---------------------------

Database Filtering
^^^^^^^^^^^^^^^^^^

* Filtering only needs to be performed once per dataset
* Degenerate filtering speeds up preparation and alignment
* You can exclude sequences during database generation that you don't want included in the final table


Degenerate Filtering
""""""""""""""""""""

**Syntax**

.. code-block:: bash

    qiime sidle filter-degenerate-sequences \
     --i-sequences [unfiltered sequences].qza \
     --p-max-degen [degenerate threshold] \
     --o-filtered-sequences [filtered sequences].qza

**Example**

.. code-block:: bash

    qiime sidle filter-degenerate-sequences \
     --i-sequences sidle-db-full-sequences.qza \
     --p-max-degen 3 \
     --o-filtered-sequences sidle-db-full-degen-filtered-sequences.qza

Taxonomic Filtering
"""""""""""""""""""

Please see the `qiime filtering tutorial`_ for more information.

**Syntax**

.. code-block:: bash

    qiime taxa filter-seqs \
     --i-sequences [unfiltered sequences].qza \
     --i-taxonomy [taxonomic descriptions].qza \
     --p-exclude [criteria to exclude] \
     --p-mode contains \
     --o-filtered-sequences [filtered sequences].qza

**Example**

.. code-block:: bash

    qiime taxa filter-seqs \
     --i-sequences sidle-db-full-degen-filtered-sequences.qza \
     --i-taxonomy ref-taxonomy.qza \
     --p-exclude "p__;,k__;" \
     --p-mode contains \
     --o-filtered-sequences filtered-defined-phylum.qza


Database Region Preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* The primers used to extract regions must be the same as the primers used to amplify your sequences in that region
* The extraction command must be re-run for each primer-pair and database
* Read preparation needs to be re-run for each primer-pair, read length, and database
* A negative trim length to ``qiime sidle prepare-extracted-region`` will trim from the reverse primer (right)


Read Extraction
"""""""""""""""

Please see the `qiime feature classifier`_ documentation for more information.

**Syntax**

.. code-block:: bash

    qiime feature-classifier extract-reads \
     --i-sequences [full length sequences] \
     --p-f-primer [forward primer] \
     --p-r-primer [reverse primer] \
     --o-reads [extracted region]

**Example**

.. code-block:: bash

    qiime feature-classifier extract-reads \
     --i-sequences filtered-defined-phylum.qza \
     --p-f-primer TGGCGGACGGGTGAGTAA \
     --p-r-primer CTGCTGCCTCCCGTAGGA \
     --o-reads filtered-defined-phylum-extract-316F-484R.qza

Regional Database Preparation
"""""""""""""""""""""""""""""

**Syntax**

.. code-block:: bash

    qiime sidle prepare-extracted-region \
     --i-sequences [extracted sequences].qza \
     --p-region [region label] \
     --p-fwd-primer [forward primer for region] \
     --p-rev-primer [reverse primer for region] \
     --p-trim-length [kmer length] \
     --o-collapsed-kmers [kmer sequences].qza \
     --o-kmer-map [kmer to database map].qza

**Example**

For forward reads (trim from the left)

.. code-block:: bash

    qiime sidle prepare-extracted-region \
     --i-sequences filtered-defined-phylum-extract-316F-484R.qza \
     --p-region "WonderWoman" \
     --p-fwd-primer TCCTACGGGAGGCAGCAG \
     --p-rev-primer TATTACCGCGGCTGCTGG \
     --p-trim-length 100 \
     --o-collapsed-kmers sidle-db-wonder-woman-100nt-kmers.qza \
     --o-kmer-map sidle-db-wonder-woman-100nt-map.qza

For reverse reads (trim from the right and in this case, reverse complement). The primers should be flipped (we'll trim from the forward primer)

.. code-block:: bash

    qiime sidle prepare-extracted-region \
     --i-sequences filtered-defined-phylum-extract-316F-484R.qza \
     --p-region "Batman" \
     --p-fwd-primer TATTACCGCGGCTGCTGG \
     --p-rev-primer TCCTACGGGAGGCAGCAG \
     --p-trim-length -100 \
     --p-reverse-complement-result \
     --o-collapsed-kmers sidle-db-batman-100nt-kmers.qzv \
     --o-kmer-map sidle-db-batman-100nt-map.qzv


Database References
+++++++++++++++++++

..  websites
.. _filtering by taxonomy: https://docs.qiime2.org/2021.2/tutorials/filtering/#taxonomy-based-filtering-of-tables-and-sequences
.. _qiime filtering tutorial: https://docs.qiime2.org/2021.2/tutorials/filtering/#taxonomy-based-filtering-of-tables-and-sequences
.. _qiime feature classifier: https://docs.qiime2.org/2021.2/tutorials/feature-classifier/#extract-reference-reads
.. _feature classifier: https://docs.qiime2.org/2021.2/tutorials/feature-classifier/#extract-reference-reads
.. _qiime2 view : https://view.qiime2.org
.. _downloading the database sequences and taxonomy: https://github.com/jwdebelius/q2-sidle/raw/main/docs/tutorial_data/database.zip

.. citations

.. [1] Fuks, C; Elgart, M; Amir, A; et al (2018) "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." *Microbiome*. **6**:17. doi: 10.1186/s40168-017-0396-x
.. [2] McDonald, D; Price, NM; Goodrich, J, et al (2012). "An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea." *ISME J*. **6**: 610. doi: 10.1038/ismej.2011.139
.. [3] Quast, C.; Pruesse, E; Yilmaz, P; et al. (2013) "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." *Nucleic Acids Research*. **41**:D560. doi: 10.1093/nar/gks1219
.. [4] Michael S Robeson II, Devon R O'Rourke, Benjamin D Kaehler, et al. "RESCRIPt: Reproducible sequence taxonomy reference database management for the masses."" bioRxiv 2020.10.05.326504; doi: 10.1101/2020.10.05.326504
