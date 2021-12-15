# Changelog

## 2022.2

* Refactors reconstruction to seperate database reconstruction from count reconstruction. Users can now propegate simple database files instead of repeatedly passing kmer maps. This refactor also significantly improves performance.
* Allows the tracking of counts that have been filtered out of the table by sample
<!-- * Introduces a pipeline to streamline reconstruction (`side_reconstruct`) which will reconstruct the database, table, and taxonomy
* Introduces a pipeline to streamline tree building (`reconstruct_tree`) which will generate representative fragments and insert them into a SEPP tree -->
<!-- * Introduces a function to to track -->