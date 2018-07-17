# Blast(Basic Local Alignment Search Tool) Example

BLAST finds regions of similarity between biological sequences. The program compares nucleotide or protein	sequences to sequence databases and calculates the statistical significance. Learn more [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
The table summarizes the different types of BLAST searches that we can perform:
  command     |query	  |database
-------|----------|-------------
blastn |nucleotide|nucleotide
blastp |protein	  | protein
blastx |nucleotide|protein
tblastn|protein	  | nucleotide

## Whats happening
  * To perform a local BLAST search we need a BLAST database and a FASTA file containing our query sequence

  * In this example, we will search the 16s ribosomal RNA sequence against the E. coli genome, both of which we have as FASTA files.

**IMPORTANT:** If you have two FASTA files, your query should be the smaller sequence (e.g. a single gene) and you should make the larger sequence (e.g. the E. coli genome) into your database

  * If you type: ``blastn -help`` it will print a list of command-line options, including all possible ways to modify the output format

## Whats contained in the Directory
Instead of performing a local BLAST search you might want to search a sequence against a very large database, containing sequences from many different species, rather than a specific genome.
 The NCBI mantains such a database of updated DNA and protien sequences.
To perform a remote BLAST search, which sends your query to the NCBI servers which then perform the search against the most recent version of the database is pretty starightforward from the command line.

  * The ``blastlocal_example_slurm.sh`` does the local BLAST search using the ``e_coli.fasta`` and the ``e_coli_16s.fasta`` files.
  * The ``blastremote_example_slurm.sh`` uses the ``query_protein1.fasta`` file to do a remote BLAST search against the NCBI database.
  * Both scripts are explained with detailed comments 
