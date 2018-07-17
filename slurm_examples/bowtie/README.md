# Bowtie Example
Bowtie is an ultrafast, memory-efficient short read aligner geared toward quickly aligning large sets of short DNA sequences (reads) to large genomes. Bowtie works best when aligning short reads to large genomes (e.g. human or mouse), though it supports arbitrarily small reference sequences and reads as long as 1024 bases.

This directory contains an example of how to use the bowtie module to efficienlty compare an unknown chromosome with the reference genome from yeast.

## Understanding the Script
  * Bowtie module is provided and is loaded when you run the script

  * Before we can start aligning reads to our reference genome, we need to create an index of the genome. We will use the chromosome I to build our index. To make the bowtie2 index, we can use the bowtie2-build command:
``bowtie2-build FILENAME INDEXNAME &> index.log``

  * It is usually easiest to name the index the same as the FASTA file. **Bowtie2** will then create several index files with the same name as the FASTA file with an extension at the end. Now we are ready to align our Illumina reads to the reference genome.

  * At a minimum, **Bowtie2** requires the index name, the name of the file that contains your Illumina sequence, and the output file. We will only be aligning single-end sequence, so we use the -U option which stands for unpaired.
``bowtie2 -x INDEX -U ILLUMINA -S OUTPUT``
It is also useful to use the --no-unal option, otherwise bowtie reports the read ID and sequence even if it was not able to align it to the genome, which is unnecessary

  * The ``yeast.fastq`` contains the FASTA sequence for a yeast chromosome that is used as the reference and the output is stored in the ``yeast.sam`` file.
