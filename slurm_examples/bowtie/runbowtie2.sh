#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --job-name=bowtie2_example          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=100                   # Real memory (RAM) required (MB)
#SBATCH --time=00:05:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

 module load bowtie2		 # load the module from the already available list of modules

bowtie2-build chrI.fa chrI.fa    # To run bowtie2, we first need to build a database (also known as an index) for our reference genome. To make the bowtie2 index, we can use the 				   #bowtie2-build command. It is usually easiest to name the index the same as the FASTA file-

 # Now we can align some Illumina reads from yeast to chromosome I from the yeast reference genome.
bowtie2 --no-unal -x chrI.fa -U yeast.fastq -S yeast.sam  
