#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --job-name=example_blast          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=100                   # Real memory (RAM) required (MB)
#SBATCH --time=00:05:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

module load blast             # load the module from the already available list of modules

#First we will need to create a BLAST database from the E. coli FASTA file. We need to specify whether the database contains nucleotide (nucl) or peptide (prot) sequence and we provide the name of the FASTA file with the -in option. 
makeblastdb -dbtype nucl -in e_coli.fasta

blastn -query e_coli_16s.fasta -db e_coli.fasta -outfmt 6 -evalue 1e-10 > output.txt
#Because we are using a nucleotide query and nucleotide database, we need to use blastn
#We need to specify the query using the -query option
#For the database, we use -db and give the name of the FASTA file that we used to create the DB. The program will automatically find the actual database files.
#We also need to specify the format for the output file. There are many options but the tab-delimited output is the easiest to parse through: -outfmt 6
#Finally, it's often useful to set an e-value threshold: -evalue 1e-10




