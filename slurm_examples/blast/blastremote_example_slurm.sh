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


#Lets use the protien example for the remote BLAST search 
#the directory contains the unknown non-redundant protein sequence that we will use as our query
#Specify the -remote option
#Set the database: -db nr for non-redundant protein  sequence database
#Change the output format so that you can tell which species the matching sequence is from: -outfmt "6 std stitle"
blastp -remote -num_alignments 20 -query query_protein1.fasta -db nr -outfmt "6 std stitle" -evalue 1e-10 > output.tx

