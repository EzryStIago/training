#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --no-requeue                 # Do not re-run job  if preempted
#SBATCH --job-name=pyspark-example   # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=4            # 4 cores
#SBATCH --mem=6000                   # Real memory (RAM) required (MB)
#SBATCH --time=00:10:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export your current env to the job env

WORKING_DIRECTORY='/home/kp807/pyspark_example' # you need to change this to the path to this example

module use /projects/community/modulefiles      # prerequisite for accessing community software
module load spark/2.3.0-kp807                   # finds spark installed in community software

cd $WORKING_DIRECTORY                           # go to the working directory because the app is written with relative paths to input and output

srun spark-submit --master local[4] run_pyspark.py   # submits the spark job as a slurm job

