#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --no-requeue                 # Do not re-run job  if preempted
#SBATCH --job-name=mnist_broadwell            # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH -c 8   # Number of CPU cores
##SBATCH --gres=gpu:1                 # Number of GPUs
#SBATCH --constraint=broadwell
#SBATCH --mem=6000                  # Real memory (RAM) required (MB)
#SBATCH --time=03:30:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.broadwell_4.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.broadwell_4.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

module purge

module use /projects/community/modulefiles
module load py-data-science-stack/5.1.0-kp807 #loads the community Python module that contains Tensorflow environment

source activate tensorflow-1.8.0 #change to 1.7.0 of you wish to run that version of tensorflow
export PYTHONPATH=$PYTHONPATH:/home/aag193/python_mnist/tf/models_1.8 
#Tells python where to look for, for required packages 
#change path to where you saved packages when you cloned the github repo
# change to models_1.7 if you wish to run that versionmnist_1.8.0.py


srun python /home/aag193/python_mnist/mnist_1.8.0.py # --export_dir '/home/aag193/tmp/' --train_epochs 2


