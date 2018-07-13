#!/bin/bash

#SBATCH --partition=gpu             # Partition (job queue)
#SBATCH --no-requeue                 # Do not re-run job  if preempted
#SBATCH --job-name=mnist_p100_gpu1            # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --gres=gpu:1                            # Number of GPUs
#SBATCH --nodelist=gpu006                
#SBATCH --mem=6000                  # Real memory (RAM) required (MB)
#SBATCH --time=03:30:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.p100_gpu1.1.8.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.p100_gpu1.1.8.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env


module purge

module use /projects/community/modulefiles
module load py-data-science-stack/5.1.0-kp807  #loads the community Python module that contains Tensorflow environment

source activate tensorflow-gpu-1.8.0 #activates tensorflow 1.8.0 with GPU 

export PYTHONPATH=$PYTHONPATH:/home/aag193/python_mnist/tf/models #Tells python where to look for, for required packages 
        # change to models_1.7 if you wish to run that version



srun python /home/$USER/python_mnist/mnist_1.8.0.py  --model_dir "/home/aag193/tmp/"
# change to mnist_1.7.0.py if you wish to run that version
#the --model_dir option creates checkpoints in the directory specified, it is imoortant to specify
#because if other people ran the job it may try to access their /tmp directory for which you dont have permission

