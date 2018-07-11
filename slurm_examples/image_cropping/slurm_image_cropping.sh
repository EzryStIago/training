#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --no-requeue                 # Do not re-run job  if preempted
#SBATCH --job-name=img-crop          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000                   # Real memory (RAM) required (MB)
#SBATCH --time=03:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env
#SBATCH --array=0-5                  # number of array tasks

module purge

module use /projects/community/modulefiles
module load imagemagick/7.0.8-kp807

WORKING_DIRECTORY='/home/kp807/projects/slurm_examples/image_cropping/'

cd $WORKING_DIRECTORY

readarray IMAGE_PATCHES < list_of_coordinates.txt

#COORD is a variable with a value 100x80+1500+1200 for a patch 100x80 with offset 1500 and 1200 for upper left corner
COORD=${IMAGE_PATCHES[$SLURM_ARRAY_TASK_ID]}     #gets the coordinate string
OUTFILE=`echo $COORD | tr [+x] [..]`               #removes + and x sign from COORD and replaces both with .

convert lena.png -crop $COORD results/lena_cropped_$OUTFILE.png

