#!/bin/bash

FILE_NAME=plot
#SBATCH -n 1
#SBATCH -t 1:00
#SBATCH -p main
#SBATCH --export=ALL
#SBATCH -o log/$FILE_NAME-%j.out
#SBATCH -e log/$FILE_NAME-%j.err

Rscript $FILE_NAME.r
