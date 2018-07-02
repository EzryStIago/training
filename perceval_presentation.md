
## Slurm commands

Interactive: 

  - salloc : allocate resources  e.g.   
  ```salloc -p main -t 10:00 --mem=20M -N 1 -n 1 -c 1```
  
  - srun :   run a job  e.g.  `srun --pty bash -i`  
  - both at the same time:  
  ```srun -p main -t 10:00 --mem=20M -N 1 -n 1 -c 1 --pty bash -i```
  
Batch:  `sbatch my_slurm_batch_script.sh`
  
Checking jour job status: `squeue -u kp807`
  
## Slurm batch script example

```
#!/bin/bash

##### DO NOT SET ANY VARIBLES DECLARATION BEFORE SBATCH INSTRUCTION !!!!! ######
### Modified by Vlad Kholodovych 01/24/2017 ###

#SBATCH --job-name=kris_hello ## Replace with your jobname 
#SBATCH --partition=main # Partition (job queue) 
#SBATCH -N 1 # Number of nodes 
#SBATCH -n 1 # Number of tasks - if not MPI make it 1
#SBATCH -c 2 # Number of cores per task
#SBATCH --mem=100 # Real memory per node required (MB) 
#SBATCH --time=00:10:00 # Total run time limit (HH:MM:SS) 
#SBATCH --output=slurm.%j.%N.out # STDOUT output file 
#SBATCH --error=slurm.%j.%N.err # STDERR error file 

### Declare job non-rerunable
#SBATCH --no-requeue

module use /projects/community/modulefiles/
module load bowtie/1.2.2-gc563

bowtie --help    #this is not parallel execution! no srun!
# alternatively: 
# srun bowtie --help        # slurm will put it in database if srun used
# time srun bowtie --help   # time your execution
```

## Right-sizing your jobs

  - they may run sooner if they are small
     + think of the job scheduler as trying to fit cubes of different sizes - if a hole opens up, a smaller cube might fit in
  - you are using your fair share needlessly if you ask for too much
     + the priority for you is pushed down if you have asked for (even if not used) a lot of resources lately
  - time a portion of the job so that you can estimate the completion time

## Using N, n, and c options in slurm
  - if in doubt, use -N 1 -n 1  (N=nodes, n=tasks per node)
  - 'N'>1 and 'n'>1 option can be used only with MPI-enabled software or such. Not a lot of software falls into that category!
    + software that does: 
       * MPI, 
       * ipython parallel (uses py4mpi), 
       * GATK 4 (uses Spark as parallelizing engine)
  - a lot of software can parallelize across the cores on a single machine so feel free to use '-c' > 1

## Debugging your jobs
  - include -e (error file) and -o (output file) = name it something that you will recognize
  - scontrol show job 123456 -dd  and what it tells you
    + command, batch, workdir
    + resources used (or asked???)
    + only available for a short time after job ends
  - interactive development i.e. run bash on a compute node:
    + IMPORTANT: use only one core at a time; and `scancel 12345` cancel your job or logout to make sure you have relinquished slurm allocation
    + IMPORTANT: without srun, the help team cannot debug properly!

## Slurm good practice
  - number of jobs not more than 5,000 at a time - queue is 10,000 
  - keep track of how much memory you are using 
      
      `sacct -o MaxRSS,job,jobName,state`
  - estimate the wall time 
  
      `sacct -o Elapsed,Start,End,job,jobName`
  - be aware of resources needed for GPUs - cuda machines need to keep a few CPUs free to use GPUs - leave at least 2 cores per GPU

## Bash commands
  - manipulating strings: `basename` and `namedir`; `cut -d '_' -f2 `
  - gotchas like `$file_i`  vs `${file}_i` ; don't leave spaces
  - list comprehensions `echo file{1..10}`
  - Bash command substitution using backticks \`pwd\`

## Data parallelism

This is if your job is "embarassingly" parallel

  - for loops and creating/submitting lots of jobs at once
  - job arrays
  - group-level software management (creating group-level software modules & best practices for installing your own software).

## Disk utilization

How to see who used how much space 

TODO: Abdullah: turn this into a table - e.g. 
|filesystem|command|explanation|
|---|---|---|
|/scratch|.....| Will show file usage on /scratch/netid |

- Will show file usage on /scratch/netid
``` mmlsquota scratch --block-size=auto   ```

- Will show usage in /home/netid

```mmlsquota home --block-size=auto     ```

- Will show individual usage in the shared folder foran (/projects/foran)

```mmlsquota home:foran --block-size=auto ```

- Will show quota and usage of the whole fileset foran.

```mmlsquota -j foran home  --block-size=auto ```

   + `du -hs /directory/to/query/*`
     

TODO: Abdullah: compare the Rmd with the ipynb and add - I think there is many more bash commands - see if they need to go to cheatsheet. 

