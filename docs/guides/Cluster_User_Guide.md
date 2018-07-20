# General Information

There are several clusters administered by OARC. Each cluster is using the resource scheduler `Slurm` to allocate users the requested resources. 

- **Amarel** - a general access cluster. [Request access](https://oarc.rutgers.edu/access/) or [become an owner](https://oarc.rutgers.edu/services/condo-model/)
- **Perceval** - similar to Amarel, but paid for, and to be used by exclusively by, NIH grants and grantees
- **NM3** - [Newark Massive Memory Machine](http://nm3wiki.newark.rutgers.edu/index.php/Main_Page)
- **Sirius** - a single big machine with no resource scheduler, acts as a big desktop with multiple concurrent users

#Cluster Resources

## Amarel
Amarel is a CentOS 7 Linux compute cluster that is actively growing through the combination of separate computing clusters into a single, shared resource.

Amarel includes the following hardware (this list may already be outdated since the cluster is actively growing):
```
52 CPU-only nodes, each with 28 Xeon e5-2680v4 (Broadwell) cores + 128 GB RAM
20 CPU-only nodes, each with 28 Xeon e5-2680v4 (Broadwell) cores + 256 GB RAM
4 28-core e5-2680v4 nodes each with 2 x Nvidia Pascal P100 GPUs onboard
2 high-memory nodes, each with 56 e7-4830v4 (Broadwell) cores + 1.5 TB RAM
53 CPU-only nodes, each with 16 Intel Xeon e5-2670 (Sandy Bridge) cores + 128 GB RAM
5 CPU-only nodes, each with 20 Intel Xeon e5-2670 (Ivy Bridge) cores + 128 GB RAM
26 CPU-only nodes, each with 24 Intel Xeon e5-2670 (Haswell) cores + 128 GB RAM
4 CPU-only nodes, each with 16 Intel Xeon e5-2680 (Broadwell) cores + 128 GB RAM
3 12-core e5-2670 nodes with 8 Nvidia Tesla M2070 GPUs onboard
2 28-core e5-2680 nodes with 4 Quadro M6000 GPUs onboard
1 16-core e5-2670 node with 8 Xeon Phi 5110P accelerators onboard
```
Default run time = 2 hours in the 'main' partition
Maximum run time = 3 days in the 'main' partition

## Perceval 

Perceval has the same setup as Amarel. The same file system is mounted on both clusters and changes made to the filesystem on one cluster will be reflected in the other cluster as well. 

Perceval includes the following hardware (this list may already be outdated since the cluster is actively growing):

132 CPU-only nodes, each with 24 Intel Xeon E5-2680 cores + 128 GB RAM
8 GPU nodes with 24 Intel Xeon E5-2680 cores + 128 GB RAM
1 Large Memory node with 48 Intel Xeon E5-2680 cores + 1.5 TB RAM

Default run time = 2 hours for all partitions
Maximum run time = 7 days in the 'main' partition, 2 days in the 'gpu' partition

# Basic operations - connecting and moving files

## Connecting to the cluster

- Amarel is currently accessed using a single hostname, `amarel.rutgers.edu`
- Perceval is currently accessed using a single hostname, `perceval.rutgers.edu`

For the example purposes, we will assume you are connecting to Amarel. When you connect to this system, your log-in session (your Linux shell) will begin on one of multiple log-in nodes, named amarel1, amarel2, etc. So, while you are logged-in to Amarel, you will see "amarel1" or "amarel2" as the name of the machine you are using.

`ssh [your NetID]@amarel.rutgers.edu`

If you are connecting from a location outside the Rutgers campus network, you will need to first connect to the campus network using the Rutgers VPN (virtual private network) service. See [here](https://oit.rutgers.edu/vpn) for details.

## Moving files

There are many different ways to this: secure copy (scp), remote sync (rsync), an FTP client (FileZilla), etc. Let’s assume you’re logged-in to a local workstation or laptop (not already logged-in to Amarel). To send files from your local system to your Amarel /home directory,  
```scp file-1.txt file-2.txt [NetID]@amarel.rutgers.edu:/home/[NetID]```

To pull a file from your Amarel /home directory to your laptop (note the “.” at the end of this command),  
```scp [NetID]@amarel.rutgers.edu:/home/[NetID]/file-1.txt  . ```

If you want to copy an entire directory and its contents using scp, you’ll need to “package” your directory into a single, compressed file before moving it:  
```tar -czf my-directory.tar.gz my-directory```

After moving it, you can unpack that .tar.gz file to get your original directory and contents:  
```tar -xzf my-directory.tar.gz```

A handy way to synchronize a local file or entire directory between your local workstation and the Amarel cluster is to use the rsync utility. First, let's sync a local (recently updated) directory with the same directory stored on Amarel:  
```rsync -trlvpz work-dir gc563@amarel.rutgers.edu:/home/gc563/work-dir```

In this example, the rsync options I'm using are t (preserve modification times), r (recursive, sync all subdirectories), l (preserve symbolic links), v (verbose, show all details), p (preserve permissions), z (compress transferred data)

To sync a local directory with updated data from Amarel:  
```rsync -trlvpz gc563@amarel.rutgers.edu:/home/gc563/work-dir work-dir```

Here, we've simply reversed the order of the local and remote locations.

For added security, you can use SSH for the data transfer by adding the e option followed by the protocol name (ssh, in this case):  
```rsync -trlvpze ssh gc563@amarel.rutgers.edu:/home/gc563/work-dir work-dir```

## OnDemand - GUI for the cluster

For users not familiar with Linux, there is an option to connect to the cluster via a web browser (either from campus or through VPN) until you get more confortable with Linux. Both connecting and moving files can be achieved through [this interface](https://ondemand.hpc.rutgers.edu). However, you are strongly encouraged to get comfortable with Linux, as your productivity will soar and a GUI is never as flexible as a command line interface.  

Functionalities of OnDemand: 

 - file upload
 - file editor
 - linux shell 
 - launch Jupyter notebook
 - launch RStudio
 - compose slurm job from a template
 - submit slurm job
 - view job queue 

# Listing available resources

Before requesting resources (compute nodes), it’s helpful to see what resources are available and what cluster partitions (job queues) to use for certain resources.

Example of using the **sinfo** command:
```
sinfo

PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
main*        up 3-00:00:00      4 drain* hal[0050-0051,0055,0093]
main*        up 3-00:00:00      5  down* slepner[084-088]
main*        up 3-00:00:00      4  drain hal[0023,0025-0027]
main*        up 3-00:00:00     86    mix gpu[003-004,006],hal[0001-0008,0017-0018,0022,0024,0028-0032,0044-0047,0054,0056-0057,0062-0068,0073-0079,0081-0092,0094-0096],mem002,pascal[001-006],slepner[010-014,016,018-023,036,042-044,046,048,071,074,076,081-082]
main*        up 3-00:00:00     84  alloc gpu005,hal[0009-0016,0019-0021,0033-0043,0048-0049,0052-0053,0058-0061,0069-0072,0080],mem001,slepner[009,015,017,024-035,037-041,045,047,054-070,072-073,075,077-080,083]
main*        up 3-00:00:00      2   down gpu[001-002]
gpu          up 3-00:00:00      8    mix gpu[003,006],pascal[001-006]
gpu          up 3-00:00:00      1  alloc gpu005
gpu          up 3-00:00:00      2   down gpu[001-002]
phi          up 3-00:00:00      1    mix gpu004
mem          up 3-00:00:00      1    mix mem002
mem          up 3-00:00:00      1  alloc mem001
```
Understanding this output:

- There are 4 basic partitions, main (traditional compute nodes, CPUs only), gpu (nodes with general-purpose GPU accelerators), mem (CPU-only nodes with 1.5 TB RAM), phi (CPU-only nodes with Xeon Phi coprocessors.
- The upper limit for a job’s run time is 3 days (72 hours).

**Term** |**Meaning**
---|---
Allocated (alloc) |nodes that are currently running jobs.
Mixed (mix)| nodes have jobs using some, but not all, CPU cores onboard.
Idle| nodes are currently available for new jobs.
Drained (drain, drng) |nodes are not available for use and may be offline for maintenance.
Slepner, Norse Mythology| "Sleipnir" 8-legged war horse (this made more sense when CPUs had 8 cores).
Hal| Hal is a dependable member of the Discovery One crew who does an excellent job of following instructions.
Pascal| French mathematician and the name of one of NVIDIA's GPU architectures.
CUDA| This is the name of a parallel computing platform and application programming interface (API) model created by Nvidia


# Loading software modules

When you first log-in, only basic system-wide tools are available automatically. To use a specific software package that is already installed, you can setup your environment using the module system.

## Lmod video walkthrough

[![Using module command](https://img.youtube.com/vi/tu9cmUrutgE/0.jpg)](https://www.youtube.com/watch?v=tu9cmUrutgE)

Commands used in the video: 
```
lmod = https://lmod.readthedocs.io = A New Environment Module System
Solves the problem of setting environment variables so you can run different software easily

Main commands: 
module avail                                  # which modules are available
module spider                                 # find information about a software
module load                                   # load a particular software
module use  /projects/community/modulefiles   # user contributed software, unsupported
module purge                                  # unload all modules
module list                                   # list currently loaded modules

gotcha: will not show modules until you load up their dependencies - e.g. R-Project

module avail               # available system-installed modules
which R                    # no results
module spider R            # too many results
module spider R-           # keyword search for specific module
module load intel/17.0.4   # after load, we will see R-Project in the list of modules
module avail               
module load R-Project/3.4.1# loading R
module list
which R                    # now R is in my path

srun -N 1 -n 1 -c 1 -t 10:00 R --no-save

module use  /projects/community/modulefiles   # user contributed software, unsupported
module avail                                  # which modules are available
module load spark/2.3.0-kp807                 # load a particular software
which spark-shell                             # shows where the software is located
```

## Longer explanation

The module avail command will show a list of the core (primary) modules available:
```
module avail

-------------------------------------------------- /opt/sw/modulefiles/Core --------------------------------------------------
   ARACNE/20110228         blat/35                  gcc/4.9.3               intel_mkl/16.0.3 (D)    mvapich2/2.2        (D)
   HISAT2/2.0.4            bowtie2/2.2.6            gcc/4.9.4               intel_mkl/17.0.0        openmpi/2.1.1
   HISAT2/2.1.0     (D)    bowtie2/2.2.9     (D)    gcc/5.3                 intel_mkl/17.0.1        pgi/16.9
   MATLAB/R2017a           bwa/0.7.12               gcc/5.4          (D)    intel_mkl/17.0.2        pgi/16.10           (D)
   MATLAB/R2017b    (D)    bwa/0.7.13        (D)    hdf5/1.8.16             java/1.7.0_79           python/2.7.11
   Mathematica/11.1        cuda/7.5                 intel/16.0.1            java/1.8.0_66           python/2.7.12
   OpenCV/2.3.1            cuda/8.0                 intel/16.0.3     (D)    java/1.8.0_73           python/3.5.0
   STAR/2.5.2a             cuda/9.0          (D)    intel/16.0.4            java/1.8.0_121          python/3.5.2        (D)
   Trinotate/2.0.2         cudnn/7.0.3              intel/17.0.0            java/1.8.0_141          samtools/0.1.19
   bamtools/2.4.0          cufflinks/2.2.1          intel/17.0.1            java/1.8.0_152   (D)    samtools/1.2
   bcftools/1.2            delly/0.7.6              intel/17.0.2            modeller/9.16           samtools/1.3.1      (D)
   bedtools2/2.25.0        gaussian/03revE01        intel/17.0.4            moe/2016.0802           trinityrnaseq/2.1.1
   blast/2.6.0             gaussian/09revD01 (D)    intel_mkl/16.0.1        mvapich2/2.1
```

Understanding this output:

The packages with a (D) are the default versions for packages where multiple versions are available.

- To see a comprehensive list of all available modules (not just the core modules) use the `module spider` command, or `module keyword` command (see CheatSheat section). 
- To be able to use community-contributed software, execute `module use /projects/community/modulefiles` before using the the above mentioned commands (see Community section). 
```
module spider

---------------------------------------------------------------------------------------------------------------
The following is a list of the modules currently available:
---------------------------------------------------------------------------------------------------------------
  ARACNE: ARACNE/20110228
    ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context

  HISAT2: HISAT2/2.0.4, HISAT2/2.1.0
    HISAT2: graph-based alignment of next generation sequencing reads to a population of genomes

  HMMER: HMMER/3.1b2
    HMMER: biosequence analysis using profile hidden Markov models

  MATLAB: MATLAB/R2017a, MATLAB/R2017b
    MATLAB: The Language of Technical Computing

  Mathematica: Mathematica/11.1
    Wolfram Mathematica: Modern Technical Computing

  NAMD: NAMD/2.10
    NAMD: Scalable Molecular Dynamics

  ORCA: ORCA/3.0.3
    ORCA: An ab initio, DFT and semiempirical SCF-MO package

  OpenCV: OpenCV/2.3.1
    OpenCV: Open Source Computer Vision
```

Loading a software module changes your environment settings so that the executable binaries, needed libraries, etc. are available for use.   

|command|explanation|usage|
|-----|--------|----------|
|module load|load a software module| `module load intel/16.0.3`|
|module unload|remove select modules|`module unload intel/16.0.3`|
|module purge| load the default version of any software package| `module load intel`|

Below are some examples.
```
module load intel/16.0.3 mvapich2/2.1
module list
Currently Loaded Modules:
  1) intel/16.0.3   2) mvapich2/2.1

module unload mvapich2/2.1
module list
Currently Loaded Modules:
  1) intel/16.0.3

module purge
module list
No modules loaded

module load intel
module list
Currently Loaded Modules:
  1) intel/16.0.3
```

If you always use the same software modules, your `~/.bashrc` (a hidden login script located in your /home directory) can be configured to load those modules automatically every time you log in. Just add your desired module load command(s) to the end of that file. You can always edit your `~/.bashrc` file to change or remove those commands later.

**PLEASE NOTE:** Software installed cluster-wide is typically configured with default or standard (basic) options, so special performance-enhancing features may not be enabled. This is because the Amarel cluster comprises a variety of hardware platforms and cluster-wide software installations must be compatible with all of the available hardware (including the older compute nodes). If the performance of the software you use for your research can be enhanced using hardware-specific options (targeting special CPU core instruction sets), you should consider installing your own customized version of that software in your /home directory.

# Running slurm jobs

## Video walkthrough

Demoing `sinfo, srun, squeue, scancel` commands: 

[![Slurm intro](https://img.youtube.com/vi/RY7yeW3X2-U/0.jpg)](https://www.youtube.com/watch?v=RY7yeW3X2-U)


## Running a serial (single-core) job

Here’s an example of a SLURM job script for a serial job. I’m running a program called “zipper” which is in my /scratch (temporary work) directory. I plan to run my entire job from within my /scratch directory because that offers the best filesystem I/O performance.
```
#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=zipx001a          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=2000                   # Real memory (RAM) required (MB)
#SBATCH --time=02:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

cd /scratch/[your NetID]

module purge
module load intel/16.0.3 fftw/3.3.1

srun /scratch/[your NetID]/zipper/2.4.1/bin/zipper < my-input-file.in
```

Understanding this job script:

- A job script contains the instructions for the SLURM workload manager (cluster job scheduler) to manage resource allocation, scheduling, and execution of your job.
- The lines beginning with #SBATCH contain commands intended only for the workload manager.
- My job will be assigned to the “main” partition (job queue).
- If this job is preempted, it will be returned to the job queue and will start again when required resources are available
- This job will only use 1 CPU core and should not require much memory, so I have requested only 2 GB of RAM — it’s a good practice to request only about 2 GB per core for any job unless you know that your job will require more than that.
- My job will be terminated when the run time limit has been reached, even if the program I’m running is not finished. It is not possible to extend this time after a job starts running.
- Any output that would normally go to the command line will be redirected into the output file I have specified, and that file will be named using the compute node name and the job ID number.

Be sure to configure your environment as needed for running your application/executable. This usually means loading any needed modules before the step where you run your application/executable.


Here’s how to run a serial batch job, loading modules and using the **sbatch** command:  
```sbatch my-job-script.sh```  

The **sbatch** command reads the contents of your job script and forwards those instructions to the SLURM workload manager. Depending on the level of activity on the cluster, your job may wait in the job queue for minutes or hours before it begins running.


## Running a parallel (multicore MPI) job

Here’s an example of a SLURM job script for a parallel job. See the previous (serial) example for some important details omitted here.
```
#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=zipx001a          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=16                  # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=124000                 # Real memory (RAM) required (MB)
#SBATCH --time=02:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

cd /scratch/[your NetID]

module purge
module load intel/16.0.3 fftw/3.3.1 mvapich2/2.1

srun --mpi=pmi2 /scratch/[your NetID]/zipper/2.4.1/bin/zipper < my-input-file.in
```
Understanding this job script:

- The srun command is used to coordinate communication among the parallel tasks of your job. You must specify how many tasks you will be using, and this number usually matches the –ntasks value in your job’s hardware allocation request.
- This job will use 16 CPU cores and nearly 8 GB of RAM per core, so I have requested a total of 124 GB of RAM — it’s a good practice to request only about 2 GB per core for any job unless you know that your job will require more than that.
- Note here that I’m also loading the module for the parallel communication libraries (MPI libraries) needed by my parallel executable.


Here’s how to run a parallel batch job, loading modules and using the **sbatch** command:  
```sbatch my-job-script.sh```

Note here that I’m also loading the module for the parallel communication libraries (MPI libraries) needed by my parallel executable.


Here’s how to run a parallel batch job, loading modules and using the **sbatch** command:  
```sbatch my-job-script.sh```

## Running array of jobs

Array job is an approach to handle multiple jobs with single job script. Here is an example to submit 500 jobs with single job script. 
```
#!/bin/bash
#SBATCH --partition=main             # Name of the partition
#SBATCH --job-name=arrayjobs         # Name of the job
#SBATCH --ntasks=1                   # Number of tasks
#SBATCH --cpus-per-task=1            # Number of CPUs per task
#SBATCH --mem=1GB                    # Requested memory
#SBATCH --array=0-499                # Array job will submit 500 jobs
#SBATCH --time=00:10:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT file
#SBATCH --error=slurm.%N.%j.err      # STDERR file 

echo  -n "Executing on the machine: " 
hostname
echo "Array Task ID : " $SLURM_ARRAY_TASK_ID 
echo " Random number : " $RANDOM
```
In the above description, the line 
``` #SBATCH --array=0-499 ```
submits 500 jobs. 

The "%" seperator is useful to limit the number of jobs in the queue at a any given time. For example, the following line would send a maximum of 100 jobs to the queue. 
```
#SBATCH --array=0-499%100
```

It is a good practice to populate the queue with less than 500 jobs at any given time. 


## Running an interactive job

An interactive job gives you an active connection to a compute node (or collection of compute nodes) where you will have a login shell and you can run commands directly on the command line. This can be useful for testing, short analysis tasks, computational steering, or for running GUI-based applications.

When submitting an interactive job, you can request resources (single or multiple cores, memory, GPU nodes, etc.) just like you would in a batch job:
```
[NetID@amarel1 ~]$ srun --partition=main --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2000 --time=00:30:00 --export=ALL --pty bash -i

srun: job 1365471 queued and waiting for resources
srun: job 1365471 has been allocated resources

[NetID@slepner045 ~]$
```
Notice that, when the interactive job is ready, the command prompt changes from NetID@amarel1 to NetID@slepner045. This change shows that I’ve been automatically logged-in to slepner045 and I’m now ready to run commands there. To exit this shell and return to the shell running on the amarel1 login node, type the exit command.

# Monitoring the status of jobs
_
The simplest way to quickly check on the status of active jobs is by using the **squeue** command:
```
squeue -u [your NetID]

  JOBID PARTITION     NAME     USER  ST       TIME  NODES NODELIST(REASON)
1633383      main   zipper    xx345   R       1:15      1 slepner36
```
Here, the state of each job is typically listed as being either PD (pending), R (running), along with the amount of allocated time that has been used (DD-HH:MM:SS).

For summary accounting information (including jobs that have already completed), you can use the **sacct** command:
```
sacct

       JobID    JobName  Partition    Account  AllocCPUS      State ExitCode 
------------ ---------- ---------- ---------- ---------- ---------- -------- 
1633383          zipper       main      statx         16    RUNNING      0:0
```
Here, the state of each job is listed as being either PENDING, RUNNING, COMPLETED, or FAILED.

For complete and detailed job info, you can use the **scontrol show job [JobID]** command:
```
scontrol show job 244348

JobId=244348 JobName=XIoT22
   UserId=gc563(148267) GroupId=gc563(148267) MCS_label=N/A
   Priority=5050 Nice=0 Account=oarc QOS=normal
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=1-04:07:40 TimeLimit=2-00:00:00 TimeMin=N/A
   SubmitTime=2017-05-14T07:47:19 EligibleTime=2017-05-14T07:47:19
   StartTime=2017-05-14T07:47:21 EndTime=2017-05-16T07:47:21 Deadline=N/A
   PreemptTime=None SuspendTime=None SecsPreSuspend=0
   Partition=main AllocNode:Sid=amarel1:22391
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=hal0053
   BatchHost=hal0053
   NumNodes=1 NumCPUs=28 NumTasks=28 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=28,mem=124000M,node=1
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=124000M MinTmpDiskNode=0
   Features=(null) Gres=(null) Reservation=(null)
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/scratch/gc563/run.STMV.CPU.slurm
   WorkDir=/scratch/gc563
   StdErr=/scratch/gc563/slurm.%N.244348.out
   StdIn=/dev/null
   StdOut=/scratch/gc563/slurm.%N.244348.out
   Power=
```


If your jobs have already completed (or have been terminated), you can see details about those jobs using the sacct command with your NetID and a start time for the list of jobs this command will produce.

```sacct --user=[NetID] --starttime=2018-07-03 --format=JobID,Partition,JobName,MaxRSS,NodeList,Elapsed,MaxDiskRead,MaxDiskWrite,State```

# Killing/ cancelling/ terminating jobs

To terminate a job, regardless of whether it is running or just waiting in the job queue, use the scancel command and specify the JobID number of the job you wish to terminate:  
```scancel 1633383```
A job can only be cancelled by the owner of that job. When you terminate a job, a message from the SLURM workload manager will be directed to STDERR and that message will look like this:  
```slurmstepd: *** JOB 1633383 ON slepner036 CANCELLED AT 2016-10-04T15:38:07 ***```

# Installing your own software

Package management systems like yum or apt-get, which are used to install software in typical Linux systems, are not available to users of shared computing resources like Amarel. Thus, most packages need to be compiled from their source code and then installed. Further, most packages are generally configured to be installed in /usr or /opt, but these locations are inaccessible to (not writeable for) general users. Special care must be taken by users to ensure that the packages will be installed in their own /home directory (/home/[NetID]).

As an example, here are the steps for installing ZIPPER, a generic example package that doesn’t actually exist:

-  Download your software package. You can usually download a software package to your laptop, and then transfer the downloaded package to your /home/[NetID] directory on Amarel for installation. Alternatively, if you have the http or ftp address for the package, you can transfer that package directly to your home directory while logged-in to Amarel using the **wget** utility:
`wget http://www.zippersimxl.org/public/zipper/zipper-4.1.5.tar.gz
`

-  Unzip and unpack the .tar.gz (or .tgz) file. Most software packages are compressed in a .zip, .tar or .tar.gz file. You can use the tar utility to unpack the contents of these files:  
```tar -zxf zipper-4.1.5.tar.gz```

-  Read the instructions for installing. Several packages come with an INSTALL or README script with instructions for setting up that package. Many will also explicitly include instructions on how to do so on a system where you do not have root access. Alternatively, the installation instructions may be posted on the website from which you downloaded the software.  
```
  cd zipper-4.1.5 
  less README
```

-  Load the required software modules for installation. Software packages generally have dependencies, i.e., they require other software packages in order to be installed. The README or INSTALL file will generally list these dependencies. Often, you can use the available modules to satisfy these dependencies. But sometimes, you may also need to install the dependencies for yourself. Here, we load the dependencies for ZIPPER:  
```module load intel/16.0.3 mvapich2/2.1```

-  Perform the installation. The next few steps vary widely but instructions almost always come with the downloaded source package. Guidance on the special arguments passed to the configure script is often available by running the **./configure -–help** command. What you see below is just a typical example of special options that might be specified.  
```  
./configure --prefix=/home/[NetID]/zipper/4.1.5 --disable-float --enable-mpi --without-x --disable-shared
make -j 4
make install
```  
  Several packages are set up in a similar way, i.e., using configure, then make, and make install. Note the options provided to the configure script – these differ from package to package, and are documented as part of the setup instructions, but the prefix option is almost always supported. It specifies where the package will be installed. Unless this special argument is provided, the package will generally be installed to a location such as /usr/local or /opt, but users do not have write-access to those directories. So, here, I'm installing software in my /home/[NetID]/zipper/4.1.5 directory. The following directories are created after installation:   
      - `/home/[NetID]/zipper/4.1.5/bin` where executables will be placed
      - `/home/[NetID]/zipper/4.1.5/lib` where library files will be placed
      - `/home/[NetID]/zipper/4.1.5/include` where header files will be placed
      - `/home/[NetID]/zipper/4.1.5/share/man` where documentation will be placed

-  Configure environment settings. The above bin, lib, include and share directories are generally not part of the shell environment, i.e., the shell and other programs don’t “know” about these directories. Therefore, the last step in the installation process is to add these directories to the shell environment:  
```
export PATH=/home/[NetID]/zipper/4.1.5/bin:$PATH
export C_INCLUDE_PATH=/home/[NetID]/zipper/4.1.5/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/home/[NetID]/zipper/4.1.5/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/home/[NetID]/zipper/4.1.5/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/[NetID]/zipper/4.1.5/lib:$LD_LIBRARY_PATH
export MANPATH=/home/[NetID]/zipper/4.1.5/share/man:$MANPATH
```
These **export** commands are standalone commands that change the shell environment, but these new settings are only valid for the current shell session. Rather than executing these commands for every shell session, they can be added to the end of your ~/.bashrc file which will result in those commands being executed every time you log-in to Amarel.

# Singularity

[Singularity](http://singularity.lbl.gov/) is a Linux containerization tool suitable for HPC environments. It uses its own container format and also has features that enable importing Docker containers.

[Docker](https://www.docker.com/) is a platform that employs features of the Linux kernel to run software in a container. The software housed in a Docker container is not standalone program but an entire OS distribution, or at least enough of the OS to enable the program to work. Docker can be thought of as somewhat like a software distribution mechanism like yum or apt. It also can be thought of as an expanded version of a chroot jail, or a reduced version of a virtual machine.

##Important differences between Docker and Singularity:

  * Docker and Singularity have their own container formats.
  * Docker containers can be imported and run using Singularity.
  * Docker containers usually run as root, which means you cannot run Docker on a  shared computing system (cluster).
  * Singularity allows for containers that can be run as a regular user. How? When importing a Docker container, Singularity removes any elements which can only run as root. The resulting containers can be run using a regular user account.

##Importing a Docker image:

If you have a pre-built Docker container, you can use Singularity to convert this container to the Singularity format. Once that's done, you can upload your Singularity container to your storage space on Amarel and run jobs using that container.

Here's an example. NOTE that most of these steps are performed on your local system, not while logged-in on Amarel.

If you need to use any of Amarel's filesystems inside your container, you will need to make sure the appropriate directories exist in your container so those filesystems can be mounted using those directories.

Start your container (in this example we will use ubuntu:latest) and create directories for mounting /scratch/gc563 and /projects/oarc. Of course, you'll need to use directories that you can access on Amarel.

```
$ sudo docker run -it ubuntu:latest bash  
root@11a87dkw8748:/# mkdir -p /scratch/gc563 /projects/oarc
```

##Exporting your Docker image

Find the name of your Docker image using the 'docker ps' command,
```
$ sudo docker ps
CONTAINER ID  IMAGE          COMMAND  CREATED        STATUS        NAMES
11a87dkw8748  ubuntu:latest  "bash"   2 minutes ago  Up 2 minutes  bendakaya_pakodi
```
In this example the name of the images is bendakaya_pakodi. Export this image to a tarball,  
```
$ sudo docker export bendakaya_pakodi > ubuntu.tar
```

##Converting to a Singularity image
You will need to have Singularity installed on your local workstation/laptop to prepare your image. The 'create' and 'import' operations of Singularity require root privileges, which you do not have on Amarel.

Create an empty singularity image, and then import the exported docker image into it,  
```
$ sudo singularity create ubuntu.img
Creating a sparse image with a maximum size of 1024MiB...
Using given image size of 1024
Formatting image (/sbin/mkfs.ext3)
Done. Image can be found at: ubuntu.img
$ sudo singularity import ubuntu.img ubuntu.tar
```

##Using Singularity containers inside a SLURM job

[Transfer](https://rutgers-oarc.github.io/amarel/#movingfiles) your new Singularity image to Amarel. The following steps are performed while logged-in to Amarel.

You can run any task/program inside the container by prefacing it with  
```singularity exec [your image name]``` 

Here is a simple example job script that executes commands inside a container,  
```
#SBATCH --partition=main             # Partition (job queue)
#SBATCH --job-name=sing2me           # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=4000                   # Real memory (RAM) required (MB)
#SBATCH --time=00:30:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file

module purge
module load singularity/.2.5.1

## Where am I running?
srun singularity exec ubuntu.img hostname

## What is the current time and date?
srun singularity exec ubuntu.img date
```

If you created directories for any Amarel filesystems, you should find they are mounted inside your container,
```
mount | grep gpfs
/dev/scratch/gc563 on /scratch/gc563 type gpfs (rw,relatime)
/dev/projects/oarc on /projects/oarc type gpfs (rw,relatime)
```
NOTE: If your container mounts Amarel directories, software inside the container may be able to destroy data on these filesystems for which you have write permissions. Proceed with caution.


# Troubleshooting/ Common Problems 

Failure to load module dependencies/prerequisites:
```
 module load R-Project/3.4.1
Lmod has detected the following error:  These module(s) exist but cannot be loaded as requested: "R-Project/3.4.1"
Try: "module spider R-Project/3.4.1" to see how to load the module(s).

```
This software module has a prerequisite module that must be loaded first. To find out what prerequisite module is required, use the 'module spider' command followed by the name of the module you're trying to load:  
```
module spider R-Project/3.4.1
    This module can only be loaded through the following modules:
      intel/17.0.4
    Help: 
      This module loads the installation R-Project 3.4.1 compiled with the Intel 17.0.4 compilers.
```
 Ah-ha, it looks like the intel/17.0.4 module must be loaded before loading R-Project/3.4.1



# Acknowledging Amarel

Please reference OARC and the Amarel cluster in any research report, journal or publication that requires citation of an author's work. Recognizing the OARC resources you used to conduct your research is important for our process of acquiring funding for hardware, support services, and other infrastructure improvements. The minimal content of a reference should include:

*Office of Advanced Research Computing (OARC) at Rutgers, The State University of New Jersey*

A suggested acknowledgement is:

*The authors acknowledge the Office of Advanced Research Computing (OARC) at Rutgers, The State University of New Jersey for providing access to the Amarel cluster and associated research computing resources that have contributed to the results reported here. URL: http://oarc.rutgers.edu*
















                    







