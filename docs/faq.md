
- **In the directory "slurm.slepner009.3695956.err" was created. Is this what should be happening?**   
Yes, if one of your options to run a slurm job included a line like `--error=slurm.%N.%j.err` then any error messages during execution will be going to a file you can inspect afterwards to
determine how to fix the job's failure. E.g. if your jobs's results are supposed to go to `results` directory, but the directory "results" was not created before running the job. Since the slurm script is accessing a nonexistent path, it aborts. You can't see the job with `squeue` because it's not running (it finished with "exit 1" i.e. failed). You can use the command `sacct` (see slurm cheatsheet) to see some details of finished jobs (but not as many details as with `scontrol show job 12345`). 

- **I'm getting this error. Why? What's happening?**

|error|explanation| fix|
|---|---|---|
|`sbatch: error: Slurm temporarily unable to accept job, sleeping and retrying.`| The Slurm queue reached its maximum and is not accepting new jobs| Wait until the queue clears up and resubmit|
|`Batch script contains DOS line breaks (\r\n)sbatch: error: instead of expected UNIX line breaks (\n).`|If you created or edited the file in Windows, it can insert Windows line breaks, which are not the same as unix line breaks.| convert between the two formats with `dos2unix` and `unix2dos` commands  - see [this explanation](https://www.geeksforgeeks.org/dos2unix-unix2dos-commands/)|
|`No such file or directory`| A path you gave is wrong or nonexistent| Check the exact path; if file exists, check permissions|

- **How do I get numpy?**

You need to load libraries that are needed by numpy first. For example Intel Math Kernel Library: 
```
module load intel_mkl/17.0.2 python/3.5.2 
```

- **When I install R package xxx, I get an error  dyn.load(file, DLLpath = DLLpath, ...):** 

A lot of R packages are actually written in other languages, very often C or C++. Depending on how the R package is packaged, it can rely on "shared objects", recognized by .so extension. These .so files are 
created by compilers, but different compiler versions will produce files which are sometimes not compatible. So either the .so is not found because its file path is missing, or it has been compiled for a 
different version of the operating system or different version of the software. Complicating this, sometimes software which is in your path (e.g. anaconda will put itself in your .bashrc so it will be 
available) has its own version of the library you want to load - e.g. R package may need a .so file that anaconda installation already contains, but is not compatible with your R package. Solutions can vary: 

- remove paths to known other softeware (such as Anaconda)
- remove other R packages that might be interfering (e.g. by moving it out of your path, not necessarily deleting)



