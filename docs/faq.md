## error and output files

**In the directory "slurm.slepner009.3695956.err" was created. Is this what should be happening?**   
Yes, if one of your options to run a slurm job included a line like `--error=slurm.%N.%j.err` then any error messages during execution will be going to a file you can inspect afterwards to
determine how to fix the job's failure. E.g. if your jobs's results are supposed to go to `results` directory, but the directory "results" was not created before running the job. Since the slurm script is accessing a nonexistent path, it aborts. You can't see the job with `squeue` because it's not running (it finished with "exit 1" i.e. failed). You can use the command `sacct` (see slurm cheatsheet) to see some details of finished jobs (but not as many details as with `scontrol show job 12345`). 

## common errors and their fixes

**I'm getting this error. Why? What's happening?**

|error|explanation| fix|
|---|---|---|
|`sbatch: error: Slurm temporarily unable to accept job, sleeping and retrying.`| The Slurm queue reached its maximum and is not accepting new jobs| Wait until the queue clears up and resubmit|
|`Batch script contains DOS line breaks (\r\n)sbatch: error: instead of expected UNIX line breaks (\n).`|If you created or edited the file in Windows, it can insert Windows line breaks, which are not the same as unix line breaks.| convert between the two formats with `dos2unix` and `unix2dos` commands  - see [this explanation](https://www.geeksforgeeks.org/dos2unix-unix2dos-commands/)|
|`No such file or directory`| A path you gave is wrong or nonexistent| Check the exact path; if file exists, check permissions|

## numpy and other hidden modules

**How do I get numpy?**

You need to load libraries that are needed by numpy first. For example Intel Math Kernel Library: 
```
module load intel_mkl/17.0.2 python/3.5.2 
```

## R DLLpath 

- **When I install R package xxx, I get an error  dyn.load(file, DLLpath = DLLpath, ...):** 

A lot of R packages are actually written in other languages, very often C or C++. Depending on how the R package is packaged, it can rely on "shared objects", recognized by .so extension. These .so files are 
created by compilers, but different compiler versions will produce files which are sometimes not compatible. So either the .so is not found because its file path is missing, or it has been compiled for a 
different version of the operating system or different version of the software. Complicating this, sometimes software which is in your path (e.g. anaconda will put itself in your .bashrc so it will be 
available) has its own version of the library you want to load - e.g. R package may need a .so file that anaconda installation already contains, but is not compatible with your R package. Solutions can vary: 

- remove paths to known other softeware (such as Anaconda)
- remove other R packages that might be interfering (e.g. by moving it out of your path, not necessarily deleting)


## email notifications

**How do I make sure that I get an email when my slurm job finally starts, finishes, or fails?**

It is possible to specify in the slurm batch script or srun command that you wish to receive an email notification. There are two options associated with this: 
`mail-type` and `mail-user`. For example, ``--mail-type=BEGIN`` and ``--mail-user=kp807@rutgers.edu`` . 

- Without mail-type specified, it doesn't send an email 
- If you have domain other than `rutgers.edu`, slurm doesn't properly send it (e.g. `kp807@oarc.rutgers.edu` won't work, but `kp807@rutgers.edu` will.)

Other options for `mail-type` are: BEGIN (send email when job starts, useful for interactive jobs), END (send email when job finishes), FAIL (send email when job finishes), ALL(send email for all 3). 

## file permissions 
 
One way to control the access to the dataset is via `getfacl` and `setfacl` commands (see https://www.computerhope.com/unix/usetfacl.htm). Here is an example of how to give permission to user user1, and only
to this user, to read one of my subdirectories, `sharedir`: 

``` setfacl -m u:user1:rx /home/kp807/sharedir
```
Here is how you would see what are the permissions to this directory: 
``` getfacl  /home/kp807/sharedir
```
- `-m` option means modify permissions
- `u`  means user
- `rx` means read and execute. You must give execute permissions to a directory if you want to execute a listing of the directory. You also have to give `x` permission to every directory that needs to be
traversed during execution and reading. E.g. you'd also have to execute `setfacl -m u:user1:x /home/kp807/`

NOTE OF CAUTION: By default, Linux umask is set to 022. This means that by default, all directories and files are readable by group and other. The only thing stopping someone else reading your home directory
is lack of `x` flag on the home directory. If you add this flag to a user, now this user will also be able to read any readable files in your home directory. So use `x` flag judiciously, or simply change
all your permissions (and user's `umask` in `.bashrc` to something more restrictive). 


