# Things affected by the outage Oct 6 - Oct 11

Please see our [system status page](https://oarc.rutgers.edu/systems-status/) for some explanation of what went on. Here we will focus on some changes and what it means, from the technical standpoint (please note we are short-staffed until Friday, Oct 19): 

### Things that have been fixed:
- intel compiler's licence path


### Things that we know are broken:
- GPU nodes are going to be rebooted to fix some things
- gcc doesn't compile on login node
- lmod version reverted to 6.0.1.

** My /scratch disappeared. Where is it? **

The filesystem has been moved to new hardware. As a result, what used to be `/scratch` before the upgrade is now `/oldscratch`. Please copy any files you need to retain to new /scratch. `/oldscratch` will be retired soon. 

** I get a message telling me someone may be spoofing, when I try to login to Perceval. Am I in danger? **

You are not in danger. During the upgrade the fingerprint of the server changed. Go to your `.ssh/authorized_hosts` or `.ssh/known_hosts` file and delete the line containing perceval information. The message you got before, that looked like the following, will disappear next time you login. 
```
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@       WARNING: POSSIBLE DNS SPOOFING DETECTED!          @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
The ECDSA host key for perceval.hpc.rutgers.edu has changed,
and the key for the corresponding IP address 172.16.94.50
is unknown. This could either mean that
DNS SPOOFING is happening or the IP address for the host
and its host key have changed at the same time.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@    WARNING: REMOTE HOST IDENTIFICATION HAS CHANGED!     @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
IT IS POSSIBLE THAT SOMEONE IS DOING SOMETHING NASTY!
Someone could be eavesdropping on you right now (man-in-the-middle attack)!
It is also possible that a host key has just been changed.
```
** What is the status of GPU nodes? ** 


To check the current status of GPU nodes with reason: check with command `sinfo -R -l -p gpu`. Here is some output: 
```
[kp807@perceval1 ~]$ sinfo -R -l -p gpu
Mon Oct 15 10:06:40 2018
REASON               USER         TIMESTAMP           STATE  NODELIST
NHC: check_cmd_outpu root(0)      2018-10-12T03:22:03 drain  cuda[001-008]
```
```
[kp807@amarel1 ~]$ sinfo -R -l -p gpu
Mon Oct 15 10:06:19 2018
REASON               USER         TIMESTAMP           STATE  NODELIST
no cnx ipmi          root(0)      2018-10-13T11:51:15 drain* gpu[001-003]
IP                   root(0)      2018-10-15T09:53:10 drng   pascal[001-004]
not part of mentat2  root(0)      2018-10-13T11:52:52 drain  gpu[005-006]
IP                   root(0)      2018-10-15T09:59:24 drain  pascal005
IP                   root(0)      2018-10-15T09:53:10 drain  pascal006
```
 
** Where are snapshots for /home stored? ""

Normally, they are in `/home/.snapshots`. After the upgrade the .snapshots was not copied. New `/home/.snapshots` will be created as soon as the backups kick in - at latest on Oct 20, but more likely midweek ~ Oct 18. Meanwhile the current version of snapshots as 10/09/2018 is available with a manual assistance  from OARC team on request.

# Usual FAQs

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

## python and conda

Backgroud:

- bash has some predefined variables that programs (or you) can use. You can get a list of those variables by typing `env` command in your bash terminal. By convention, the names are uppercase.
- to set variables and get variables use `=` (without spaces) and `echo` with `${NAME_OF_VARIABLE}`. E.g. `PWD=/home/kris` to set the `PWD` variable, or `echo $PATH` to print the value of variable `PATH`
- `PATH` is a bash variable that tells the systems in which directories to look for the programs/commands you invoke from the bash terminal. You can see the current path by `echo $PATH`
- .bashrc file is the file that always gets called when you start up a terminal, and sets PATH variable to some default

*How to install anaconda with python 2, but not set it as default anaconda version*

1. install anaconda2 (https://www.anaconda.com/download/) - follow instructions for your OS
2. anaconda puts itself automatically in your .bashrc file. You need to find the line that anaconda added and comment it out or delete from your .bashrc file so that anaconda2 is not your default every time you log in. The lines to comment out look like this: 
```
# added by Anaconda2 installer
export PATH="/home/kris/anaconda2/bin:$PATH"
```
3. When you want to use anaconda2 instead of anaconda3, you reset your PATH variable manually, by issuing this command in your bash terminal:
```export PATH="/home/$USER/anaconda2/bin:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin"```
This will be valid ONLY in that terminal, so remember to launch any python2 scripts in that terminal, otherwise python will default to your other default
4. Verify which conda and python version are you using by this command: ```which python```
5. When you are sure that you are using conda from anaconda2, create new conda environment with a descriptive name, e.g.
```conda create --name mytensorflow-1.5```
6. Activate this environment:
```source activate mytensorflow-1.5```
7. Install tensorflow-1.5 into this environment:
```pip install tensorflow==1.5.0```
8. After having finished your work, ```source deactivate``` gets you out of the mytensorflow-1.5 environment.

Now you have installed your environment. Next time, these are the commands to get you in the python2, mytensorflow-1.5 environment ready to do your python2 work:
```
export PATH="/home/$USER/anaconda2/bin:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin"
source activate mytensorflow-1.5
```
