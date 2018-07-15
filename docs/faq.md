
- **In the directory "slurm.slepner009.3695956.err" was created. Is this what should be happening?**   
Yes, if one of your options to run a slurm job included a line like `--error=slurm.%N.%j.err` then any error messages during execution will be going to a file you can inspect afterwards to
determine how to fix the job's failure. E.g. if your jobs's results are supposed to go to `results` directory, but the directory "results" was not created before running the job. Since the slurm script is accessing a nonexistent path, it aborts. You can't see the job with `squeue` because it's not running (it finished with "exit 1" i.e. failed). You can use the command `sacct` (see slurm cheatsheet) to see some details of finished jobs (but not as many details as with `scontrol show job 12345`). 

- **I'm getting this error. Why? What's happening?**

|error|explanation| fix|
|---|---|---|
|`sbatch: error: Slurm temporarily unable to accept job, sleeping and retrying.`| The Slurm queue reached its maximum and is not accepting new jobs| Wait until the queue clears up and resubmit|
|`Batch script contains DOS line breaks (\r\n)sbatch: error: instead of expected UNIX line breaks (\n).`|If you created or edited the file in Windows, it can insert Windows line breaks, which are not the same as unix line breaks.| convert between the two formats with `dos2unix` and `unix2dos` commands  - see [this explanation](https://www.geeksforgeeks.org/dos2unix-unix2dos-commands/)|
|`No such file or directory`| A path you gave is wrong or nonexistent| Check the exact path; if file exists, check permissions|



