# Bash cheatsheet (command line): 

|command | description | usage example |
|:-----------|:--------|:-------------|
|ls -ltra mydir| list all files in mydir in descending order of creation, with permissions | ls -ltra .|
|du | disk usage, e.g. how much space does your directory occupy, -h human-readable | du -h mydir|
|which <command>| see where command is installed | which python|
|pwd| which directory I'm in | pwd |
|man <command>| manual page for command | man cut|
|grep <pattern>| filter for lines which fit pattern | cat myfile &#124; grep GATK |
|cut -d<delimiter> -f<number>| split line by delimiter and get field number 3| cat myfile &#124; cut -d'_' -f3 |
|sort <file>| sort lines, often used with `uniq` | sort myfile &#124; uniq |
|uniq| suppress repeated lines, works only if sorted | see above example |
|less | paginated output | less myfile |
| >| redirect output (e.g. list files and save filenames in aaa.txt) | ls > aaa.txt |
|>>| append output to existing file | echo "blah" >> aaa.txt |
| find| find files with some properties e.g. display all files recursively from current directory| `find .`|
| chmod| change permissions on a file or directory, eg. make myscript.sh executable for user |chmod u+x myscript.sh|
|top| display most intensive processes | `top`|
|ps auxw| list processes | `ps auxw` |
|time myscript.sh| measure how much time does myscript.sh take to finish execution 
|`echo` $variable| output the value of the variable  |echo $USER |

For a more complete cheatsheet including ifs, loops and functions, see [this website](https://devhints.io/bash)
Here is a sample: 
```
dirname /home/kp807/projects/cluster_reports/cluster.csv       # everything but the last part
basename /home/kp807/projects/cluster_reports/cluster.csv      # last part of path
CURRENT_DIR=`pwd` ; echo $CURRENT_DIR                          # backtick for execution of bash command
echo 'projects_cluster_file.csv' | cut -d '_' -f1              # split name to retain a part of file
echo $((1 + 2))                                                # double parenthesis for arithmetic expressions
echo filename_fly{5..10}.csv                                   # list comprehensions   {start..end}

# for-loop: 
for file in filename_fly{5..10}.csv; do  echo $file ; done     

# if-statement - 2 examples
if [ 1 -gt 2 ]; then  echo '1 > 2' ; else echo '1<2' ; fi      
if [ -d "newdir" ]; then  echo 'directory exists' ; else echo 'this directory doesnt exist' ; fi                        
#variable assignment
a=10; echo $a        #good - no spaces
b = 10; echo $b      #bad - spaces around =
```



# Environment modules (lmod) cheatsheet

**Command**|**Description**
---|---
`module avail`|  show a list of the core (primary) modules available
`module load *module_name*`| loads the named software module 
`module use` /projects/community/projectsmodule`| Now if you do module avail, it also shows modules created by other users 
`module spider` | shows a comprehensive list of all available modules or type name after spider to show details about specific module
`module keyword` | looks for a keyword in the description of a module
`module purge`| removes all loaded module
`ml`| lists modules loaded 

To learm more about environments and modules, see [lmod documenation](https://lmod.readthedocs.io/en/latest/). 

# Slurm cheatsheet

Command | Description
---|---
`sinfo`                                |info about utilization and resources
`sbatch my_slurm_script.sh`            |submit a slurm script
`srun -N 1 -c 2 -t 1:00:00 --pty bash` |run interactive shell on a compute node for 1 hour
`scancel  123456`                      |cancel job 123456 that's currently running or queued
`sacct -u kp807 --start=2018-05-02`    |show all jobs by user kp807 since some date
`scontrol show job 123456 -dd`         |details about job that's running 
`scontrol show node hal0001`           |details about node `ha0001`
`scontrol --help`                      |generic way to invoke quick help on any command
`sacct -o MaxRSS,job,jobName,state`    | keep track of how much memory you are using
`sacct -o Elapsed,Start,End,job,jobName` | estimate the wall time

To learn more about slurm, see [slurm documentation](https://slurm.schedmd.com/)

# Check quota cheatsheet

NOTE: `/home` and `/scratch` are two different filesystems, with different backup policies. 

Command | Description of which file usage
---|---
`mmlsquota scratch --block-size=auto` | file usage on `/scratch/netid` for the user's netid
`mmlsquota home --block-size=auto`    | file usage on `/home/netid` for the user's netid
`mmlsquota home:foran --block-size=auto` | file usage in the shared folder foran (/projects/foran)
`mmlsquota -j foran home  --block-size=auto` | quota and usage of the whole fileset foran
`du -hs /directory/to/query/*`          | human-readable sizes of all 1st-level subdirectories of `/directory/to/query/`

# Linux basic commands

## File tree

| **Command**  | **Meaning** |
|------------:|:----:|
|*ls*     |list files and directories|
|*ls -a*	|list all files and directories|
|*mkdir*  | make a directory|
|*cd directory* |	change to named directory|
|*cd*|    change to home-directory|
|*cd ~* |	change to home-directory|
|*cd ..* |	change to parent directory|
|*pwd*	|     display the path of the current directory |

## Moving and viewing files

|**Command**|**Meaning**|
|------------:|:----------|
|cp file1 file2 |	copy file1 and call it file2|
|mv file1 file2|	move or rename file1 to file2|
|rm file	|remove a file|
|rmdir directory |	remove a directory (only if empty)|
|cat file |	display a file|
|less file|	display a file a page at a time|
|head -19 file	|display the first 19 lines of a file|
|tail -19 file	|display the last 19 lines of a file|
|grep Finally myfile.txt	|search myfile.txt for word Finally|
|wc file |	count number of lines/words/characters in file

## Piping

|**Command**	| **Description**
---------|--------------
command > file |	redirect standard output to a file
command >> file	|append standard output to a file
command < file|	redirect standard input from a file
*command1 &#x7c; command2*	| pipe the output of command1 to the input of command2
cat file1 file2 > file0	| concatenate file1 and file2 to file0
sort |	sort data
who|	list users currently logged in
*	|match any number of characters
?	|match one character
man command|	read the online manual page for a command
whatis command	|brief description of a command
apropos keyword|	match commands with keyword in their man pages

## Permissions

**Symbol** |**Description**
---|---
u |user
g |group
o |other
a |all
r |read
w |write (and delete)
x |execute (and access directory)
u+x |add executing permission for user
u-x |take away executing permission

Output of `ls -ltra`  
``` drwxr-x--x  12 kp807 kp807      4096 Jun 26 23:33 .rstudio ```

- d = it's a directory
- rwx = first group of permissions, readable, writeable and executable by user
- r-x = second group of permissions, readable and executable by group, but not writeable by group
- --x = third group of permissions, executable by others, but not readable and writeable by others
- .rstudio = hidden directory

## Linux on the web

- [Galen's tutorial](http://www.rci.rutgers.edu/~gc563/linux/index.html)
- [Excellent short lessons](http://www.ee.surrey.ac.uk/Teaching/Unix/)
- [Alexei's course he runs for engineering students every semester](http://linuxcourse.rutgers.edu/)
- [Cornell virtual workshop](https://cvw.cac.cornell.edu/Linux/) - interactive if you have xsede login; many HPC topics available
- [Software carpentry](https://software-carpentry.org/lessons/) - list of introductory programming courses


