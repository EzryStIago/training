# Navigating file tree

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
|*rm* | removes files only|
|*rmdir* | removes directory(make sure its empty)|

# Viewing files

|**Command**|**Meaning**|
|------------:|:----------|
|cp file1 file2 |	copy file1 and call it file2|
|mv file1 file2|	move or rename file1 to file2|
|rm file	|remove a file|
|rmdir directory |	remove a directory|
|cat file |	display a file|
|less file|	display a file a page at a time|
|head file	|display the first few lines of a file|
|tail file	|display the last few lines of a file|
|grep 'keyword'| file	search a file for keywords|
|wc file -l/w |	count number of lines/words/characters in file

# Useful commands: 

|command | description | usage example |
|:-----------|:--------|:-------------|
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
| chmod| change permissions on a file or directory |chmod u+x myscript.sh|
|top| display most intensive processes | `top`|
|ps auxw| list processes | `ps auxw` |


# Redirecting

|**Command**	| **Meaning**
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

# Permissions

**Symbol** |**Meaning**
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

# Modules
**Command**|**Meaning**
---|---
module avail|  show a list of the core (primary) modules available
module load *module_name*| loads the named software module 
module use `/projects/community/projectsmodule`| Now if you do module avail, it also shows modules created by other users 
module spider | shows a comprehensive list of all available modules or type name after spider to show details about specific module
module keyword | looks for a keyword in the description of a module
module purge| removes all loaded module
ml| lists modules loaded 

To learm more about environments and modules, see [lmod documenation](https://lmod.readthedocs.io/en/latest/). 
