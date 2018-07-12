*Useful resources contributed or setup by members of the Amarel Research Computing Community for the benefit of fellow members*

## Community-Contributed Software Modules
___
Users who wish to share software with the Amarel community of researchers can do so using the guidelines and suggestions presented here.

## Log-in to Amarel or Perceval and check out /projects/community

Have you created a particularly fast build of your favorite simulation tool? Have you compiled a packaged that's notoriously hard to compile? Want to share some software with collaborators within the Amarel user community?

OARC has setup a dedicated repository for community-contributed software and associated modulefiles in /projects/community. All Amarel and Perceval users can contribute software to this repository. We have established some basic guidelines to help keep things organized.

If you would like to see the community-contributed software packages when you run the 'module avail' or 'module spider' command, add /projects/community/modulefiles to your MODULEPATH,  
`export MODULEPATH=$MODULEPATH:/projects/community/modulefiles`    
or add that line to your ~/.bashrc file for that setting to persist.

## A Note About Governance
Please follow the guidelines presented here. The OARC research support team will remove or edit contributions that do not conform to these guidelines. If you find any problems or if you have questions, please let us know by sending a message to <help@hpc.rutgers.edu>

## User Support
Contributed software is officially unsupported. However, individual users who setup contributed software may be willing to answer basic questions about usability, performance, or selected build options. We certainly want to encourage users to contribute software, but at the same time, we do not want to burden busy researchers with support expectations. Be considerate if you contact a user about contributed software.

## Permissions
We've adjusted permissions (and set the sticky bit) so Amarel and Perceval users can write to `/projects/community` and its subdirectories.
```
# chmod 777 /projects/community
# chmod o+t /projects/community
```
This means that only a file's owner, a directory's owner, or root can rename or delete those files or directories, thus preventing others from modifying your contributions. You can adjust permissions of the files and directories you create, which may be necessary to enable others to use the software you place there.

## Naming Conventions

New software packages should be added using the following conventions.

Software Packages| Directory Structure
---|---
Modulefiles|  `/projects/community/modulefiles/compiler-name/version/mpi-name/version/CUDA/version/pkg-name/version/NetID`
Core packages| `/projects/community/pkg-name/version/NetID`
Compiler-dependent packages|  `/projects/community/compiler-name/version/pkg-name/version/NetID`
Compiler-dependent packages using CUDA|`/projects/community/compiler-name/version/cuda/version/pkg-name/version/NetID`
MPI-Compiler-dependent packages|`/projects/community/compiler-name/version/mpi-name/version/pkg-name/version/NetID`
MPI-Compiler-dependent packages using CUDA|`/projects/community/compiler-name/version/mpi-name/version/cuda/version/pkg-name/version/NetID`

## Good Housekeeping and Reproducibility in Research
 
(1) All contributed software packages MUST have a README.OARC file created and stored in the directory containing the installed software. This file must include your contact info and any useful information about the software. This file should include all commands used to build the package since this information is important for anyone wishing to use contributed software for research purposes.

There is a long but useful example README.OARC file in `/projects/community/gcc/7.3.0/gc563`

(2) All contributed software packages MUST have an associated [Lmod modulefile](http://lmod.readthedocs.io/en/latest) so users can easily access the software. There are details for doing this and examples below.

## An Example: GCC-7.3.0

March 23, 2018: I want to build GCC, The GNU Compiler Collection, version 7.3.0 since a project I'm working on requires it, but GCC-7.3.0 is not already setup on Amarel. I'll begin by "bootstrapping" using the core/default tools that come with the version of CentOS installed on Amarel. That version (at the time of this writing) is
```
$ cat /etc/*release | head -n 1
CentOS Linux release 7.4.1708 (Core)
```
plus the additional components and tools that come from the Enterprise Linux repository. This version of CentOS is based on the following kernel and GCC version:
```
$ cat /proc/version
Linux version 3.10.0-693.21.1.el7.x86_64 (builder@kbuilder.dev.centos.org) (gcc version 4.8.5 20150623 (Red Hat 4.8.5-16) (GCC) ) #1 SMP Wed Mar 7 19:03:37 UTC 2018
```
Installing GCC requires a few prerequisite packages: these must be installed and configured for use before installing GCC. Those prerequisites are: GNU Multiple Precision Library (GMP) version 4.3.2 (or later) MPFR Library version 2.4.2 (or later) MPC Library version 0.8.1 (or later) I'll download, test, and install a recent version of each of these software packages. Since I'm just trying to get setup for installing GCC here, I'll skip the description of these steps because I'll discuss those steps in detail later (below).
```
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2
tar -jxf gmp-6.1.2.tar.bz2
./configure --prefix=/projects/community/gmp/6.1.2/gc563 --enable-cxx --enable-fft
make -j 4
make check
make install
cd .. ; rm -rf gmp-6.1.2*
```
**Note:** the following environment settings are needed for installing the next package. These environment settings have not yet been set for my shell session, so that's why I'm not prepending them to existing path settings here:
```
export C_INCLUDE_PATH=/home/gc563/gmp/6.1.2/include
export CPLUS_INCLUDE_PATH=/projects/community/gmp/6.1.2/gc563/include
export LIBRARY_PATH=/projects/community/gmp/6.1.2/gc563/lib
export LD_LIBRARY_PATH=/projects/community/gmp/6.1.2/gc563/lib
```
```
wget http://www.mpfr.org/mpfr-current/mpfr-4.0.1.tar.bz2
tar -jxf mpfr-4.0.1.tar.bz2 
cd mpfr-4.0.1
./configure --prefix=/projects/community/mpfr/4.0.1/gc563 --enable-thread-safe
make -j 4
make install
cd .. ; rm -rf mpfr-4.0.1*
```
**Note:** the following environment settings have already been set, so that's why I'm prepending additional segments to existing path settings here:
```
export C_INCLUDE_PATH=/home/gc563/mpfr/4.0.1/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/projects/community/mpfr/4.0.1/gc563/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/projects/community/mpfr/4.0.1/gc563/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/projects/community/mpfr/4.0.1/gc563/lib:$LD_LIBRARY_PATH
```
```
wget https://ftp.gnu.org/gnu/mpc/mpc-1.1.0.tar.gz
tar -zxf mpc-1.1.0.tar.gz
cd mpc-1.1.0
./configure --prefix=/projects/community/mpc/1.1.0/gc563
make
make install
cd .. ; rm -rf mpc-1.1.0*
```
```
export C_INCLUDE_PATH=/home/gc563/mpc/1.1.0/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/projects/community/mpc/1.1.0/gc563/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/projects/community/mpc/1.1.0/gc563/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/projects/community/mpc/1.1.0/gc563/lib:$LD_LIBRARY_PATH
```
Now, I can finally start building GCC. I'll download the latest stable supported release of GCC, open the tarball, and get ready to install it:
```
wget http://mirrors.concertpass.com/gcc/releases/gcc-7.3.0/gcc-7.3.0.tar.gz
tar -zxf gcc-7.3.0.tar.gz
cd gcc-7.3.0
```
Like the 3 prerequisite software packages I just installed, the GCC compiler suite uses the traditional configure/make/make-install installation procedure. The first step is to run the Bash script named 'configure'. The configure script inspects the existing hardware and/or software configuration of the machine where you are about to install this new software. If something is needed and missing (e.g., a needed library or access to specific tools), the configure script should let you know. The result of running this script is the creation of a file named 'Makefile' that contains the instructions for compiling your new software based on the findings from the configure script. In addition, when you run the configure script, you have the opportunity to select important options for how your new software will be setup. Most importantly, you can select where the software will be installed. This is accomplished using the '--prefix=[some location]' option. In my case, I'll specify a location in my /home directory. There are many other options and those options vary depending on what you're installing. To see a summary of the available options, use './configure --help' and for details about what those options mean, see the documentation for your new software.
```
$ ./configure --prefix=/projects/community/gcc/7.3.0/gc563--with-  
mpc=/projects/community/mpc/1.1.0/gc563 --with-
mpfr=/projects/community/mpfr/4.0.1/gc563 --with-
gmp=/projects/community/6.1.2/gc563 --disable-multilib
```
Next, the 'make' utility will be used to compile your software. It requires the file named 'Makefile' that you created with the configure step. The Makefile indicates the sequence that an options to be used for building various components of your new software. The Makefile uses labels (i.e., names for different sections of the procedures contained therein), so entire sections of the Makefile can be skipped or used in a particular order. Running 'make' to compile your code can take a long time for some packages. You can parallelize this step to some extent using the '-j [n]' option where n is the number of tasks you wish to run simultaneously. This isn't quite the same as running a parallel or multithreaded program, but it can help get a large number of compiling tasks done in a shorter time.  
`make -j 8`  

Finally, I need to copy my newly created executables and/or libraries to their final destinations (the location I specified with '--prefix=' in the configure step, above). One of the labels present in the Makefile is named 'install' and I can instruct make to run the commands under that label as follows:
```
 make install
 cd .. ; rm -rf gcc-7.3.0*
```

## Modulefiles for Community-Contributed Packages

All contributed software packages MUST have an associated [Lmod modulefile](http://lmod.readthedocs.io/en/latest) so users can easily access the software.

__Directory structure for modulefiles (omit parts that aren't used, like /CUDA/version or /mpi-name/version):__
`/projects/community/modulefiles/compiler-name/version/mpi-name/version/CUDA/version/pkg-name/version/NetID`

Here's an example Lmod modulefile for the GCC-7.3.0 example (above):
```
help(
[[
This module loads the GNU Compiler Collection version 7.3.0. 
The GNU Compiler Collection includes front ends for C, C++, Fortran, as well as libraries for these languages.
]])

whatis("Description: GCC: the GNU Compiler Collection")
whatis("URL: https://gcc.gnu.org")

conflict("gcc")

load("gmp/6.1.2-gc563")
load("mpfr/4.0.1-gc563")
load("mpc/1.1.0-gc563")

local base = pathJoin("/projects/community", myModuleName(), "7.3.0", "gc563")
prepend_path("PATH", pathJoin(base, "bin"))
prepend_path("C_INCLUDE_PATH", pathJoin(base, "include"))
prepend_path("CPLUS_INCLUDE_PATH", pathJoin(base, "include"))
prepend_path("LIBRARY_PATH", pathJoin(base, "lib64"))
prepend_path("LD_LIBRARY_PATH", pathJoin(base, "lib64"))
prepend_path("MANPATH", pathJoin(base, "share/man"))
```
Once created, this file should be named "7.3.0-gc563.lua" because that's the name that will appear in the list of modules when a user runs the 'module avail' or 'module spider' command.

__Understanding this file:__

A 'conflict' can be specified to prevent loading a potentially conflicting module (e.g., loading 2 different versions of GCC at the same time)

In this example, 3 prerequisite modules are automatically loaded when this GCC module is loaded. Alternatively, you can simply specify a prerequisite using 'prereq("module/version")' to notify a user of a prerequisite without automatically trying to load it.

The 'local base' statement establishes the general path for your software's various subdirectories

The 'prepend_path' statements define the specific path additions for those subdirectories











