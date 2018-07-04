# Using R

Generally, there are 2 approaches for accessing R on Amarel:  
1. use one of the pre-installed R modules named R-Project/*version* (these versions come bundled with a very broad range of common and useful tools).   
2. install your own custom build of R in your /home directory or in a shared directory (e.g. `/projects/[group]` or `/projects/community`).

##Using pre-installed R modules

Start by finding which module you wish to use with the 'module spider R-Project' command:
```
module spider R-Project
--------------------------------------------------
  R-Project:
--------------------------------------------------
    Description:
      R: The R Project for Statistical Computing
     Versions:
        R-Project/3.2.2
        R-Project/3.2.5
        R-Project/3.3.3
        R-Project/3.4.1
--------------------------------------------------
  To find detailed information about R-Project please enter the full name.
  For example:
     $ module spider R-Project/3.4.1
--------------------------------------------------
```
Next, use 'module spider' again to see how to load the module you wish to use (e.g., are there any prerequisites that must be loaded first?):
```
 module spider R-Project/3.4.1
--------------------------------------------------
  R-Project: R-Project/3.4.1
--------------------------------------------------
    Description:
      R: The R Project for Statistical Computing
    This module can only be loaded through the following modules:
      intel/17.0.4
    Help:    
      This module loads the installation R-Project 3.4.1 compiled with the Intel 17.0.4 compilers.
```
Load the R-Project module of your choice:
```
module load intel/17.0.4 R-Project/3.4.1
which R
/opt/sw/packages/intel-17.0.4/R-Project/3.4.1/bin/R
```
What R packages are already installed?
```R
pkgs <- installed.packages ()
pkgs[,c("Package", "Version")]

                     Package                Version    
base                 "base"                 "3.4.4"    
BH                   "BH"                   "1.66.0-1" 
Biobase              "Biobase"              "2.38.0"   
BiocGenerics         "BiocGenerics"         "0.24.0"   
BiocInstaller        "BiocInstaller"        "1.28.0"   
BiocParallel         "BiocParallel"         "1.12.0"   
Biostrings           "Biostrings"           "2.46.0"   
bitops               "bitops"               "1.0-6"    
boot                 "boot"                 "1.3-20"   
BSgenome             "BSgenome"             "1.46.0"
...
```
It's very common for uers to need additional or custom packages for a base R installation. On a large, shared computing system, users are unable to install (write) to the usual places where R places new packages by default (/usr/local or /usr/lib). Therefore, managing your own local package/library installation location is necessary. In the example below, I'll demonstrate how I did this for my Amarel user account.

First, I'll create a directory where I can store my locally-installed R packages. This can have any name and it can be located anywhere you have access:  
`mkdir ~/my.R.libs`

Next, to ensure that my new private R packages directory is searched when I try to load a library that's installed there, I need to make an environment setting that will point R to the right location. I'll create a new file in my /home directory named .Renviron (note the leading "." in that name) and I'll add the following line to that file:  
`export R_LIBS=~/my.R.libs`  
Now, every time I start any version of R, my ~/my.R.libs directory will be the first location to be searched when loading a library.

** Some important notes about installing packages:**

There are a variety of different ways to install packages in R. The most straightforward way is to use the built-in 'install.packages()' function while R is running. Using this approach gives you the flexibility to install the latest version of a package or you can specify an older version of a package. To install a specifc version of a package, you'll need the URL (web address) for the tarball (*.tar.gz or *.tgz file) containing the source code for that version.

For example, I want to load the following list of packages, and I need the specifc versions listed here: 'kernlab' version 0.9-24 'ROCR' version 1.0.7 'class' version 7.3.14 'party' version 1.0.25 'e1071' version 1.6.7 'randomForest' version 4.6.12 I can use a web search to find the source tarballs for these packages. For example, to find 'kernlab' version 0.9-24, I search for "kernlab" and find the [website](https://cran.r-project.org/web/packages/kernlab/index.html). At that site, I see that 0.9-25 is the current version (not what I want), but there is "kernlab archive" link there that takes me to a long list of previous versions. I see a link for version 0.9-24 there, so I copy that URL and use that URL in my install.packages() command:  
```R
install.packages("https://cran.r-project.org/src/contrib/Archive/kernlab/kernlab_0.9-24.tar.gz", lib="~/my.R.libs")
```

The other packages I need can be found in the same way. While installing them, I find that 'ROCR-1.0.7' requires 'gplots' and 'party-1.0-25' requires 6 other prerequisites. So, I have to also install those prerequisite packages. In the end, my install.packages() commands are as follows:
```R
install.packages("gplots", lib="~/my.R.libs")
install.packages("https://cran.r-project.org/src/contrib/ROCR_1.0-7.tar.gz", lib="~/my.R.libs")
install.packages("https://cran.r-project.org/src/contrib/class_7.3-14.tar.gz", lib="~/my.R.libs")
install.packages(c("mvtnorm","modeltools","strucchange","coin","zoo","sandwich"), lib="~/my.R.libs")
install.packages("https://cran.r-project.org/src/contrib/Archive/party/party_1.0-25.tar.gz", lib="~/my.R.libs")
install.packages("https://cran.r-project.org/src/contrib/Archive/e1071/e1071_1.6-7.tar.gz", lib="~/my.R.libs")
install.packages("https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz", lib="~/my.R.libs")
```
Once all of my package installs have completed successfully, those packages can be loaded normally and they will be available every time I log-in to the cluster.


## Installing your own build of R

For some users or groups, installing and customizing or even modifying the latest version (or a specific version) of R is necessary. For those users, I'll demonstrate how to install a version of R below.

Here are the commands to use for installing R-3.4.4 from source:
```
wget https://cran.r-project.org/src/base/R-3/R-3.4.4.tar.gz
tar -zxf R-3.4.4.tar.gz
cd R-3.4.4
module load gcc/5.4 java/1.8.0_161
./configure --prefix=/home/gc563/R/3.4.4 --enable-java
make -j 4
make install
cd ..
rm -rf R-3.4.4*
```
Here, I have loaded the GCC compiler suite and Java before installing R. This is an optional step and I did it because there might be specific tools I will use with R that require these extra software packages.

When I configured my installation, I specified (with --prefix=) that I want R to be installed in my /home directory. I prefer to use a [package]/[version] structure because that enables easy organization of multiple verisons of any software package. It's a personal preference.

At the end of my install procedure, I remove the downloaded install package and tarball, just to tidy-up.


Since I've installed R in my /home directory, I can add packages using the default library directory since that too will be in my /home directory.


Before using my new R installation, I'll need to set some environment variables and load needed modules (the same ones I used for building my R installation). This can be done from the command line (but the settings won't persist after you log-out) or by adding these commands to the bottom of your ~/.bashrc file (so the settings will persist):

```
module load gcc/5.4 java/1.8.0_161
export PATH=/home/gc563/R/3.4.4/bin:$PATH
export LIBRARY_PATH=/home/gc563/R/3.4.4/lib64
export LD_LIBRARY_PATH=/home/gc563/R/3.4.4/lib64
export MANPATH=/home/gc563/R/3.4.4/share/man
```
If you're adding these lines to the bottom of your ~/.bashrc file, log-out and log-in again, then verify that the settings are working:

```
$ module list
Currently Loaded Modules:
  1) gcc/5.4   2) java/1.8.0_161
$ which R
~/R/3.4.4/bin/R
```
Now that my new R installation is setup, I can begin adding R packages. Since this is my own installation of R and not one of the preinstalled versions available on the cluster, my default packages/libraries directory is /home/gc563/R/3.4.4/lib64/R/library  
```
> .libPaths()  
[1] "/home/gc563/R/3.4.4/lib64/R/library"
```

Install a package:
```R
install.packages("rJava")
library(rJava)
```
**Saving figures/ plots from R (without a display):**  
Need to save a PDF, PostScript, SVG, PNG, JPG, or TIFF file in your working directory? Normally, writing a graphics file from R requires a display of some kind and the X11 protocol. That's often not convenient for batch jobs running on the cluster. Alternatively, you can use the Cairo graphics device/library for R. Cairo enables you to write bitmap or vector graphics directly to a file. Here's an example:
```R
$ R --no-save
png('my-figure.png', type='cairo')
plot(rnorm(10),rnorm(10))
dev.off()
q()
```
# Using Python

Generally, there are 2 approaches for using Python and its associated tools: (1) use one of the pre-installed Python modules (version 2.7.x or 3.5.x) which come bundled with a very broad range of common and useful tools (you can add or update packages if needed) or (2) install your own custom build of Python in your /home directory or in a shared directory (e.g., /projects/[group] or /projects/community).

## Using pre-installed Python modules 

With the pre-installed Python modules, you can add or update Python modules/packages as needed if you do it using the '--user' option for pip. This option will instruct pip to install new software or upgrades in your ~/.local directory. Here's an example where I'm installing the Django package:
```
module load python/3.5.2
pip install --user Django
```
Note: if necessary, pip can also be upgraded when using a system-installed build of Python, but be aware that the upgraded version of pip will be installed in ~/.local/bin. Whenever a system-installed Pytyon module is loaded, the PATH location of that module's executables (like pip) will precede your ~/.local/bin directory. To run the upgraded version of pip, you'll need to specify its location because the previous version of pip will no longer work properly:
```
$ which pip
/opt/sw/packages/gcc-4.8/python/3.5.2/bin/pip
$ pip --version
pip 9.0.3 from /opt/sw/packages/gcc-4.8/python/3.5.2/lib/python3.5/site-packages (python 3.5)
$ pip install -U --user pip
Successfully installed pip-10.0.1
$ which pip
/opt/sw/packages/gcc-4.8/python/3.5.2/bin/pip
$ pip --version
Traceback (most recent call last):
  File "/opt/sw/packages/gcc-4.8/python/3.5.2/bin/pip", line 7, in 
    from pip import main
ImportError: cannot import name 'main'

$ .local/bin/pip --version
pip 10.0.1 from /home/gc563/.local/lib/python3.5/site-packages/pip (python 3.5)
$ .local/bin/pip install --user Django
```

## Building your own Python installation

Using this approach, I must specify that I want Python to be installed in my /home directory. This is done using the '--prefix=' option. Also, I prefer to use a [package]/[version] naming scheme because that enables easy organization of multiple verisons of Python (optional, it's just a personal preference).  
At the end of my install procedure, I remove the downloaded install package and tarball, just to tidy-up.

```
wget https://www.python.org/ftp/python/3.6.5/Python-3.6.5.tgz
tar -zxf Python-3.6.5.tgz
cd Python-3.6.5
./configure --prefix=/home/gc563/python/3.6.5
make -j 4
make install
cd ..
rm -rf Python-3.6.5*
```
Before using my new Python installation, I'll need to set or edit some environment variables. This can be done from the command line (but the settings won't persist after you log-out) or by adding these commands to the bottom of your ~/.bashrc file (so the settings will persist):
```
export PATH=/home/gc563/python/3.6.5/bin:$PATH
export LD_LIBRARY_PATH=/home/gc563/python/3.6.5/lib
export MANPATH=/home/gc563/python/3.6.5/share/man
```
If you're adding these lines to the bottom of your ~/.bashrc file, log-out and log-in again, then verify that the settings are working:  
`which python3  
~/python/3.6.5/bin`


# Running GROMACS 

Here is a simple example procedure that demonstrates how to use GROMACS 2016 on Amarel. In this example, we’ll start with a downloaded PDB file and proceed through importing that file into GROMACS, solvating the protein, a quick energy minimization, and then an MD equilibration. This example is not intended to teach anyone how to use GROMACS. Instead, it is intended to assist new GROMACS users in learning to use GROMACS on Amarel.  

1. Download a PDB file.  
    ```wget https://files.rcsb.org/view/5EWT.pdb```  
2. Load the GROMACS software module plus any needed prerequisites.  
    ```
    module purge  
    module load intel/17.0.1 mvapich2/2.2 gromacs/2016.1
    ```
3. Import the PDB into GROMACS, while defining the force field and water model to be used for this system.  
   ```gmx_mpi pdb2gmx -f 5EWT.pdb -ff charmm27 -water tip3p -ignh -o 5EWT.gro -p 5EWT.top -i 5EWT.itp``` 
4. Increase the size of the unit cell to accommodate a reasonable volume of solvent around the protein.  
    ```gmx_mpi editconf -f 5EWT.gro -o 5EWT_newbox.gro -box 10 10 10 -center 5 5 5```
5.  Now add water molecules into the empty space in the unit cell to solvate the protein.  
    ```gmx_mpi solvate -cp 5EWT_newbox.gro -p 5EWT.top -o 5EWT_solv.gro```
6. Prepare your SLURM job script(s). The 2 **mdrun** commands in the following steps can be executed from within an interactive session or they can be run in batch mode using job scripts. If your mdrun commands/job might take more than a few minutes to run, it would be best to run them in batch mode using a job script. Here’s an example job script for a GROMACS MD simulation. To run the 2 **mdrun** commands below, simply replace the example **mdrun** command in this script with one of the mdrun commands from the steps below and submit that job after preparing the simulation with the appropriate **grompp** step.
```
#!/bin/bash
#SBATCH --partition=main                # Partition (job queue)
#SBATCH --job-name=gmdrun               # Assign an 8-character name to your job
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=16                     # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1               # Threads per process (or per core)
#SBATCH --mem=124000                    # Memory per node (MB)
#SBATCH --time=00:20:00                 # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out        # combined STDOUT and STDERR output file
#SBATCH --export=ALL                    # Export you current env to the job env

module purge
module load intel/17.0.1 mvapich2/2.2 gromacs/2016.1
srun --mpi=pmi2 gmx_mpi mdrun -v -s 5EWT_solv_prod.tpr \
                -o 5EWT_solv_prod.trr -c 5EWT_solv_prod.gro \
                -e 5EWT_solv_prod.edr -g 5EWT_solv_prod.md.log
```
7. Perform an inital, quick energy minimization. Here, we’re using a customized MD parameters file named em.mdp, which contains these instructions:
```
integrator     = steep
nsteps         = 200
cutoff-scheme  = Verlet
coulombtype    = PME
pbc            = xyz
emtol          = 100
```  
These are the commands (both the grompp step and the mdrun step) used to prepare and run the minimization:  
```
gmx_mpi grompp -f em.mdp -c 5EWT_solv.gro -p 5EWT.top -o 5EWT_solv_mini.tpr -po 5EWT_solv_mini.mdp

gmx_mpi mdrun -v -s 5EWT_solv_mini.tpr -o 5EWT_solv_mini.trr -c 5EWT_solv_mini.gro -e 5EWT_solv_mini.edr -g 5EWT_solv_mini.md.log
```
8. Perform a quick MD equilibration (same syntax/commands for a regular MD run). Here, we’re using a customized MD parameters file named equil.mdp, which contains these instructions:
```
integrator               = md
dt                       = 0.002
nsteps                   = 5000
nstlog                   = 50
nstenergy                = 50
nstxout                  = 50
continuation             = yes
constraints              = all-bonds
constraint-algorithm     = lincs
cutoff-scheme            = Verlet
coulombtype              = PME
rcoulomb                 = 1.0
vdwtype                  = Cut-off
rvdw                     = 1.0
DispCorr                 = EnerPres
tcoupl                   = V-rescale
tc-grps                  = Protein  SOL
tau-t                    = 0.1      0.1
ref-t                    = 300      300
pcoupl                   = Parrinello-Rahman
tau-p                    = 2.0
compressibility          = 4.5e-5
ref-p                    = 1.0
```
These are the commands (both the grompp step and the mdrun step) used to prepare and run the equilibration:
```
gmx_mpi grompp -f equil.mdp -c 5EWT_solv_mini.gro -p 5EWT.top -o 5EWT_solv_equil.tpr -po 5EWT_solv_equil.mdp

gmx_mpi mdrun -v -s 5EWT_solv_equil.tpr -o 5EWT_solv_equil.trr -c 5EWT_solv_equil.gro -e 5EWT_solv_equil.edr -g 5EWT_solv_equil.md.log
```
# Running TensorFlow with a GPU

TensorFlow has two versions of its python package: `tensorflow` and `tensorflow-gpu`, but confusingly the command to use it is the same in both cases: `import tensorflow as tf`  (and not `import tensorflow-gpu as tf` in case of the GPU version). This means that it really matters which package is installed in your environment. 

- you can control your environment using Singularity image (but the problem arises if you need a package not included in the prebuilt image, in which case you need to build the image yourself)
- you can control your environment using conda environments (or virtual-env). 

## Using Singularity

To do this, you can use the [Singularity](http://singularity.lbl.gov/) container manager and a Docker image containing the TensorFlow software. 

Running Singularity can be done in batch mode using a job script. Below is an example job script for this purpose (for this example, we'll name this script `TF_gpu.sh`)
```
#!/bin/bash
#SBATCH --partition=main             # Partition (job queue)
#SBATCH --no-requeue                 # Do not re-run job  if preempted
#SBATCH --job-name=TF_gpu            # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=2            # Cores per task (>1 if multithread tasks)
#SBATCH --gres=gpu:1                 # Number of GPUs
#SBATCH --mem=16000                  # Real memory (RAM) required (MB)
#SBATCH --time=00:30:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional)
#SBATCH --export=ALL                 # Export you current env to the job env

module purge
module load singularity/.2.5.1

srun singularity exec --nv docker://tensorflow/tensorflow:1.4.1-gpu python 
```
Once your job script is ready, submit it using the sbatch command:  
```$ sbatch TF_gpu.sh```
Alternatively, you can run Singularity interactively:
```
$ srun --pty -p main --gres=gpu:1 --time=15:00 --mem=6G singularity shell --nv docker://tensorflow/tensorflow:1.4.1-gpu

Docker image path: index.docker.io/tensorflow/tensorflow:1.4.1-gpu
Cache folder set to /home/user/.singularity/docker
Creating container runtime...
Importing: /home/user/.singularity/docker/sha256:054be6183d067af1af06196d7123f7dd0b67f7157a0959bd857ad73046c3be9a.tar.gz
Importing: /home/user/.singularity/docker/sha256:779578d7ea6e8cc3934791724d28c56bbfc8b1a99e26236e7bf53350ed839b98.tar.gz
Importing: /home/user/.singularity/docker/sha256:82315138c8bd2f784643520005a8974552aaeaaf9ce365faea4e50554cf1bb44.tar.gz
Importing: /home/user/.singularity/docker/sha256:88dc0000f5c4a5feee72bae2c1998412a4b06a36099da354f4f97bdc8f48d0ed.tar.gz
Importing: /home/user/.singularity/docker/sha256:79f59e52a355a539af4e15ec0241dffaee6400ce5de828b372d06c625285fd77.tar.gz
Importing: /home/user/.singularity/docker/sha256:ecc723991ca554289282618d4e422a29fa96bd2c57d8d9ef16508a549f108316.tar.gz
Importing: /home/user/.singularity/docker/sha256:d0e0931cb377863a3dbadd0328a1f637387057321adecce2c47c2d54affc30f2.tar.gz
Importing: /home/user/.singularity/docker/sha256:f7899094c6d8f09b5ac7735b109d7538f5214f1f98d7ded5756ee1cff6aa23dd.tar.gz
Importing: /home/user/.singularity/docker/sha256:ecba77e23ded968b9b2bed496185bfa29f46c6d85b5ea68e23a54a505acb81a3.tar.gz
Importing: /home/user/.singularity/docker/sha256:037240df6b3d47778a353e74703c6ecddbcca4d4d7198eda77f2024f97fc8c3d.tar.gz
Importing: /home/user/.singularity/docker/sha256:b1330cb3fb4a5fe93317aa70df2d6b98ac3ec1d143d20030c32f56fc49b013a8.tar.gz
Importing: /home/user/.singularity/metadata/sha256:b71a53c1f358230f98f25b41ec62ad5c4ba0b9d986bbb4fb15211f24c386780f.tar.gz
Singularity: Invoking an interactive shell within container...

Singularity tensorflow:latest-gpu:~> 
```
Now, you're ready to execute commands:

```

Singularity tensorflow:latest-gpu:~> python -V
Python 2.7.12
Singularity tensorflow:latest-gpu:~> python3 -V
Python 3.5.2
```
Please remember to exit from your interactive job after you are finished with your calculations.  

There are several Docker images available on Amarel for use with Singularity. The one used in the example above, tensorflow:1.4.1-gpu, is intended for python 2.7.12. If you want to use Python3, you'll need a different image, docker://tensorflow/tensorflow:1.4.1-gpu-py3, and the Python command will be `python3` instead of `python` in your script.

## Using conda

You can either install your own version of Anaconda in your home directory, or you can use a community module. 
```
module use /projects/community/modulefiles  #loads community-contributed software packages
module keyword anaconda                     #search packages with 'anaconda' in description 
```
output: 
```
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
The following modules match your search criteria: "anaconda"
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  py-data-science-stack: py-data-science-stack/5.1.0-kp807
    Sets up anaconda 5.1.0 in your environment
```
So here are commands with which you can load tensorflow package: 
```
module load py-data-science-stack/5.1.0-kp807
conda env list                                 # be patient
source activate tensorflow-gpu-1.7.0           #
```

Note that if you try to import tensorflow on a node without a gpu, the import will fail, because it will try to load cuda driver that is not installed (because there is no GPU on the machine). 
So you need to include this line in a slurm script, where you request a GPU resource in slurm: 
```#SBATCH --gres=gpu:1   ```












 

















