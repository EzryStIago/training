
# Transition from sirius

On Feb 25, 2019, sirius machine was retired. This section documents information for users who are transitionining from sirius to Amarel. 

## Similarities

- most of the same software is available as on the sirius machine - and many more. If there is a special program you use, and it is not available, it may be installed upon request. 
- to see what software is there, try `module use /projects/community/modulefiles` and `module avail` (as on sirius). A comparison is provided in the last section
- similar graphical interface based on the Linux desktop is available as the interface to the cluster

## Differences

- sirius3 is not a single isolated machine - it is one "compute" node on the computer cluster of some 400+ machines
- there is a scheduler called slurm, which makes sharing of the resources fairer between users - making it impossible to monopolize resources
- sirius3 is intended to be used for graphical applications and for interactive code development. As such, code for more intense computations may be developed and tested on small scale on this machine. However, to perform a more intense computation, you should prepare a "batch" script (see cluster guide and video walktroughs for OpenOnDemand job templates). To ensure that there are enough resources for many users at once, each user on sirius3 is **limited to 2 cores and 4Gb of memory for a maximum of 24 hours**. 
- you can submit sbatch jobs, from the sirius3 desktop, to the wider Amarel cluster

## Instructions to connect to OpenOnDemand - Sirius3 Desktop

1. Open your browser at https://ondemand.hpc.rutgers.edu/  (make sure you are connected through VPN)
2. Go to Interactive Apps -> Sirius3 Deskop
3. Choose number of hours that your session will be active (max 24) and number of cores (max 2) and click Launch
4. In maximum 30 seconds, a blue button "Launch VNC in New Tab" will appear; clicking it will open Linux Desktop on sirius3

## Video walkthroughs

We have prepared some videos to illustrate the functionalities of OpenOnDemand graphical interface to the Amarel cluster: 

- TODO:  lauching Sirius3 desktop
- https://youtu.be/GPpDG1eAiKo - 4:24 mins - Amarel ondemand intro - available functionalities
- https://youtu.be/zNrpS5CUa-8 - 4:24 mins - launching Amarel desktop
- https://youtu.be/znmnbli0d5Y - 6:33 mins - using the Job Composer to keep track of slurm jobs, create new slurm jobs from templates, and create new templates.

## FastX users

**IMPORTANT NOTE**: When you login through FastX on Amarel, you are on the "login" node. This is a shared machine between tens of users at a time. It is essential that you don't run any intensive computations on this machine (like you did on sirius). Instead, you need to either ssh to sirius3, or submit a computational job through slurm (the scheduler). 

TODO: rest - how to connect from login node to sirius3. 

## Software comparison between sirius and Amarel

Please let us know if any of this software is not available on Amarel, but is essential to your work. We will consider adjusting our priorities to your needs: 

### Software that exists on Amarel: 

- AmberTools
- Bcl2fastq    
- Bcl2fastq2   
- Bedtools2  
- Bowtie2  
- Cufflinks  
- Clc-workbench
- Eclipse 
- FastQC   
- IGV      
- Ncbi-blast  - blast/2.6.0 on Amarel. This version 2.6, it is called ncbi-blast, and it should be called ncbi-blast 
- gdb
- Glibc  
- intel
- Java
- MATLAB
- Moe
- openbabel
- python
- RStudio  - available through OnDemand; R available on command line
- Rna-star    - this is called STAR/2.5.2a on Amarel
- Samtools
- Trimmomatic 
- Trinityrnaseq 
- Trinotate   
- Tophat 

### In question, but important software and is a work in progress: 

- Goldsuite - Vlad in communication with CCDC vendor
- SAS - licence will be ported from existing sirius machine

### Will install on request: 

- Weka 
- Knime
- Coot  -- crystalographic package
- Phenix  - crystallographic packag
- Mgltools  -  autodock related
- Autodock - community
- Autodock_vina -- community
- Chimera -- installed in community, needs module
- Vmd --  community

### Installation will progress as time permits 

- Blasr   
- Emboss  
- HMMER   
- Ngsplot 
- Qiime  
- PhyloCSF 
- PLINK   
- SNP2HLA 
- RSEM  

### These will not be considered for installation unless there are objections. 

- CHANCE
- MACS
- MACS2
- MEME
- MISO
- pandaseq
- abyss
- Transabyss
- WaveCNV

### Definitely not installing: 

- Sybyl-x
- FASTX-Toolkit
- Llvm
- Spartan
- binutils

