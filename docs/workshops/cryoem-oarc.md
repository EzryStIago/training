## As a part of ongoing collaboration between Rutgers [**Office of Advanced Research Computing (OARC)**](http://oarc.rutgers.edu) and Rutgers [**Institute for Quantitative Biomedicine**](http://iqb.rutgers.edu/) we are proud to announce that OARC clusters were chosen to be a primary computational platform to support the 6th Annual Interdisciplinary Quantitative Biology Winter Boot Camp for  Single-Particle Cryo-Electron Microscopy January 7 - 11, 2019 https://iqb.rutgers.edu/bootcampwinter2019

### We are proud of being at the frontier of science and helping to spread the knowledge of Rutgers scientists with the  national and international scientific community.

# Setup

## Connect to the OARC cluster 
The boot camp participants are expected to bring their own laptops with a **Chrome web browser** installed.
The preferred method to connect to the cluster is through remote desktop (fastX) in a web browser.
Each participant will be granted an access to an individual compute node for the duration of the boot camp. 
The compute node will be accessible only on the Rutgers network or if you are connected via [VPN](https://ssl-vpn.rutgers.edu/dana-na/auth/url_default/welcome.cgi). The access port is unique for each participant and should not be shared.


 **Remote desktop via FastX**: 
 In your browser, go to https://cryoem-oarc.hpc.rutgers.edu:port <p>
 Enter your Rutgers netID and Password and click on the **Log in** button.
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/fastX_login.png" width="500" height="500" center> <p>
 
 In the new window, click on **Launch Session** on the left, select **XFCE desktop** and then click on **Launch** button on the right.<p>
<img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_sessions.png"> 
  
<img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_XFCE.png"> 
<p> 
 Ignore the warning<p> 
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_popup.png" width="600" height="150" center> 
<p>
 
 **Reconnect to the FastX remote desktop session**:
 - Open the web browser and go to https://cryoem-oarc.hpc.rutgers.edu:port <br>
 - Login into your account with your Rutgers netID and a password.<br>
 - Click on the **XFCE Desktop** session in the My Session list. <br>
 
  
# Overview





## CRYO-EM Software

This is a list of software to install for the workshop

|Software|Description|URL|
|------|-----------------------------|-----------------------|
|Scipion|An image processing framework for obtaining 3D models of macromolecular complexes using Electron Microscopy (3DEM)|[[link]](https://github.com/I2PC/scipion/wiki)|
|Scipion packages|Scipion integrated packages: Bsoft, cryoFF, CTFFIND, Chimera, EMAN, ETHAN, Frealign, Gautomatch, gCTF, gEMpicker, IMAGIC-4D, Localrec, Motioncorr/dosefgpu, Motioncorr2, Relion, ResMap, SIMPLE, SPIDER, Xmipp, Unblur & summovie|[[link]](https://github.com/I2counting the  reads mapped on to genomics featurePC/scipion/wiki/Integrated-Packages)|
|Chimera|A program for interactive visualization and analysis of molecular structures and related data, including density maps|[[link]](https://www.cgl.ucsf.edu/chimera/)|
|fastQC|widely used for sequencing read QC|[[link]](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) [[download]](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.6.zip)|
|RSeQC-2.6.4|An RNA-seq quality control package, multiple functions|[[link]](http://rseqc.sourceforge.net/)|
|trimmomatic|fastq quality trim and adaptor removal|[[link]](http://www.usadellab.org/cms/?page=trimmomatic) [[download]](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)|




## DO NOT RUN ANY COMPUTATIONAL JOBS ON THE LOGIN NODE 
## Get resources on the compute node 

When you login to the cluster you are on the login node. Jobs are not allowed to be run on the loging node, intstead you need to request a resource on the compute node for your job. This means that you will not impede other users who are also using the login node, and will be placed on a machine which you share with only a few people. You can do so by running the following command in your terminal:    
```
srun  -p main --x11 --reservation=genomics -N 1 -c 2 -n 1 -t 01:40:00 --pty /bin/bash -i
```  
 or just run this script
```
node_request.sh
```
Notice that the name in your terminal will change from `amarel` to a node name like `hal0025` or `slepner086`.  The following table explains the parts of this command: 

|command part| meaning|
|----|----|
|`srun`| `slurm` run, i.e. allocate resources and run via `slurm` scheduler | 
|`-p main` | on the main partition, one of several queues on the cluster|
|`--x11` | it allows the graphical output from a compute node, e.g. GUI of the program|
|`--reservation=genomics`| we reserved some compute nodes for this workshop to avoid waiting in the queue|
|`-N 1`| ask for one compute node|
|`-c 2`| ask for two cpu cores|
|`-n 1`| tells slurm that the job will be run as  1 task ( for parallel MPI jobs it could be more than 1|
|`-t 01:40:00`| run this for a maximum time of 1 hour 40 minutes|
|`--pty /bin/bash -i`| run the terminal shell in an interactive mode|

### Understanding your data space

You have two main spaces on the Amarel cluster. These are: 

- your home directory (100Gb) - `/home/netid/` 
- your scratch directory (500Gb)- `/scratch/netid/` 

  They differ in how they are backed up (scratch is not backed up) and by read/write speed. <br> 
  So we will install programs in  `/home/$USER/Genomics_Workshop`, <br>
  while the data and computational output will be held in `/scratch/$USER/Genomics_Workshop`. 
  

## 1. Install programs and create a workspace for the workshop

   Each program has slightly different installation instructions. 
   You do not need to install programs manually.  Instead just run the following scirpt:
   **/projects/oarc/Genomics_Workshop/RNA-Seq_analysis/Labs/lab_PartI.sh**<br>
   It will install neccessary programs and creates folders for this workshop.
   
   For curious one, here is the content of the script
   ########################################################
 ```  
 #!/bin/bash

mkdir -p /home/$USER/Genomics_Workshop/
mkdir -p /scratch/$USER/Genomics_Workshop/scripts
mkdir -p /scratch/$USER/Genomics_Workshop/download
mkdir -p /scratch/$USER/Genomics_Workshop/untreated
mkdir -p /scratch/$USER/Genomics_Workshop/dex_treated
mkdir -p /scratch/$USER/Genomocs_Workshop/Reference

echo "Copying files... Please wait"

cp  -r /projects/oarc/Genomics_Workshop/Programs/ /home/$USER/Genomics_Workshop/

echo '## Genomics_Workshop specific settings 07/16/2018' >> ~/.bashrc
echo 'export PATH=$HOME/Genomics_Workshop/Programs/seqtk:$PATH' >> ~/.bashrc
echo 'export PATH=$HOME/Genomics_Workshop/Programs/sratoolkit.2.8.2-centos_linux64/bin:$PATH' >> ~/.bashrc
echo 'export PATH=$HOME/Genomics_Workshop/Programs/FastQC:$PATH' >> ~/.bashrc
source ~/.bashrc 


module load intel/17.0.2 python/2.7.12
pip install HTSeq --user
wait
pip install bx-python==0.7.3 --user
wait
pip install RSeQC --user
```
   ########################################################<br>
   


