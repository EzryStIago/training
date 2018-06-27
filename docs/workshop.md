
# Genomic Software

This is a list of software to install for the workshop

|software|description|link|
|------|-----------------------------|-----------------------|
|Seqtk|very handy and fast for processing fastq/a files|[link](https://github.com/lh3/seqtk)|
|sratoolkit|downloading and processing data from GEO/SRA database|[link](http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-centos_linux64.tar.gz)|
|htseq-count|counting the  reads mapped on to genomics feature|[link](http://htseq.readthedocs.io/en/master/install.html)|
|fastQC|widely used for sequencing read QC|[link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.6.zip)|
|RSeQC-2.6.4|An RNA-seq quality control package, multiple functions|python package|
|trimmomatic|fastq quality trim and adaptor removal|[link](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)|

This is a list of software already available on the cluster and the command you need to execute to load it in your environment: 

|software|description|load it on the cluster|
|------|-----------------------------|-----------------------|
|Samtools||`module load samtools`|
|Bedtools||`module load bedtools2./2.25.0`|
|bowtie2| alignment software|`module load bowtie2`|
|tophat2| alignment software|`module load mvapich2/2.1  boost/1.59.0  tophat2/2.1.0`|
|R|language for statistical analysis|`module load intel/17.0.4 R-Project/3.4.1`|

This is a list of other software you might find useful
               | GSEA     | genome set enrichment analysis| [link] (http://software.broadinstitute.org/gsea/index.jsp) |
               | IGV      | Interactive Genome Viewer      | [link] (http://software.broadinstitute.org/software/igv/)|
               | Cytoscape| Network visualization softwar   | [link] (http://www.cytoscape.org/)|

# Setup

1. Connect to the cluster login node by one of the following methods: 
  - *via a terminal*: if you have a Mac or Linux, terminal is part of your standard apps. If you have Windows, install an SSH client such as `putty` or `moba-xterm`. Then from your terminal connect to the cluster by executing the following command: 
``` ssh -X <your net id>@amarel.hpc.rutgers.edu```
   This is the preferred way, as your copy-pasting will most likely work best. 
  - *via FastX*: in your browser, go to `https://amarel.hpc.rutgers.edu:3443`

2. Get resources on the compute node of the cluster to execute your computations by running the following command in your terminal: 
```srun  -p main --reservation=genomics -N 1 -c 2 -n 1 -t 01:40:00 --export=ALL --pty /bin/bash```
Notice that the name in your terminal will change from `amarel` to node name like `hal0025` or `slepner086`. This means that you will not impede other users who are also using the login node, and will be placed on a machine which you share with only a few people. This explains the parts of this command: 

|command part| meaning|
|----|----|
|`srun`| `slurm` run, i.e. allocate resources and run via `slurm` scheduler | 
|`-p main` | on the main partition, one of several queues on the cluster|
|`--reservation=genomics`| we reserved some compute nodes for this workshop to not wait long for resources|
|`-N 1`| ask for one node|
|`-c 2`| ask for two cores|
|`-n 1`| this will be 1 most times|
|`-t 01:40:00`| run this for a maximum time of 1 hour 40 minutes|
|`--pty /bin/bash`| run the terminal shell in an interactive mode|

3. *Prepare some directories for the data* 
   You have two main spaces on the Amarel cluster. These are `/home/netid/` (e.g. `/home/kp807/` for my netid) and '/scratch/netid/'. They differ in how often they are backed up and by size (100Gb for `/home` and 500Gb for `/scratch`). So we will install programs in `/home`, while the data and output will be in `/scratch`. Execute these commands: 
```
                cd ~                      # change directory to your home directory
                mkdir Genomics_Workshop
                cd Genomics_Workshop
                mkdir Programs            # download and install programs here
```

```
                mkdir -p  /scratch/$USER/Genomics_Workshop/
                cd /scratch/$USER/Genomics_Workshop/
                mkdir untreated  
                mkdir dex_treated
```

4. Install programs
   Each program will have slightly different installation instructions. Here is a handy sequence of commands that will install them: 
```
##We are going to do some modifications to a system file .bashrc, be careful doing it and make sure that you created a copy of your .bashrc file
                cd
                cp .bashrc .bashrc_20180118
                nano .bashrc
##At the end of the file add the line  “##  Genomics_Workshop 06/27/2018 settings”
##Exit nano (ctrl+x)

#Seqtk:   https://github.com/lh3/seqtk   ##  very handy and fast for processing fastq/a files
                cd
                cd Genomics_Workshop/Programs/
                git clone https://github.com/lh3/seqtk.git 
                cd seqtk
                make
                echo ‘export PATH=$HOME/Genomics_Workshop/Programs/seqtk:$PATH’ >>  ~/.bashrc
                source ~/.bashrc

#sratoolkit     https://www.ncbi.nlm.nih.gov/books/NBK158900/
                        https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
##for downloading and processing data from GEO/SRA database
                cd
                cd Genomics_Workshop/Programs/
                wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-centos_linux64.tar.gz
                tar zxvf sratoolkit.2.8.2-centos_linux64.tar.gz
                echo ‘export PATH=$HOME/Genomics_Workshop/Programs/sratoolkit.2.8.2-centos_linux64/bin:$PATH’ >> ~/.bashrc
                source ~/.bashrc 

#htseq-count    http://htseq.readthedocs.io/en/master/install.html ##for counting the  reads mapped on to genomics feature

                cd
                module load intel/17.0.2 python/2.7.12
                pip install HTSeq --user

#fastQC          #widely used for sequencing read QC

                cd
                cd  Genomics_Workshop/Programs/
                wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.6.zip
                unzip  fastqc_v0.11.6.zip
                echo ‘export PATH=$HOME/Genomics_Workshop/Programs/FastQC:$PATH’ >> ~/.bashrc
                source ~/.bashrc


#FASTX-toolkit    http://hannonlab.cshl.edu/fastx_toolkit/   (##also a tool kit for fastq processing, quality trim, adaptor removal, etc. try if time allows)

#RSeQC-2.6.4     ##An RNA-seq quality control package, multiple functions

                cd
                cd  Genomics_Workshop/Programs/
                module load python/2.7.12
                module load intel/17.0.2
                pip install RSeQC --user


#trimmomatic             ##for fastq quality trim and adaptor removal
              
                cd
                cd  Genomics_Workshop/Programs/
                wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
                unzip Trimmomatic-0.36.zip 
                java -jar ~/Genomics_Workshop/Programs/Trimmomatic-0.36/trimmomatic-0.36.jar #-h                    

```
Moreover, execute these commands to load system-installed sofware so the system knows where to find it (i.e. `samtools command will work`): 
```
module load samtools       
module load bedtools2./2.25.0
module load bowtie2
module load mvapich2/2.1  boost/1.59.0  tophat2/2.1.0
module load intel/17.0.4 R-Project/3.4.1
```

# Download data

#
