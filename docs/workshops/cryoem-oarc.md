## As a part of ongoing collaboration between Rutgers [**Office of Advanced Research Computing (OARC)**](http://oarc.rutgers.edu) and Rutgers **Institute for Quantitative Biomedicine** we are proud to announce that OARC clusters were chosen to be a primary computational platform to support the 6th Annual Interdisciplinary Quantitative Biology Winter Boot Camp for  Single-Particle Cryo-Electron Microscopy January 7 - 11, 2019 https://iqb.rutgers.edu/bootcampwinter2019

### We are proud of being at the frontier of science and helping to spread the knowledge of Rutgers scientists with the  national and international scientific community.

# Setup

## Connect to the OARC cluster 
The boot camp participants are expected to bring their own laptops with a **Chrome web browser** installed.
The preferred method to connect to the cluster is through remote desktop (fastX) in a web browser.
Each participant will be granted an access to an individual compute node for the duration of the boot camp. 
The compute node will be accessible only on the Rutgers network or if you are connected via [VPN](https://ssl-vpn.rutgers.edu/dana-na/auth/url_default/welcome.cgi). The access port is unique for each participant and should not be shared.


 **Remote desktop via FastX**: 
 - in your browser, go to https://cryoem-oarc.hpc.rutgers.edu:port
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/fastX_login.png" width="500" height="500" center> 
 Enter your Rutgers netID and Password and click on the Login button.
 
 - in the new window 
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_sessions.png"> 
  Click on Launch Session on the left, in a new window select XFCE desktop and then click on Launch button on the right.
  
- in the new window
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_XFCE.png"> 
  Click on Launch Session, then  Select XFCE desktop and in the new.
 
  - ignore the warning 
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_popup.png" width="500" height="300" center> 

 
 
 
 |Login page|Click on Launch Session and Select XFCE desktop|Open new terminal, copy text from your computer into clipboard window and click on Send to Remote. Put cursor into a terminal and click Edit->Paste. The text from a clipboard appears in the terminal. Hit Return to submit a job to a compute node. Notice that the jobid number has been assigned and resources have been allocated. The command prompt in the treminal will also change from a login node netid@amarel to a compute node netid@slepner or netid@hal. 
|------|-----------------------------|-----------------------|
| ![](https://github.com/rutgers-oarc/training/blob/master/workshops/fastX_web.png "FastX web interface")     | ![](https://github.com/rutgers-oarc/training/blob/master/workshops/FX_newsession.png "FastX new session")| ![](https://github.com/rutgers-oarc/training/blob/master/workshops/FX_srun.png "FastX terminal srun") |
 
On rare occasions, especially if the user has a modified .bashrc file,  FastX doesn't work. If it is the case, you can still use a terminal but the connection will be much slower for graphical outputs.
  - **via a terminal**: if you have a Mac or Linux, a terminal is part of your standard apps. If you have Windows, install an SSH client such as  [mobaXterm](https://mobaxterm.mobatek.net/). Then from your terminal connect to the cluster by executing the following command:   
``` ssh -X <your netid>@amarel.hpc.rutgers.edu```   
  
# Overview





## Genomic Software

This is a list of software to install for the workshop

|Software|Description|URL|
|------|-----------------------------|-----------------------|
|seqtk|very handy and fast for processing fastq/a files|[[link]](https://github.com/lh3/seqtk)|
|sratoolkit|downloading and processing data from GEO/SRA database|[[link]](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) [[download]](http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-centos_linux64.tar.gz)|
|htseq-count|counting the  reads mapped on to genomics feature|[[link]](http://htseq.readthedocs.io/en/master/install.html)|
|fastQC|widely used for sequencing read QC|[[link]](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) [[download]](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.6.zip)|
|RSeQC-2.6.4|An RNA-seq quality control package, multiple functions|[[link]](http://rseqc.sourceforge.net/)|
|trimmomatic|fastq quality trim and adaptor removal|[[link]](http://www.usadellab.org/cms/?page=trimmomatic) [[download]](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)|

This is a list of software already available on the cluster and the command you need to execute to load it in your environment: 

|Software|Description|Access on the cluster|
|------|-----------------------------|-----------------------|
|samtools|utilities for short DNA seq alignments [[link]](http://samtools.sourceforge.net/) |`module load samtools/1.3.1`|
|bedtools|tools for a wide-range of genomics analysis [[link]](http://bedtools.readthedocs.io/en/latest/)|`module load bedtools2/2.25.0`|
|bowtie2| alignment software [[link]](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)|`module load bowtie2/2.2.9`|
|tophat2| a fast splice junction mapper for RNA-Seq reads [[link]](https://ccb.jhu.edu/software/tophat/index.shtml)|`module load mvapich2/2.1  boost/1.59.0 tophat2/2.1.0`|
|R|language for statistical analysis|`module load intel/17.0.4 R-Project/3.4.1`|

This is a list of other software you might find useful:

|Software|Description|URL|
|------|-----------------------------|-----------------------|
| GSEA     | Genome set enrichment analysis| [[link]](http://software.broadinstitute.org/gsea/index.jsp) |
| IGV      | Interactive Genome Viewer      | [[link]](http://software.broadinstitute.org/software/igv/)|
| Cytoscape| Network visualization softwar   | [[link]](http://www.cytoscape.org/)|



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
   

## 2.  Download data

We will download human RNA-seq data with [GEO accession GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778). The samples  are in NCBI's short read archive format (SRA). 

We will use **sratoolkit** programs to download data but first we need to configure a location where all data files will be stored. <br>
`vdb-config` is a configuration subprogram for `sratoolkit`. We will use it to specify the directory where `sratoolkit` fetches data. You will need to type in the followingt path, but remember to replace `netID` with your own `Rutgers netid`  `/scratch/your_netID/Genomics_Workshop/download`.  Do not copy blindly!
Remember that  **sratoolkit**  is not designed to handle complex downloads. All data will be placed in one folder. You will need to move downloaded files for further analysis/manipulation into different locations manually.
```
                vdb-config   --interactive-mode textual     ### dash-dash before interactive-mode 
```
```
## Now in the program:
                
                         Your choice > 4
## type new path in
                        /scratch/your_netID/Genomics_Workshop/download
                        Your choice > Y
## Hit Enter and exit the program
```

Then execute the following commands to get the data. Both `prefetch` and `fastq-dump` are part of sratools.<br>
Downloading can be time consuming, it takes about 3 minutes per file. Thus we already downloaded files and placed them in your /scratch/$USER/Genomics_Workshop/download/sra folder. Skim through this section to understand how we did it and continue with 'fast-dump' below.  <br>

```
# fetch the SRA data
# Untreated
prefetch -v SRR1039508  
prefetch -v SRR1039512  
prefetch -v SRR1039516  
#DEX treated  
prefetch -v SRR1039509  
prefetch -v SRR1039513  
prefetch -v SRR1039517                          
```
Next step is to unpack sra files and convert them to more suitable *fastq* format with  `fastq-dump` program <br>
To unpack the original sequence files can be a bit tricky at first. We need to put them into different directories:<br>
*508,512,516* into /scratch/$USER/Genomics_Workshop/download/untreated and <br> 
 *509,513,517* into /scratch/$USER/Genomics_Workshop/download/dex_treated <br>
 Luckily `fastq-dump` can do processing and output results into specified folders at the same time. 

```
cd /scratch/$USER/Genomics_Workshop/download/sra
fastq-dump --outdir /scratch/$USER/Genomics_Workshop/untreated --split-files  SRR10395{08,12,16}.sra
fastq-dump --outdir /scratch/$USER/Genomics_Workshop/dex_treated --split-files  SRR10395{09,13,17}.sra
```
It takes a while to convert sra files. To save time, files are already converted for you. 
Run the following command to copy files into your /scratch/ directories.
```
sra-fastq.sh
```

# Running bioinformatics jobs

## 3. FastQC - raw data QC 

FastQC performs a quality control checks on raw sequence data and produces various graphical outputs for visual analysis.
```
        cd /scratch/$USER/Genomics_Workshop/untreated         
        module load java  ## fastqc is written in java; we need to load java before using fastqc
        mkdir fastqc      ## create a folder to store the QC output 
        fastqc -o fastqc SRR1039508_1.fastq SRR1039508_2.fastq
```
FastQC produces html pages in `fastqc/SRR1039508_1(2)_fastqc.html`, with different kinds of views of data (and Phred scores). You can open this file in Firefox browser. 
```
       firefox fastqc/SRR1039508_1_fastqc.html
```
To learn more about FastQC, see this pdf file <br>
/projects/oarc/Genomics_Workshop/RNA-Seq_analysis/misc/FastQC_details.pdf <br>
Close Firefox when you are done.

## 4. Trimmomatic - quality trim/adaptor removal

For demonstration purpose, we will take a small subset data using `seqtk` program
```
        cd /scratch/$USER/Genomics_Workshop/untreated
        seqtk sample -s100  SRR1039508_1.fastq 10000 > SRR1039508_1_10k.fastq 
        seqtk sample -s100  SRR1039508_2.fastq 10000 > SRR1039508_2_10k.fastq 
```   
More details and examples how to use `seqtk` can be found in <br> 
/projects/oarc/Genomics_Workshop/RNA-Seq_analysis/misc/Seqtk_Examples

Now, run `trimmomatic` to trim the read quality, and remove adaptor<br>

**NOTE:**  trimmomatic command starting with `java -jar` is a one line command, move the slider to the right to see the whole line. 
```
        module load java    ### needed for trimmomatic
        java -jar /home/$USER/Genomics_Workshop/Programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog trim.log SRR1039508_1_10k.fastq SRR1039508_2_10k.fastq SRR1039508_1.paired.fastq SRR1039508_1.unpaired.fastq SRR1039508_2.paired.fastq SRR1039508_2.unpaired.fastq ILLUMINACLIP:/home/$USER/Genomics_Workshop/Programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35
``` 
For your convenience we put this command into a bash script, thus you may just type
```
run_trimmo.sh
``` 

Once it started run, you should see the following:
```
        TrimmomaticPE: Started with arguments:
        -phred33 -trimlog trim.log SRR1039508_1_10k.fastq SRR1039508_2_10k.fastq SRR1039508_1.paired.fastq SRR1039508_1.unpaired.fastq SRR1039508_2.paired.fastq SRR1039508_2.unpaired.fastq ILLUMINACLIP:/home/yc759/Programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35
        Multiple cores found: Using 2 threads
        Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
        ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
        Input Read Pairs: 100000 Both Surviving: 96596 (96.60%) Forward Only Surviving: 1542 (1.54%) Reverse Only Surviving: 1467 (1.47%) Dropped: 395 (0.40%)
        TrimmomaticPE: Completed surccessfully 
```

View the output, the trim.log file, .e.g.  length=63 55 1 56 7 (the original read length 63, now 55 after trim, 1 base from left end and 7 bases from the right end were trimmed off, 56 bases in middle remained)
 
Alternatively, you may also use `fastx_quality_stats` from *the FASTX—toolkit* ( not covered in this workshop).


## 5. FastQC - Run on cleaned reads
You may run FastQC again on the cleaned by trimmomatic reads and compare new results with results for raw data, step 3 above.
```
     module load java
     fastqc -o fastqc SRR1039508_1.paired.fastq SRR1039508_2.paired.fastq
```

## 6. Download reference and reference indexing 

Indexing of Human genome may take hours. Thus we have the reference genome prepared and ready for use in <br>
`/projects/oarc/Genomics_Workshop/Reference/ `
<br>
For in class practice we demonstrate how to do genome indexing  with E.coli genome

```
        cd /scratch/$USER/Genomics_Workshop/Reference

        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/latest_assembly_versions/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz

        gunzip GCA_000005845.2_ASM584v2_genomic.fna.gz
        module load bowtie2
        bowtie2-build GCA_000005845.2_ASM584v2_genomic.fna GCA_000005845.2_ASM584v2_genomic
```
For your convenience we put these commands into a bash script, thus you may just type
```
run_bowtie2.sh
``` 


## 7. Mapping with tophat2, (STAR, HISAT2)
Before we start mapping job, make sure that you WON'T start a job on a login node.
tophat2 jobs are computationally intensive and may require hours to be completed. <br>

- **NOTE: We do not expect you to complete mapping in the class. You may start it to get familiar with the procedure, but do not wait till the completion. You may however do it as a homework for a preparation for the next session ( Genomics Workshop part II).** <br>

Also tophat2 is a multithreaded program, which means it can utilize more than one cpu-core, thus it is better to request more resources with srun command ( so far  we used `srun -c 2` (two cpu-cores).  For tophat2 jobs it is advisable to request more cpu-cores. For the convenience of this workshop we will use `-c 7`, however you may request more cpu-cores for your research.<br>
Exit the current interactive session: type ```exit``` in the terminal wondow. You should see that the prompt in the terminal changed from a compute node to a login node.
`netid@hal0011` to `netid@amarel`
Start a new interactive session requesting more cpu-cores (-c 7) <br>
```
srun  -p main --x11 --reservation=genomics -N 1 -c 7 -n 1 -t 01:40:00 --pty /bin/bash -i
```
For your convenience we put this command into a bash script, thus you may just type
```
tophat2_node_request.sh
``` 
Notice that the prompt changed again to a compute node, e.g. `netid@hal0011` 

 Go to your folder with data.<br>
 Remember that you should match the `tophat2 -p` option  to be consistent with the number of cores `srun -c` that you requested, thus we use `tophat2 -p 7`, because we requested `srun -c 7`. <br>
**NOTE:**  tophat2 command starting with `tophat2 -p 7` is a one line command, move the slider to the right to see the whole line.        
```
        cd  /scratch/$USER/Genomics_Workshop/untreated
        module load mvapich2/2.1  boost/1.59.0  tophat2/2.1.0
        module load samtools  bowtie2 
        mkdir tophat_out
        tophat2 -p 7 --library-type fr-unstranded  -o tophat_out/untreated_SRR1039508_10k --transcriptome-index /projects/oarc/Genomics_Workshop/Reference/hg20_transciptome/GRCh38.78 /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.dna.toplevel SRR1039508_1.paired.fastq SRR1039508_2.paired.fastq
```
For your convenience we put these commands into a bash script, thus you may just type
```
run_tophat2.sh
``` 

You should see a similar output:
```
   [2018-03-30 11:48:57] Beginning TopHat run (v2.1.0)
-----------------------------------------------
[2018-03-30 11:48:57] Checking for Bowtie
                  Bowtie version:        2.2.9.0
[2018-03-30 11:48:58] Checking for Bowtie index files (transcriptome)..
[2018-03-30 11:48:58] Checking for Bowtie index files (genome)..
[2018-03-30 11:48:58] Checking for reference FASTA file
[2018-03-30 11:48:58] Generating SAM header for /projects/oarc/Genomics_Workshop/Referen
ce/hg20/Homo_sapiens.GRCh38.dna.toplevel
[2018-03-30 11:49:09] Reading known junctions from GTF file
[2018-03-30 11:49:27] Preparing reads
         left reads: min. length=35, max. length=63, 96592 kept reads (4 discarded)
        right reads: min. length=35, max. length=63, 96594 kept reads (2 discarded)
[2018-03-30 11:49:29] Using pre-built transcriptome data..
[2018-03-30 11:49:35] Mapping left_kept_reads to transcriptome GRCh38.78 with Bowtie2
[2018-03-30 11:49:49] Mapping right_kept_reads to transcriptome GRCh38.78 with Bowtie2
[2018-03-30 11:50:03] Resuming TopHat pipeline with unmapped reads
[2018-03-30 11:50:03] Mapping left_kept_reads.m2g_um to genome Homo_sapiens.GRCh38.dna.t
oplevel with Bowtie2
[2018-03-30 11:50:16] Mapping left_kept_reads.m2g_um_seg1 to genome Homo_sapiens.GRCh38.dna.toplevel with Bowtie2 (1/2)
[2018-03-30 11:50:18] Mapping left_kept_reads.m2g_um_seg2 to genome Homo_sapiens.GRCh38.dna.toplevel with Bowtie2 (2/2)
[2018-03-30 11:50:20] Mapping right_kept_reads.m2g_um to genome Homo_sapiens.GRCh38.dna.toplevel with Bowtie2
[2018-03-30 11:50:23] Mapping right_kept_reads.m2g_um_seg1 to genome Homo_sapiens.GRCh38.dna.toplevel with Bowtie2 (1/2)
[2018-03-30 11:50:25] Mapping right_kept_reads.m2g_um_seg2 to genome Homo_sapiens.GRCh38.dna.toplevel with Bowtie2 (2/2)
……………………………………………………….
```
Here is another example of mapping with tophat. <br>
The transcriptome index was built  by pointing to a gtf file first (the gtf file is prepared earlier).  The following would be the command to generate the transcriptome index while running tophat alignment.<br>
**NOTE:**  tophat2 command starting with `tophat2 -p 7` is a one line command, move the slider to the right to see the whole line. 
```
   tophat2 -p 7 --library-type fr-unstranded  -o tophat_out/untreated_SRR1039508 –GTF /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.78.gtf --transcriptome-index /projects/oarc/Genomics_Workshop/Reference/hg20_transciptome/GRCh38.78 /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.dna.toplevel SRR1039508_1.fastq.pairedOut.fastq SRR1039508_2.fastq.pairedOut.fastq
```
After the completion the output folder `tophat_out/untreated_SRR1039508/` shall contain the following files/folders:  
```
 -rw-rw-r-- 1 yc759 oarc 2174796848 Jan 16 21:57 accepted_hits.bam
 -rw-rw-r-- 1 yc759 oarc        565 Jan 16 21:57 align_summary.txt
 -rw-rw-r-- 1 yc759 oarc    1921529 Jan 16 21:57 deletions.bed
 -rw-rw-r-- 1 yc759 oarc    2239884 Jan 16 21:57 insertions.bed
 -rw-rw-r-- 1 yc759 oarc   14181618 Jan 16 21:57 junctions.bed
 drwxrwsr-x 2 yc759 oarc       4096 Jan 16 21:57 logs
 -rw-rw-r-- 1 yc759 oarc        184 Jan 16 21:57 prep_reads.info
 -rw-rw-r-- 1 yc759 oarc   42846571 Jan 16 21:57 unmapped.bam
```

## 8. Read counts using htseq-count

Tophat mapping is resource intense and time consuming. If you managed to produce the bam file yourself, it's great. If not, let's move forward and copy files that we prepared for you.

``` 
     cd /scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508 
     cp /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/SRA_data/untreated/tophat_out/untreated_SRR1039508/accepted_hits.bam accepted_hits.bam   
```

Next, we will use samtools to sort the bam file by name:  because htseq-count accepts bam file sorted by **name** as default, but tophat generates bam sorted by **coordinates** as default. 
```
     module load samtools intel/17.0.2 python/2.7.12   #loads the software
     samtools sort -n  accepted_hits.bam | samtools view | htseq-count -m intersection-nonempty -t exon -i gene_id -s no --additional-attr=gene_name - /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.78.gtf > untreated08.txt
```
In the same way, generate the counts file `untreated12.txt`, `untreated16.txt`,`dex09.txt`, `dex13.txt`, `dex17.txt`. <br>
To save time and for your convenience we prepared files already, thus you may just type
```
run_htseq_count.sh
```

## 9. Perform Mapping QC using RSeQC

Now,  quality control using RSeQC - a few examples are given here, please go to the [rseqc website](http://rseqc.sourceforge.net/) for more functions.<br> 
- Example 1 <br>
Execute:   
```
cd /scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508 
module load python/2.7.12
module load intel/17.0.4
     
read_distribution.py -i accepted_hits.bam -r /projects/oarc/Genomics_Workshop/Reference/Homo_sapiens.GRCh38.79.bed
```
Provided a BAM/SAM file and reference gene model, `read_distribution.py` will calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions).<br>
Output: 
```
        processing/projects/oarc/Genomics_Workshop/Reference/Homo_sapiens.GRCh38.79.bed ... Done
        processing accepted_hits.bam ... Finished

        Total Reads                   43474036
        Total Tags                    54438789
        Total Assigned Tags           53991382
        =====================================================================
        Group               Total_bases         Tag_count           Tags/Kb
        CDS_Exons           103371993           43264842            418.54
        5'UTR_Exons         5217678             583447              111.82
        3'UTR_Exons         29324747            8145122             277.76
        Introns             1500197093          1805034             1.20
        TSS_up_1kb          33306654            18893               0.57
        TSS_up_5kb          148463534           41165               0.28
        TSS_up_10kb         265823549           55644               0.21
        TES_down_1kb        35215293            50954               1.45
        TES_down_5kb        152556214           113325              0.74
        TES_down_10kb       268614580           137293              0.51
        =====================================================================
```
- Example 2 <br>
Execute:  
```
bam_stat.py -i accepted_hits.bam
```
