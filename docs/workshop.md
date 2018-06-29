# Overview

These are notes from the workshop on genomics. What's covered in these notes: 
1. Setup


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

This is a list of other software you might find useful:

|software|description|link|
|------|-----------------------------|-----------------------|
| GSEA     | genome set enrichment analysis| [link](http://software.broadinstitute.org/gsea/index.jsp) |
| IGV      | Interactive Genome Viewer      | [link](http://software.broadinstitute.org/software/igv/)|
| Cytoscape| Network visualization softwar   | [link](http://www.cytoscape.org/)|

# Setup

### Connect to the cluster login node 
Do this by one of the following methods: 

  - **via a terminal**: if you have a Mac or Linux, terminal is part of your standard apps. If you have Windows, install an SSH client such as `putty` or `moba-xterm`. Then from your terminal connect to the cluster by executing the following command:   
``` ssh -X <your net id>@amarel.hpc.rutgers.edu```   
   This is the preferred way, as your copy-pasting will most likely work best. 
  - **via FastX**: in your browser, go to `https://amarel.hpc.rutgers.edu:3443`

### Get resources on the compute node 

You get to the cluster to execute your computations by running the following command in your terminal:    
```srun  -p main --reservation=genomics -N 1 -c 2 -n 1 -t 01:40:00 --export=ALL --pty /bin/bash```   
Notice that the name in your terminal will change from `amarel` to node name like `hal0025` or `slepner086`. This means that you will not impede other users who are also using the login node, and will be placed on a machine which you share with only a few people. The following table explains the parts of this command: 

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

### Prepare some directories for the data

You have two main spaces on the Amarel cluster. These are: 

- your home directory - `/home/netid/`  (e.g. `/home/kp807/` for my netid) -
- your scratch directory - `/scratch/netid/`

  They differ in how often they are backed up and by size (100Gb for `/home` and 500Gb for `/scratch`). So we will install programs in `/home`, while the data and output will be in `/scratch`. Execute these commands: 
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

## Install programs

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
Moreover, execute the following commands to load system-installed sofware so the system knows where to find it (i.e. `samtools` command will work if you execute `module load samtools`): 
```
module load samtools       
module load bedtools2./2.25.0
module load bowtie2
module load mvapich2/2.1  boost/1.59.0  tophat2/2.1.0
module load intel/17.0.4 R-Project/3.4.1
```

# Download data

We will download human RNA-seq data with [GEO accession GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778). The samples we download are in NCBI's short read archive format (SRA). To unpack the original sequence files can be a bit tricky at first. Please put them in different directories:
```        
                mkdir -p  /scratch/$USER/Genomics_Workshop/
                cd /scratch/$USER/Genomics_Workshop/
                mkdir untreated  
                mkdir dex_treated
```

We will use sratoolkit programs to download data but first we need to configure a location where all data files will be stored. `sratoolkit` will be in your home directory, under `Programs`, and the `vdb-config` might be under the `bin` directory. You will enter `/scratch/your_netID/Genomics_Workshop/download` for the path - NOTE you have to replace `your_netID` with your true netId, e.g. `kp807`. Do not copy blindly! So your downloads will always go to this directory and you will need to move it out to wherever you want to have them. 
```
                vdb-config   --interactive-mode textual     ### dash-dash before interactive-mode
                         Your choice > 4
## type new path in
                        /scratch/your_netID/Genomics_Workshop/download
                        Your choice > Y
```

Then execute the following commands to get the data. Both `prefetch` and `fastq-dump` are part of sratools. Downloading can take some time! [TODO: check how much time for these files!]
```
prefetch -v SRR1039508                           # fetches the SRA data
fastq-dump --gzip --split-files SRR1039508       # ???? 
```
You have to pay attention to where you are putting your data. So these two commands will actually be several: 
```
                cd  untreated                       # now you are in /scratch/..../Genomics_Workshop/untreated
                prefetch -v SRR1039508
                mv /scratch/$USER/Genomics_Workshop/download/sra/SRR1039508.sra .  # moving from download to actual directory 
                fastq-dump --gzip --split-files SRR1039516
``` 

The commands above showed how to do it for one sample. You need to do it for 6 samples total. 
```
                SRR1039508  SRR1039512 SRR1039516   (untreated)
                SRR1039509  SRR1039513  SRR1039517  (dex_treated)
```

# Running bioinformatics jobs

### FastQC - raw data QC and quality trim/adaptor removal

Explain what is fastqc is doing here - TODO
```
        cd /scratch/$USER/Genomics_Workshop/untreated         
        module load java  ## fastqc is written in java; we need to load java before using fastqc
        mkdir fastqc      ## create a folder to store the QC output 
        fastqc -o fastqc SRR1039508_1.fastq SRR1039508_2.fastq
```
FastQC produces an html page as output, `fastqc/SRR1039508_1_fastqc.html`, with different kinds of views of data (and Phred scores). You can download this file to your local machine and open it in browser. It is also possible to open browser on the cluster, but the cluster is not really designed for that. To see more about FastQC, see this pdf file - /projects/oarc/Genomics_Workshop/Labs/FastQC_details.pdf

### Trimmomatic - quality trim/adaptor removal

        ##for demonstration purpose, we will take a small subset data using seqtk
        cd /scratch/$USER/Genomics_Workshop/untreated
        seqtk sample -s100  SRR1039508_1.fastq 10000 > SRR1039508_1_10k.fastq 
        seqtk sample -s100  SRR1039508_2.fastq 10000 > SRR1039508_2_10k.fastq 
        ## /projects/oarc/Genomics_Workshop/Labs/Seqtk_Examples.docx
        ## This file contains useful examples how to use seqtk 

        ##now, run trimmomatic to trim the read quality , and remove adaptor
        module load java    ### because trimmomatic
        java -jar /home/$USER/Genomics_Workshop/Programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog trim.log SRR1039508_1_10k.fastq SRR1039508_2_10k.fastq SRR1039508_1.paired.fastq SRR1039508_1.unpaired.fastq SRR1039508_2.paired.fastq SRR1039508_2.unpaired.fastq ILLUMINACLIP:/home/$USER/Programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35

**NOTE:**  the above is a one line command, illustrated as the following:
```
        java -jar trimmomatic-0.36.jar PE \
        -phred33 -trimlog trim.log \
        input_1.fq  input_2.fq \
        output_1_paired.fq  output_1_unpaired.fq \
        output_2_paired.fq  output_2_unpaired.fq \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35 / 

        ## Once it started run, you shall see the following:
        TrimmomaticPE: Started with arguments:
        -phred33 -trimlog trim.log SRR1039508_1_10k.fastq SRR1039508_2_10k.fastq SRR1039508_1.paired.fastq SRR1039508_1.unpaired.fastq SRR1039508_2.paired.fastq SRR1039508_2.unpaired.fastq ILLUMINACLIP:/home/yc759/Programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35
        Multiple cores found: Using 2 threads
        Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
        ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
        Input Read Pairs: 100000 Both Surviving: 96596 (96.60%) Forward Only Surviving: 1542 (1.54%) Reverse Only Surviving: 1467 (1.47%) Dropped: 395 (0.40%)
        TrimmomaticPE: Completed successfully

        ##view the output, the trim.log file, .e.g.  length=63 55 1 56 7 (the original read length 63, now 55 after trim, 1 base from left end and 7 bases from the right end were trimmed off, 56 bases in middle remained)
 
   ##you may also try fastx_quality_stats from the FASTX—toolkit

```

### FastQC - Run on cleaned reads, compare result
```
     module load java
     fastqc -o fastqc SRR1039508_1.paired.fastq SRR1039508_2.paired.fastq

     ## /projects/oarc/Genomics_Workshop/Labs/FastQC_details.pdf , helpful in viewing and interpreting the output
```

### Download reference and reference indexing 

Human genome indexing will take hours, we have the reference pre-prepared. Stored at  `/projects/oarc/Genomics_Workshop/Reference/ `
For in class practice, we will do this on E.coli genome

```
       cd /scratch/$USER/Genomics_Workshop/
        mkdir Reference
        cd Reference

        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/latest_assembly_versions/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz

        gunzip GCA_000005845.2_ASM584v2_genomic.fna.gz
        module load bowtie2
        bowtie2-build GCA_000005845.2_ASM584v2_genomic.fna GCA_000005845.2_ASM584v2_genomic

###if download from ENSEMBLE
        wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-38/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz
```

### Mapping with tophat2, (STAR, HISAT2)

Now, go to your data folder
        `cd  /scratch/$USER/Genomics_Workshop/untreated`

```
cd  /scratch/$USER/Genomics_Workshop/untreated
        module load mvapich2/2.1  boost/1.59.0  tophat2/2.1.0
        module load samtools   #bowtie2 is loaded already
        mkdir tophat_out
        tophat2 -p 10 --library-type fr-unstranded  -o tophat_out/untreated_SRR1039508_1
        0k --transcriptome-index /projects/oarc/Genomics_Workshop/Reference/  hg20_transciptome/GR
        Ch38.78 /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.dna.toplevel
        SRR1039508_1.paired.fastq SRR1039508_2.paired.fastq
## you shall modify the -p value to be consistent with the -c value you requested in the beginning
```
You shall see something like:
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
The transcriptome index was built  by pointing to gtf file first,  here we have it prepared already, just so we can save time.  The following would be the command to generate the transcriptome index while running tophat alignment.
```
   tophat2 -p 10 --library-type fr-unstranded  -o tophat_out/untreated_SRR1039516 –GTF /projects/oarc/Genomics_Workshop/Reference/hg20/ Homo_sapiens.GRCh38.78.gtf --transcriptome-index /projects/oarc/Genomics_Workshop/Reference/hg20_transciptome/GRCh38.78 /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.dna.toplevel SRR1039516_1.fastq.pairedOut.fastq SRR1039516_2.fastq.pairedOut.fastq

```
The output folder `tophat_out/untreated_SRR1039508/` shall contain the following files/folder (in blue): 
```
     cd  /projects/oarc/Genomics_Workshop/SRA_data/untreated/tophat_out/untreated_SRR1039508
        $ ll
   total 2183632
 -rw-rw-r-- 1 yc759 oarc 2174796848 Jan 16 21:57 accepted_hits.bam
 -rw-rw-r-- 1 yc759 oarc        565 Jan 16 21:57 align_summary.txt
 -rw-rw-r-- 1 yc759 oarc    1921529 Jan 16 21:57 deletions.bed
 -rw-rw-r-- 1 yc759 oarc    2239884 Jan 16 21:57 insertions.bed
 -rw-rw-r-- 1 yc759 oarc   14181618 Jan 16 21:57 junctions.bed
 drwxrwsr-x 2 yc759 oarc       4096 Jan 16 21:57 logs
 -rw-rw-r-- 1 yc759 oarc        184 Jan 16 21:57 prep_reads.info
 -rw-rw-r-- 1 yc759 oarc   42846571 Jan 16 21:57 unmapped.bam
```

### Read counts using htseq-count

GO TO WHERE YOUR ALIGNMENT OUTPUT FOLDER IS, FOR EXAMPLE: 

``` 
  cd /scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508 
  ln –s /projects/oarc/Genomics_Workshop/SRA_data/untreated/tophat_out/untreated_SRR1039508/accepted_hits.bam accepted_hits.bam   

  ##make a soft link to the full bam file we already prepared, if you didn’t have the bam ready yet

     module load samtools intel/17.0.2 python/2.7.12 
     samtools sort -n  accepted_hits.bam | samtools view | htseq-count -m intersection-nonempty -t exon -i gene_id -s no --additional-attr=gene_name  -/projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.78.gtf > untreated08.txt
```
Use samtools to sort the bam file by name:  because htseq-count accepts bam file sorted by name as default, but tophat generates bam sorted by coordinates by default
The same way to generate the counts file *untreated12.txt*, *untreated16.txt*,*dex09.txt*, *dex13.txt*, *dex17.txt*

### Perform Mapping QC using RSeQC

Now,  quality control using RSeQC –a few examples here, please go to the website for more functions `http://rseqc.sourceforge.net/`
```
  cd /scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508 
        module load python/2.7.12
        module load intel/17.0.4
     
        $ read_distribution.py -i accepted_hits.bam -r /projects/oarc/Genomics_Workshop/
Reference/Homo_sapiens.GRCh38.79.bed
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

 $ bam_stat.py -i accepted_hits.bam
        Load BAM file ...

        Done

        #==================================================
        #All numbers are READ count
        #==================================================

        Total records:                          52528699

        QC failed:                              0
        Optical/PCR duplicate:                  0
        Non primary hits                        9054663
        Unmapped reads:                         0
        mapq < mapq_cut (non-unique):           2684801

        mapq >= mapq_cut (unique):              40789235
        Read-1:                                 20414530
        Read-2:                                 20374705
        Reads map to '+':                       20393901
        Reads map to '-':                       20395334
        Non-splice reads:                       30860931
        Splice reads:                           9928304
        Reads mapped in proper pairs:           32386536
        Proper-paired reads map to different chrom:312

```
The script does genebody coverage calculation requires the input bam files to be sorted and indexed (we will do it using samtools). The calculation and plot will require R
Go to one of the tophat_out sample folder
```
$ module load intel/17.0.4  R-Project/3.4.1    
        $ module load samtools
        $ samtools sort accepted_hits.bam –o accepted_hits.sorted.bam 
        ##this may take a while, you may use the one already prepared for you by making a soft link
        ln –s /projects/oarc/Genomics_Workshop/SRA_data/untreated/tophat_out/untreated_SRR1039508/accepted_hits.sorted.bam accepted_hits.sorted.bam


        $ samtools index accepted_hits.sorted.bam
 
        $ geneBody_coverage.py -r /projects/oarc/Genomics_Workshop/Reference/hg38.housekeepingGenes.bed -i accepted_hits.sorted.bam -o test
@ 2018-01-14 13:17:33: Read BED file (reference gene model) ...
@ 2018-01-14 13:17:33: Total 3802 transcripts loaded
@ 2018-01-14 13:17:33: Get BAM file(s) ...
        accepted_hits.sorted.bam
@ 2018-01-14 13:17:33: Processing accepted_hits.sorted.bam ...
        3800 transcripts finished


        Sample  Skewness
        accepted_hits.sorted    -3.61577607436
@ 2018-01-14 13:28:59: Running R script ...
null device
          1

```
output files:        test.geneBodyCoverage.r
                                test.geneBodyCoverage.txt
                log.txt
                test.geneBodyCoverage.curves.pdf

download the gene.bed files :  `https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/ `  
Be careful that the genome version, be consistent between reference genome used in mapping and now. For now, you may use what’s provided 

if you want to find out whether the sequencing read is strand specific or not, do:
```
        $ infer_experiment.py -r  /projects/oarc/Genomics_Workshop/Reference/ Homo_sapiens.GRCh38.79.bed -i accepted_hits.bam

    Reading reference gene model /projects/oarc/Genomics_Workshop/Reference/Homo_sapiens.GRCh38.79.bed ... Done
    Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.1406
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4302
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4292
```





