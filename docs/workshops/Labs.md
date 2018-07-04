# LAB 1: Visualization using IGV
___
 **IGV**  --- Focus on visualization, best for validation and confirmation of the analysis result, Not good for primary analysis
 The mapping file is in bam format, located under the folder of tophat_out, they shall be sorted and indexed using the following command

```
  cd /scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508
  module load samtools                                                             
  samtools sort accepted_hits.bam -o accepted_hits.sorted.bam  ##this step  takes about 10 minutes to complete 
  samtools index accepted_hits.sorted.bam ## it takes about 30 seconds
```
The resulting files: accepted_hits.sorted.bam    
                     accepted_hits.sorted.bam.bai  
                     are the files to be uploaded to IGV

You need to repeat these steps for every sample 

## Shortcut Lab 1

We have prepared 4 sets of such files (dex09, dex13, untreated08 and untreated12), located at: `/projects/oarc/Genomics_Workshop/Bam_for_IGV/` . Make a soft link (see the following command), or copy them into your scratch  folder, then we use IGV to analyze them. 

```
cd /scratch/$USER/Genomics_Workshop/
ln -s /projects/oarc/Genomics_Workshop/Bam_for_IGV  Bam_for_IGV  ## this step was  done when you ran lab_PartII.sh
#start IGV 
module load java
/projects/oarc/Genomics_Workshop/IGV_2.4.6/igv.sh 
```  
 Practice and get familiar with:  

   * How to Load genome and data track
   * How to navigate
   * How and what to visualize:
   * Examine coverage
   * Low mapping quality
   * Mis-alignment
   * Translocation
   * Novel genes/transcript
   * Alternative splicing
   * Inversion
   * Look for SNPs
   * CNV, ChipSeq, RNASeq, WGS alignmentSNP

 More detailed explanation [here](http://software.broadinstitute.org/software/igv/book/export/html/37)

 The following is a [sample](HTTP://SOFTWARE.BROADINSTITUTE.ORG/SOFTWARE/IGV/) snap shot of the above two files loaded to IGV.  
  CLOSE your interactive session on a node  when done with IGV by typing exit


# LAB 2:  Data processing and expression analysis using edgeR
___
 All the htseq-count output files should be present in one folder. Here we created the folder read_counts.
```
            read_counts/dex09.txt
                        dex13.txt
                        dex17.txt
                        untreated08.txt
                        untreated12.txt
                        untreated16.txt
```
 Go to this folder containing all above six counts output files  
 `cd /scratch/$USER/Genomics_Workshop/read_counts`

 We need to combine all counts into one file, which will be imported into R for further analysis.
 The following shows how this can be done with bash commands on linux ( you may do it in excel too)
 
```
 paste dex09.txt dex13.txt dex17.txt untreated08.txt untreated12.txt untreated16.txt > merged_counts.txt    
 ## merge files by columns
 cut -f1,2,4,6,8,10,13 merged_counts.txt  > merged_counts_clean.txt     ## extract the relevant columns  
 head -n -5 merged_counts_clean.txt > merged_counts_clean2.txt      ##remove the last 5 line stats summary
```

We also need to prepare a file containing group/treatment information. This file Targets.txt is a tab delimited text file. You may construct in excel. 
 **The file should contain the following info**

 label  | sample |group | treatment
------  |------- |------|--------------
1       |dex09 |  dex   | dex_treated
2       |dex13 | dex    |dex_treated
3       |dex17 | dex    |dex_treated
4       |untreated08    ||control untreated
5       |untreated12    ||control untreated
6       |untreated16    ||control untreated

This file has been prepared for you. You will need to copy it into your folder, see later Lab2.2:
                        
## Lab2.2: Pre-processing the data in R

```
cd /scratch/$USER/Genomics_Workshop/
mkdir  DE_analysis  ##set working directory to run differential expression analysis
cd DE_analysis 

## Copy the needed files here

cp /projects/oarc/Genomics_Workshop/SRA_data/DE_analysis/Targets.txt $PWD                           ##This file denotes the experimental group
cp /projects/oarc/Genomics_Workshop/SRA_data/read_counts/merged_counts_clean2.txt $PWD              ##This file contains read counts on genes for all samples
cp /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.78.gtf $PWD     ##annotation file needed to calculate exonic gene length -- needed for FPKM calculation
```
##  HOMEWORK CATCH-UP
 From your homework assignment you should have the following packages to be installed already :).
If you didn't, install them now.

```
##Open a terminal for  amarel2 login node.
ssh -X amarel2.hpc.rutgers.edu

##On the login node start R
module load intel/17.0.4
module load R-Project/3.4.1
R

##then in the R workspace do the following:
source("https://bioconductor.org/biocLite.R")
biocLite("MKmisc")
Would you like to use a personal library instead?  (y/n)  y
Would you like to create a personal library
~/R/x86_64-pc-linux-gnu-library/3.4 to install packages into?  (y/n) y

##Wait till it finishes.

biocLite("Heatplus")
biocLite("affycoretools")
biocLite("flashClust")
biocLite("affy")
biocLite("GenomicFeatures")
quit()   ###quit R, no save
Save workspace image? [y/n/c]: n
```

## Starting the Job

Now,start a new interactive job on the compute node or switch to another terminal if you still have an interactive job running  
`srun --x11 -p main --reservation=genomics_2 -N 1 -c 2 -n 1 -t 02:00:00 --pty /bin/bash -i`     
Go to your working directory  
 `cd /scratch/$USER/Genomics_Workshop/DE_analysis/`
 Then, start R  
 ```
 module load intel/17.0.4
 module load R-Project/3.4.1 
 ```  
 `getwd()` Check which directory you are in  
 You should be in the directory where the files merged_counts_clean2.txt and "Targets.txt", annotation gtf file are. If not there, set your working directory :
 `setwd("/scratch/<your_netID>/Genomics_Workshop/DE_analysis/")`. Set the working directory

 Now load up the libraries needed for the analysis

        library(edgeR)
        library(MKmisc)
        library(affy)
        library(flashClust)
        library(affycoretools)
        library(Heatplus)
        library(GenomicFeatures)

```
raw.data <- read.delim("merged_counts_clean2.txt", header=F, row.names=1)  
##import count data to R
head(raw.data)  ## check the beginning portion of the imported file, now an object
dim(raw.data)   ##check the dimention of this object
class(raw.data)  ## check what class of this object
apply(raw.data, 2, summary) ## check the range of counts per sample
range(rowSums(raw.data))      ## check the range of counts per gene
colnames(raw.data) <- c("dex09","dex13","dex17","untreated08","untreated12","untreated16")
##add column header

raw.data2 <- raw.data[rowSums(raw.data) != 0 ,] ##remove genes with zero counts for all samples
dim(raw.data2)

cpm.values <- cpm(raw.data2) #calculate counts per million mapped reads without any other normalization

above1cpm <- apply(cpm.values, 1, function(x) sum(x >=1)) ##How many samples/genes had at least 1 cpm 
counts.use <- raw.data2[above1cpm >= 3,] ##we have three replicates in each group. If a gene can be reliably detected, it should be detected in all 3 replicates. So, a gene to be included for further analysis shall have 1 cpm in at least 3 samples. (The 3 samples are irrespective of group)
dim(counts.use)
colSums(counts.use) / colSums(raw.data)  ##the % of total counts kept after filtering
nrow(counts.use) / nrow(raw.data)  ##the % of genes kept after filtering

targets <- readTargets()      ##import targets file that contains experiment group info.
targets$GpF <- factor(targets$group)   ##change character to factor
targets$GpN <- as.numeric(targets$GpF) ##change factor to numeric (optional)

ls()  #see that objects have been loaded
save.image("RNASeqDemo.RData")   ##save the work workspace
savehistory("RNASeqDemo.Rhistory")  ##save command history

```
## Lab 2.3: To calculate expression values as fpkm

First, compute the gene length as described in Lab6--partI
```
 library(GenomicFeatures)  ##may skip, because we already loaded at start
 gtfdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.78.gtf",format="gtf")
 exons.list.per.gene <- exonsBy(gtfdb,by="gene")
 exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
 class(exonic.gene.sizes)
 Hg20_geneLength <-do.call(rbind, exonic.gene.sizes)
 colnames(Hg20_geneLength) <- paste('geneLength')    

 Hg20_geneLength2 <- data.frame(Hg20_geneLength[rownames(counts.use),]) ## to extract the gene length file to contain 
 the same number genes in the same order as in the filtered read counts file
 colnames(Hg20_geneLength2) <- paste('geneLength')  ## to change column name, make    it neat
 fpkm.data <- cpm(counts.use) / (Hg20_geneLength2$geneLength/1000) ## compute fpkm
 min.count <- apply(fpkm.data,1,min)
 write.csv(fpkm.data,file="fpkm_values.csv")  #### To output FPKM data  

```

## Lab 2.4:  Analysis QC ---or sample diagnosis

```
## density distribution
plotDensity(log2(raw.data+0.1),col=targets$GpF,lty=1,lwd=2)
legend("topright", legend=levels(targets$GpF),fill=1:4)
```

Change data from raw.data to raw.data2, to CPM, FPKM,.. to see the effect of filtering and normalization

### Clustering

```
hc.raw <- flashClust(dist(t(log2(raw.data2+0.1))),method="average")
plot(hc.raw,hang = -1, main = "RNASeqDemo, raw values", sub = "", xlab = "",cex=0.9, labels=targets$sample)
##change data from raw.data to raw.data2, to CPM, logCPM, FPKM,.. to see the effect of filtering and normalization
 
#####PCA#######
plotPCA(fpkm.data, pch=16, col=targets$GpF,groupnames=levels(targets$GpF), addtext=rep(1:3,4),main="PCA on FPKM")

###do it after edgeR analysis, otherwise some values not existing yet###

plotPCA(logCPM), pch=16, col=targets$GpF,groupnames=levels(targets$GpF), addtext=rep(1:3,4),main="PCA on logCPM")
 

#####heatmap####
test <- topTags(eR.dex_Ctl,n=Inf,sort.by="PValue")$table ####sort gene list according to P values
test2 <- test[1:500,]  ###Take the top 500 significant genes
logCPM2 <- logCPM[rownames(test2),]  ###Extract logCPM values of these 500 genes
meanFC <- logCPM2 - rowMeans(logCPM2)
color.meanFC <- heatmapCol(data = meanFC, lim = 3, col =colorRampPalette(c("blue","white","red"))(128))
heatmap_2(meanFC, col=color.meanFC, legend=3, scale="none")

 

 #####Volcano plot####
with(eR.dex_Ctl.detailed, plot(logFC, -log10(PValue), pch=20, main="Volcano plot", xlim=c(-2.5,2),ylim=c(0,25)))
with(subset(eR.dex_Ctl.detailed,FDR<0.05), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(eR.dex_Ctl.detailed, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(eR.dex_Ctl.detailed, FDR<0.05&abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))

```

## Lab 2.5: Differential expression analysis  with edgeR

Create the design matrix
```
  groups <- factor(targets$group)
  design <- model.matrix(~0+groups)
  colnames(design) <- levels(groups)
  rownames(design) <- targets$sample   
  design
```

```
                                        control dex
dex09             0   1
dex13             0   1
dex17             0   1
untreated08       1   0
untreated12       1   0
untreated16       1   0
attr(,"assign")
[1] 1 1
attr(,"contrasts")
attr(,"contrasts")$groups
[1] "contr.treatment"
```

### Create contrast matrix
` cont.matrix <- makeContrasts(dex_Ctl= dex - control, levels=design)`

### Create DGEList object
```
 d <- DGEList(counts=counts.use, lib.size=colSums(counts.use), group=targets$GpF)
 class(d)
 names(d)    ## the names of the items in the list
 d$counts[1:5,]    ## The counts are stored in the $counts:
 d$samples   ## The group info and library sizes stored in $samples
 d <- calcNormFactors(d)  ## an additional normalization factor using a TMM method
 d$samples   ## now the norm.factors are no longer 1

 ###this can be useful when diagnose problem)###
 apply(d$counts,2,function(x) sum(sort(x/sum(x),decreasing=T)[1:20])) * 100
 ##the proportions of total counts for the top 20 genes in each sample, 10-20% is OK. 

 logCPM <-  cpm(d$counts, log = TRUE)   ## modified logCPM values in edgeR, can be used for clustering, heatmap

```
Now,the DE test! The term "estimating dispersions" in edgeR describes a method to account for the variance among 
replicates.

```
  d <- estimateGLMCommonDisp(d, design, verbose=TRUE)   
  Disp = 0.06037 , BCV = 0.2457 
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design)
  plotBCV(d)   ## the relationship between the overall abundance and the tagwise dispersion estimates


  fit.edgeR <- glmFit(d, design)  ## Estimate model coefficients from count data and design matrix
  names(fit.edgeR)
```
Specify contrasts of interest, do empirical Bayes "shrinkage" of  variances and calculate test statistics. Both of these are performed with same function glmLRT in edgeR (Genewise Negative Binomial Generalized Linear Models)

```
 eR.dex_Ctl <- glmLRT(fit.edgeR, contrast=cont.matrix[,1])
  
 summary(decideTestsDGE(eR.dex_Ctl)) [,1]  ## Correct for multiple tests and extract relvant data (1 means sig up, -1 means sig down, and 0 means NS)

 eR.dex_Ctl.detailed <- topTags(eR.dex_Ctl,n=Inf,sort.by="none")$table ## Get detailed output for a single contrast
 
 eR.dex_Ctl.detailed[1:5,]
                      logFC   logCPM         LR    PValue       FDR
 ENSG00000000003 -0.36354885 5.227841 2.37143785 0.1235732 0.5501740
 ENSG00000000419  0.18626588 4.674100 0.79016266 0.3740509 0.8571700
 ENSG00000000457  0.04304803 3.838702 0.03317218 0.8554789 0.9905954
 ENSG00000000460 -0.05998067 1.609803 0.02021703 0.8869326 0.9916274
 ENSG00000000971  0.35401120 8.007586 1.55551391 0.2123233 0.7070242

 ## Fold change is simply group A-B (if on the log scale), or A/B (if on   raw scale).
 ## logCPM = the average log2-counts-per-million 
 ## LR = likelihood ratio statistics 
 ##PValue = the two-sided p-value 
 ## FDR = adjusted p-value 

 #### To back-translate logFC to regular FC with down-reg as -FC
 eR.dex_Ctl.detailed$FC <- 2^abs(eR.dex_Ctl.detailed$logFC) * sign(eR.dex_Ctl.detailed$logFC)

 write.csv(eR.dex_Ctl.detailed,file="Demo_eR.dex_Ctl_results.csv")  ##output the file

```

# Lab 3  Running enrichment analysis using GSEA 
___
The command to start the gsea:
```
srun --x11 -p main --reservation=genomics_2 -N 1 -c 2 -n 1 -t 2:00:00 --pty /bin/bash -i  ##get onto a reserved compute node
module load java
java -jar /scratch/$USER/Genomics_Workshop/gsea-3.0.jar      

```
Prepare files required to run GSEA  
For detailed file format, see [here]( https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29)

**Files**| **Format**
---|---
Expression data files| Gene cluster text file format (.gct)
Gene set files|gene matrix file format(.gmx), gene matrix transposed format(.gmt)
Phenotype data files|Categorical class file (.cls) (defining experimental group)

A set of following sample files are prepared for GSEA analysis practice, which are located at  `/projects/oarc/Genomics_Workshop/GSEA/`
```
fpkm_values.ready.gct    gct file (expression fpkm values)
fpkm_values.ready.cls     cls file  (defining experimental group)
Mouse_Human_NCI_Nature_November_01_2017_symbol.gmt     gmt file (gene set file biological function set)

```
 In practice analysis, use online broad C2 geneset instead of the above .gmt file  
 Results are located at `/home/Net_ID/gsea_home/output/<jun27>/my_analysis.Gsea.nnnnnnnnnnnnn/ `  

  To view your result:
```
  cd ~/gsea_home/output/<jun27>/my_analysis.Gsea.nnnnnnnnnnnnn/  
  ## (same as: cd /home/Net_ID/gsea_home/output/<jun27>/my_analysis.Gsea.nnnnnnnnnnn/)
  firefox index.html
```
## Additional gene set database downloading source:

  * <http://software.broadinstitute.org/gsea/msigdb/index.jsp>
  * <http://download.baderlab.org/EM_Genesets/>
  * <http://www.go2msig.org/cgi-bin/prebuilt.cgi>
  * or build your [own](http://www.go2msig.org/cgi-bin/go2msig.cgi)


# Lab 4 Running the GO term analysis
___
Open [this](http://amigo.geneontology.org/amigo) link  
Gene list from de-analysis of our downloaded data (selected based on FDR and FC):

```
FDR<0.05, FC < -2.5     FDR<0.05, FC < -2.5     FDR<0.05, FC > 2.5      FDR<0.05, FC > 2.5
ENSG00000146006 ENSG00000123610 ENSG00000139220 ENSG00000235927
ENSG00000108700 ENSG00000124766 ENSG00000136478 ENSG00000099860B 
ENSG00000162692 ENSG00000176771 ENSG00000112218 ENSG00000153904
ENSG00000105989 ENSG00000196517 ENSG00000157510 ENSG00000243244
ENSG00000188176 ENSG00000132622 ENSG00000071282 ENSG00000169218
ENSG00000141469 ENSG00000126016 ENSG00000171385 ENSG00000163513
ENSG00000116991 ENSG00000128342 ENSG00000108960 ENSG00000187498
ENSG00000119714 ENSG00000116711 ENSG00000116962 ENSG00000148175
ENSG00000214814 ENSG00000025423 ENSG00000111859 ENSG00000108924
ENSG00000126878 ENSG00000125848 ENSG00000145675 ENSG00000180139
ENSG00000172061 ENSG00000163394 ENSG00000140511 ENSG00000245812
ENSG00000184564 ENSG00000272841 ENSG00000110756 ENSG00000158813
ENSG00000122877 ENSG00000181634 ENSG00000162616 ENSG00000068383
ENSG00000131771 ENSG00000243742 ENSG00000278621 ENSG00000221869
ENSG00000165272 ENSG00000103742 ENSG00000135678 ENSG00000213626
ENSG00000145777 ENSG00000172497 ENSG00000241399 ENSG00000149591
ENSG00000013293 ENSG00000254726 ENSG00000165507 ENSG00000131386
ENSG00000146250 ENSG00000131389 ENSG00000267669 ENSG00000164442
ENSG00000143494 ENSG00000016391 ENSG00000179862 ENSG00000261490
ENSG00000154864 ENSG00000157368 ENSG00000147119 ENSG00000072571
ENSG00000163491 ENSG00000099194 ENSG00000134121 ENSG00000156675
ENSG00000183508 ENSG00000049246 ENSG00000168621 ENSG00000171793
ENSG00000128165 ENSG00000028277 ENSG00000048540 ENSG00000174437
ENSG00000123689 ENSG00000107562 ENSG00000133142 ENSG00000163171
ENSG00000136999 ENSG00000146592 ENSG00000179820 ENSG00000172260
ENSG00000128606 ENSG00000100784 ENSG00000175471 ENSG00000161647
ENSG00000128510 ENSG00000139269 ENSG00000151726 ENSG00000137265
ENSG00000178695 ENSG00000168398 ENSG00000135362 ENSG00000162878
ENSG00000177614 ENSG00000235109 ENSG00000162998 ENSG00000198431
ENSG00000138316 ENSG00000196932 ENSG00000106617 ENSG00000137959
ENSG00000108830 ENSG00000148848 ENSG00000035664 ENSG00000131979
FDR<0.05, FC < -2.5     FDR<0.05, FC < -2.5     FDR<0.05, FC > 2.5      FDR<0.05, FC > 2.5
ENSG00000168811 ENSG00000147883 ENSG00000270885 ENSG00000162493
ENSG00000177570 ENSG00000183876 ENSG00000146122 ENSG00000162772
ENSG00000117600 ENSG00000131242 ENSG00000172403 ENSG00000116675
ENSG00000160145 ENSG00000100302 ENSG00000164647 ENSG00000154930
ENSG00000134253 ENSG00000126950 ENSG00000137880 ENSG00000196569
ENSG00000130513 ENSG00000182010 ENSG00000103175 ENSG00000145244
ENSG00000089041 ENSG00000105516 ENSG00000167191 ENSG00000169738
ENSG00000168918 ENSG00000235513 ENSG00000169031 ENSG00000211448
ENSG00000070808 ENSG00000007237 ENSG00000154856 ENSG00000237697
ENSG00000134363 ENSG00000162643 ENSG00000163110 ENSG00000157214
ENSG00000278727 ENSG00000135472 ENSG00000142871 ENSG00000116194
ENSG00000106484 ENSG00000138669 ENSG00000213160 ENSG00000095637
ENSG00000225783 ENSG00000160097 ENSG00000280143 ENSG00000169715
ENSG00000276600 ENSG00000054938 ENSG00000100242 ENSG00000119138
ENSG00000013297 ENSG00000138135 ENSG00000197943 ENSG00000149218
ENSG00000126861 ENSG00000186198 ENSG00000280099 ENSG00000185950
ENSG00000106976 ENSG00000185745 ENSG00000128262 ENSG00000137672
ENSG00000172738 ENSG00000127824 ENSG00000100206 ENSG00000138829
ENSG00000129682 ENSG00000158806 ENSG00000246430 ENSG00000166741
ENSG00000134259 ENSG00000123612 ENSG00000184307 ENSG00000163661
ENSG00000122966 ENSG00000149256 ENSG00000259426 ENSG00000253368
ENSG00000112837 ENSG00000143786 ENSG00000081052 ENSG00000267480
ENSG00000102524 ENSG00000170989 ENSG00000070404 ENSG00000165072
ENSG00000132321 ENSG00000223949 ENSG00000137801 ENSG00000165899
ENSG00000133216 ENSG00000129467 ENSG00000154736 ENSG00000176928
ENSG00000100739 ENSG00000196155 ENSG00000119139 ENSG00000067798
ENSG00000143320 ENSG00000111728 ENSG00000127083 ENSG00000162614
ENSG00000183496 ENSG00000117461 ENSG00000108387 ENSG00000143869
ENSG00000227268 ENSG00000103647 ENSG00000137393 ENSG00000163251
ENSG00000092969 ENSG00000272168 ENSG00000174944 ENSG00000163017
ENSG00000223764 ENSG00000137872 ENSG00000170873 ENSG00000150907
ENSG00000088756 ENSG00000167992 ENSG00000139132 ENSG00000197381
ENSG00000166592 ENSG00000166793 ENSG00000185432 ENSG00000205364
ENSG00000107611 ENSG00000100292 ENSG00000134243 ENSG00000167549
ENSG00000213420 ENSG00000137266 ENSG00000261685 ENSG00000060718
ENSG00000110900 ENSG00000165891 ENSG00000158246 ENSG00000102554
ENSG00000258947 ENSG00000164619 ENSG00000122035 ENSG00000141401
ENSG00000101825 ENSG00000246763 ENSG00000119508 ENSG00000159167
ENSG00000102984 ENSG00000173114 ENSG00000140807 ENSG00000145390
ENSG00000064309 ENSG00000230417 ENSG00000125148 ENSG00000116285

FDR<0.05, FC < -2.5             FDR<0.05, FC > 2.5      FDR<0.05, FC > 2.5
ENSG00000164761         ENSG00000177283 ENSG00000268913
ENSG00000149633         ENSG00000230018 ENSG00000103196
ENSG00000012048         ENSG00000261468 ENSG00000162630
ENSG00000126860         ENSG00000197301 ENSG00000247311
ENSG00000092621         ENSG00000154734 ENSG00000046653
ENSG00000154263         ENSG00000169271 ENSG00000167641
ENSG00000079462         ENSG00000124440 ENSG00000135821
ENSG00000182580         ENSG00000099998 ENSG00000136237
ENSG00000167771         ENSG00000120162 ENSG00000099337
ENSG00000205208         ENSG00000126803 ENSG00000120129
ENSG00000172986         ENSG00000068831 ENSG00000004799
ENSG00000272341         ENSG00000123685 ENSG00000221866
                ENSG00000128045 ENSG00000157150
                ENSG00000101342 ENSG00000102760
                ENSG00000096060 ENSG00000198624
                ENSG00000128917 ENSG00000179094
                ENSG00000163083 ENSG00000179300
                ENSG00000173838 ENSG00000136383
                ENSG00000143127 ENSG00000189221
                ENSG00000163884 ENSG00000174697
                ENSG00000168309 ENSG00000112936
                ENSG00000152583 ENSG00000165995
                ENSG00000127954 ENSG00000157514
                ENSG00000250978 ENSG00000233117
                ENSG00000109906 ENSG00000157152
                ENSG00000179593 ENSG00000187193
                ENSG00000101347 ENSG00000152779
                ENSG00000211445 ENSG00000170214
```


# Lab 5 ID mapping and conversion
___
 Learn about gene identifiers, g:profiler, Synergizer and BioMart

Use the above gene list:    
1. Convert Gene IDs to Entrez Gene, gene name: Use g:Profiler 
  Explore more functions, what the site can do for you   
2. Get gene name, GO annotation + evidence codes Use Ensembl BioMart  
3. Do it again with your own gene list  

