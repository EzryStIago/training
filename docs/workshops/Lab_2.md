# LAB 2:  data processing and expression analysis using edgeR

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

We also need to prepare a file containing group/treatment information. This file Targets.txt is a tab delimited text file.  You may construct in excel. 
 **The file should contain the following info**

 label  | sample |group | treatment
------  |------- |------|--------------
1       |dex09 |  dex   | dex_treated
2       |dex13 | dex    |dex_treated
3       |dex17 | dex    |dex_treated
4       |untreated08    ||control untreated
5       |untreated12    ||control untreated
6       |untreated16    ||control untreated

This file has been prepared for you. You will need to Copy it into your folder, see later Lab2_2:
                        
## Lab2.2: Pre-processing the data in R

```
cd /scratch/$USER/Genomics_Workshop/
        mkdir  DE_analysis  ##set working directory to run differential expression analysis
        cd DE_analysis
## Copy the needed files here

        cp /projects/oarc/Genomics_Workshop/SRA_data/DE_analysis/Targets.txt $PWD                           ##This file denotes the experimental group
        cp /projects/oarc/Genomics_Workshop/SRA_data/read_counts/merged_counts_clean2.txt $PWD              ##This file contains read counts on genes for all samples
        cp /projects/oarc/Genomics_Workshop/Reference/hg20/Homo_sapiens.GRCh38.78.gtf $PWD            ###annotation file needed to calculate exonic gene length -- needed for FPKM calculation
```
##  HOMEWORK CATCH-UP
 From your homework assignment you should have the following packages to be installed already :).
If you didn't, install them now.

```
Open a terminal for  amarel2 login node.

ssh -X amarel2.hpc.rutgers.edu
On the login node start R
        module load intel/17.0.4
        module load R-Project/3.4.1
        R

 #then in the R workspace do the following:
        source("https://bioconductor.org/biocLite.R")
        biocLite("MKmisc")
        Would you like to use a personal library instead?  (y/n)  y

        Would you like to create a personal library
        ~/R/x86_64-pc-linux-gnu-library/3.4
        to install packages into?  (y/n) y

 Wait till it finishes.

         biocLite("Heatplus")
         biocLite("affycoretools")
         biocLite("flashClust")
         biocLite("affy")
         biocLite("GenomicFeatures")
         quit()   ###quit R, no save
Save workspace image? [y/n/c]: n
```

## Starting the Job

Now,start a new interactive job on the compute node or switch to another terminal if you still have an interactive job running)  
`srun --x11 -p main --reservation=genomics_2 -N 1 -c 2 -n 1 -t 02:00:00 --pty /bin/bash -i`  
Go to your working directory
 `cd /scratch/$USER/Genomics_Workshop/DE_analysis/`
 Then, start R
 ```
 module load intel/17.0.4
 module load R-Project/3.4.1 
 ```  
 `>getwd()` #check which directory you are in  
 You should be in the directory where the files merged_counts_clean2.txt and "Targets.txt", annotation gtf file are. If not there, set your working directory :
 `>setwd("/scratch/<your_netID>/Genomics_Workshop/DE_analysis/")` ##set the working directory

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
 ## density distribution##
  plotDensity(log2(raw.data+0.1),col=targets$GpF,lty=1,lwd=2)
  legend("topright", legend=levels(targets$GpF),fill=1:4)
```

change data from raw.data to raw.data2, to CPM, FPKM,.. to see the effect of filtering and normalization

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
Specify contrasts of interest,  do empirical Bayes "shrinkage" of  variances and calculate test statistics. Both of these are performed with same function glmLRT in edgeR (Genewise Negative Binomial Generalized Linear Models)

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

    




