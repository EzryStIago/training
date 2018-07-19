       Practice Lab partII

##Getting the files ready:      

##By now, you are assumed to have completed the tophat mapping, the output alignment files shall be located at

 /scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508 (untreated_SRR1039512  untreated_SRR1039516) 

/scratch/$USER/Genomics_Workshop/dex_treated/tophat_out/dex_SRR1039509 (dex_SRR1039513, dex_SRR1039517);

You are also expected to have the read counts ready, located at:   /scratch/$USER/Genomics_Workshop/read_counts

SHORTCUT 
## To ensure the workshop progress, we have those files prepared already. 
 Please run the following code to have those files copied and linked into your folders:

	 /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/Labs/lab_PartII.sh

Start an interactive session on the compute node
       srun --x11 -p main --reservation=genomics -N 1 -c 2 -n 1 -t 2:00:00 --pty /bin/bash -i
##The above command will land you on a reserved compute node

#############################################################
Lab 1: Visualization using IGV
#############################################################   
### IGV  --- Focus on visualization, best for validation and confirmation of the analysis result,  Not good for primary analysis

###The mapping file is in bam format, located under the folder of tophat_out, they shall be sorted and indexed using the following command ### 

	cd /scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508
	module load samtools							         
	samtools sort accepted_hits.bam -o accepted_hits.sorted.bam  ##sort the bam file 
	samtools index accepted_hits.sorted.bam ## index the sorted bam

### The resulting files: accepted_hits.sorted.bam    
  			 accepted_hits.sorted.bam.bai 
				are the files to be uploaded to IGV

You need to repeat these steps for every sample 

	SHORTCUT LAB_1
	
	### We have prepared 4 sets of such files (dex09, dex13, untreated08 and untreated12), located at: /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/Bam_for_IGV/ . Make a soft link (see the following command), or copy them into your scratch  folder, then we use IGV to visualize them.
    
	cd /scratch/$USER/Genomics_Workshop/
	## ln -s /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/Bam_for_IGV  Bam_for_IGV  

###start IGV 
	module load java
	/projects/oarc/Genomics_Workshop/Programs/IGV_2.4.6/igv.sh 

###Practice and get familiar with:
	How to Load genome and data track
	How to navigate
	How and what to visualize:
	Examine coverage
	Low mapping quality
	Mis-alignment
	Translocation
	Novel genes/transcript
	Alternative splicing
	Inversion
	Look for SNPs

More detailed explanation here:  http://software.broadinstitute.org/software/igv/book/export/html/37
 
		
https://software.broadinstitute.org/software/igv/ (##DETAILED SOFTWARE DOCUMENTATION) 

CLOSE your interactive session on a node  when done with IGV by typing
	exit

#####################################################################
Lab 2:  data processing and expression analysis using edgeR
#####################################################################

#####Lab2_1). Prepare the read_counts file to be imported into R for analysis  #####

## Here we created a folder named as read_counts, where all 6 read counts files are located --- this has been done for you via lab_PartII.sh in the beginning

            read_counts/dex09.txt
			dex13.txt
			dex17.txt
			untreated08.txt
			untreated12.txt
			untreated16.txt

###go to this folder containing all above six counts output files
cd /scratch/$USER/Genomics_Workshop/read_counts
 
### We need to combine all counts into one file, which will be imported into R for further analysis.
 The following shows how this can be done with bash commands on linux ( you may do it in excel too)
	
	 paste dex09.txt dex13.txt dex17.txt untreated08.txt untreated12.txt untreated16.txt > merged_counts.txt    	## merge files by columns

	 cut -f1,2,4,6,8,10,13 merged_counts.txt  > merged_counts_clean.txt     ## extract the relevant columns  

	 head -n -5 merged_counts_clean.txt > merged_counts_clean2.txt      ##remove the last 5 line stats summary


## We also need to prepare a file containing group/treatment information. This file Targets.txt is a tab delimited text file.  You may construct in excel. 
### The file should contain the following info:
label	sample	group	treatment
1	dex09	dex	dex_treated
2	dex13	dex	dex_treated
3	dex17	dex	dex_treated
4	untreated08	control	untreated
5	untreated12	control	untreated
6	untreated16	control	untreated
			
This file has been prepared for you. You will need to Copy it into your folder, see later Lab2_2:


##### Lab2_2). pre-processing the data in R #######

	cd /scratch/$USER/Genomics_Workshop/
	mkdir  DE_analysis  ##set working directory to run differential expression analysis
	cd DE_analysis

## Copy the needed files here

	cp /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/SRA_data/DE_analysis/Targets.txt $PWD         	##This file denotes the experimental group
	cp /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/SRA_data/read_counts/merged_counts_clean2.txt $PWD    ##This file contains read counts on genes for all samples
	cp /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/SRA_data/DE_analysis/Hg20_geneLength $PWD      	### gene length file needed for FPKM calculation (This file shall be generated from lab_PartI, lab10)  


#### From your homework assignment you should have the following packages to be installed already :).
If you didn't, install them now following instruction in lab_PartI.txt, lab10
########## 

##now,start a new interactive job on the compute node
or  switch to another terminal if you still have an interactive job running)
	srun --x11 -p main --reservation=genomics -N 1 -c 2 -n 1 -t 02:00:00 --pty /bin/bash -i
# go to your working directory
	cd /scratch/$USER/Genomics_Workshop/DE_analysis/
##then start R by doing the following:
	module load intel/17.0.4
	module load R-Project/3.4.1   
	R 
	
	>getwd()  ##check which directory you are in
## you should be in the directory where the files merged_counts_clean2.txt and "Targets.txt", annotation gtf file are.
If not there,  set your working directory :

        >setwd("/scratch/$USER/Genomics_Workshop/DE_analysis")

##import count data to R 
> raw.data <- read.delim("merged_counts_clean2.txt", header=F, row.names=1)
> head(raw.data)  ## check the beginning portion of the imported file, now an object
> dim(raw.data)   ##check the dimention of this object
> class(raw.data)  ## check what class of this object
> apply(raw.data, 2, summary) ## check the range of counts per sample
> range(rowSums(raw.data))      ## check the range of counts per gene
##add column header 
> colnames(raw.data) <- c("dex09","dex13","dex17","untreated08","untreated12","untreated16")

>raw.data2 <- raw.data[rowSums(raw.data) != 0 ,] ##remove genes with zero counts for all samples
> dim(raw.data2)

> cpm.values <- cpm(raw.data2) #calculate counts per million mapped reads without any other normalization

>above1cpm <- apply(cpm.values, 1, function(x) sum(x >=1)) ##How many samples/genes had at least 1 cpm
> counts.use <- raw.data2[above1cpm >= 3,] ##we have three replicates in each group. If a gene can be reliably detected, it should be detected in all 3 replicates. So, a gene to be included for further analysis shall have 1 cpm in at least 3 samples.

>dim(counts.use)
>colSums(counts.use) / colSums(raw.data)  ##the % of total counts kept after filtering
>nrow(counts.use) / nrow(raw.data)  ##the % of genes kept after filtering

> targets <- readTargets()      ##import targets file that contains experiment group info.
> targets$GpF <- factor(targets$group)   ##change character to factor
> targets$GpN <- as.numeric(targets$GpF) ##change factor to numeric (optional)

###density distribution####
> plotDensity(log2(raw.data+0.1),col=targets$GpF,lty=1,lwd=2)
> legend("topright", legend=levels(targets$GpF),fill=1:4)

##change data from raw.data to raw.data2, to CPM, FPKM,.. to see the effect of filtering and normalization (When those values generated)


> ls()  #see that objects have been loaded
>save.image("RNASeqDemo.RData")   ##save the work workspace
>savehistory("RNASeqDemo.Rhistory")  ##save command history


##### Lab2_3). To calculate expression values as fpkm #######

##First, compute the gene length as described in lab10--partI, the code is pasted as following: (This is the same as lab_PartI.txt, lab10, and the result file Hg20_geneLength is provided already, you may import it to R:  

> Hg20_geneLength <- read.delim(file="Hg20_geneLength", header=T, row.names=1)

				> library(GenomicFeatures)  ##may skip, because we already loaded it at the start
				> gtfdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.78.gtf",format="gtf")
				> exons.list.per.gene <- exonsBy(gtfdb,by="gene")
				> exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
				> class(exonic.gene.sizes)
				>Hg20_geneLength <-do.call(rbind, exonic.gene.sizes)
				> colnames(Hg20_geneLength) <- paste('geneLength')

>Hg20_geneLength2 <- data.frame(Hg20_geneLength[rownames(counts.use),]) ## to extract the gene length file to contain the same number genes in the same order as in the filtered read counts file
>colnames(Hg20_geneLength2) <- paste('geneLength')  ## to change column name, make it neat
> fpkm.data <- cpm(counts.use) / (Hg20_geneLength2$geneLength/1000) ## compute fpkm
>  min.count <- apply(fpkm.data,1,min)
> write.csv(fpkm.data,file="fpkm_values.csv")  #### To output FPKM data


##### Lab2_4). differential expression analysis  with edgeR #######

#### Create the design matrix
> groups <- factor(targets$group)
> design <- model.matrix(~0+groups)
> colnames(design) <- levels(groups)
> rownames(design) <- targets$sample
> design
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


 #### Create contrast matrix
>cont.matrix <- makeContrasts(dex_Ctl= dex - control, levels=design)


###Create DGEList object

>  d <- DGEList(counts=counts.use, lib.size=colSums(counts.use), group=targets$GpF)
> class(d)
> names(d)    ## the names of the items in the list
> d$counts[1:5,]    ## check the counts which are stored in the $counts:
> d$samples   ## The group info and library sizes stored in $samples
> d <- calcNormFactors(d)  ## an additional normalization factor using a TMM method
> d$samples   ## now the norm.factors are no longer 1

###this can be useful when diagnose problem)###
> apply(d$counts,2,function(x) sum(sort(x/sum(x),decreasing=T)[1:20])) * 100
##the proportions of total counts for the top 20 genes in each sample, 10-20% is OK.

>logCPM <-  cpm(d$counts, prior.count=0.01,  log = TRUE)   ## modified logCPM values in edgeR, can be used for clustering, heatmap (prior.count: average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.)

#####Now,the DE test! The term "estimating dispersions" in edgeR describes a method to account for the variance among replicates.

> d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
Disp = 0.06037 , BCV = 0.2457
> d <- estimateGLMTrendedDisp(d, design)
> d <- estimateGLMTagwiseDisp(d, design)
> plotBCV(d)   ## the relationship between the overall abundance and the tagwise dispersion estimates

> fit.edgeR <- glmFit(d, design)  ## Estimate model coefficients from count data and design matrix
> names(fit.edgeR)

####Specify contrasts of interest,  do empirical Bayes "shrinkage" of  variances and calculate test statistics. Both of these are performed with same function glmLRT in edgeR (Genewise Negative Binomial Generalized Linear Models)
 > eR.dex_Ctl <- glmLRT(fit.edgeR, contrast=cont.matrix[,1])

> summary(decideTestsDGE(eR.dex_Ctl)) [,1]  ## Correct for multiple tests and extract relvant data (1 means sig up, -1 means sig down, and 0 means NS)

> eR.dex_Ctl.detailed <- topTags(eR.dex_Ctl,n=Inf,sort.by="none")$table ## Get detailed output for a single contrast


>  eR.dex_Ctl.detailed[1:5,]
                      logFC   logCPM         LR    PValue       FDR
ENSG00000000003 -0.36354885 5.227841 2.37143785 0.1235732 0.5501740
ENSG00000000419  0.18626588 4.674100 0.79016266 0.3740509 0.8571700
>  eR.dex_Ctl.detailed[1:5,]
                      logFC   logCPM         LR    PValue       FDR
ENSG00000000003 -0.36354885 5.227841 2.37143785 0.1235732 0.5501740
ENSG00000000419  0.18626588 4.674100 0.79016266 0.3740509 0.8571700
ENSG00000000457  0.04304803 3.838702 0.03317218 0.8554789 0.9905954
>  eR.dex_Ctl.detailed[1:5,]
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
> eR.dex_Ctl.detailed$FC <- 2^abs(eR.dex_Ctl.detailed$logFC) * sign(eR.dex_Ctl.detailed$logFC)

> write.csv(eR.dex_Ctl.detailed,file="Demo_eR.dex_Ctl_results.csv")  ##output the file


##### Lab2_5). Analysis QC ---or sample diagnosis #######

######Clustering#####
hc.raw <- flashClust(dist(t(log2(raw.data2+0.1))),method="average")
plot(hc.raw,hang = -1, main = "RNASeqDemo, raw values", sub = "", xlab = "",cex=0.9, labels=targets$sample)
##change data from raw.data to raw.data2, to CPM, logCPM, FPKM,.. to see the effect of filtering and normalization

#####PCA#######
>  plotPCA(fpkm.data, pch=16, col=targets$GpF,groupnames=levels(targets$GpF), addtext=rep(1:3,4),main="PCA on FPKM")

###do it after edgeR analysis, otherwise some values not existing yet###

>  plotPCA(logCPM), pch=16, col=targets$GpF,groupnames=levels(targets$GpF), addtext=rep(1:3,4),main="PCA on logCPM")


#####heatmap####
>test <- topTags(eR.dex_Ctl,n=Inf,sort.by="PValue")$table ####sort gene list according to P values
>test2 <- test[1:500,]  ###Take the top 500 significant genes
>logCPM2 <- logCPM[rownames(test2),]  ###Extract logCPM values of these 500 genes
> meanFC <- logCPM2 - rowMeans(logCPM2)
> color.meanFC <- heatmapCol(data = meanFC, lim = 3, col =colorRampPalette(c("blue","white","red"))(128))
> heatmap_2(meanFC, col=color.meanFC, legend=3, scale="none")
 
 
 #####Volcano plot####
> with(eR.dex_Ctl.detailed, plot(logFC, -log10(PValue), pch=20, main="Volcano plot", xlim=c(-2.5,2),ylim=c(0,25)))
> with(subset(eR.dex_Ctl.detailed,FDR<0.05), points(logFC, -log10(PValue), pch=20, col="red"))
> with(subset(eR.dex_Ctl.detailed, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
> with(subset(eR.dex_Ctl.detailed, FDR<0.05&abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))



#############################################################
Lab 3:  Running enrichment analysis using GSEA 
#############################################################
#####The command to start the gsea:
	srun --x11 -p main --reservation=genomics -N 1 -c 2 -n 1 -t 2:00:00 --pty /bin/bash -i  ##get onto a reserved compute node
	module load java
	java -jar /scratch/$USER/Genomics_Workshop/gsea-3.0.jar      

##Prepare files required to run GSEA
##For detailed file format: see https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29

Expression data files	
	Gene cluster text file format (.gct)

Gene set files
		gene matrix file format  (.gmx)
		 
gene matrix transposed format (.gmt)


Phenotype data files
      		Categorical class file (.cls)  (defining experimental group)


## a set of following sample files are prearaed for GSEA analysis practice, which are located at /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/GSEA_analysis/

	fpkm_values.ready.gct    gct file (expression fpkm values)
	fpkm_values.ready.cls     cls file  (defining experimental group)
	Human_GO_AllPathways_no_GO_iea_July_01_2018_symbol.gmt     gmt file (gene set file biological function set)

** in practice analysis, use online broad C2 geneset instead of the above .gmt file
** Results are located at (by default)  /home/$USER/gsea_home/output/<july19>/my_analysis.Gsea.nnnnnnnnnnnnn/, but you can set a location to name your own output, we will practice this in class, lets set it at your GSEA_analysis folder 
  ## To view your result:
  ## go to the folder where you set your GSEA output file, here assume your GSEA_analysis folder,  then
  firefox index.html
  ##if you did not set your GSEA output folder, then you need to go to its default location to view the file:
  cd ~/gsea_home/output/<july19>/my_analysis.Gsea.nnnnnnnnnnnnn/  ## (same as: cd /home/$USER/gsea_home/output/<july19>/my_analysis.Gsea.nnnnnnnnnnn/)
  
** Additional gene set database downloading source:
	
 http://software.broadinstitute.org/gsea/msigdb/index.jsp
		http://download.baderlab.org/EM_Genesets/
		http://www.go2msig.org/cgi-bin/prebuilt.cgi
		or build your own:
			http://www.go2msig.org/cgi-bin/go2msig.cgi
	
#############################################################
Lab 4: Running the GO term analysis   
#############################################################

        cd /scratch/$USER/Genomics_Workshop/
	mkdir GO_analysis
	cd GO_analysis

###start Cytoscape
        module load java
        /projects/oarc/Genomics_Workshop/Programs/Cytoscape_v3.6.1/cytoscape.sh


Gene list from de-analysis of our downloaded data (selected based on FDR and FC):

FDR<0.05, FC > 2		
ENSG00000179593	ENSG00000247311	ENSG00000165072	ENSG00000280099	ENSG00000108960	ENSG00000025708	ENSG00000097096	ENSG00000113460
ENSG00000109906	ENSG00000162630	ENSG00000267480	ENSG00000197943	ENSG00000171385	ENSG00000169750	ENSG00000102804	ENSG00000164125
ENSG00000250978	ENSG00000103196	ENSG00000253368	ENSG00000100242	ENSG00000071282	ENSG00000156804	ENSG00000118689	ENSG00000196507
ENSG00000127954	ENSG00000268913	ENSG00000163661	ENSG00000280143	ENSG00000157510	ENSG00000163378	ENSG00000077943	ENSG00000182552
ENSG00000152583	ENSG00000126803	ENSG00000166741	ENSG00000213160	ENSG00000112218	ENSG00000275482	ENSG00000256235	ENSG00000057657
ENSG00000168309	ENSG00000120162	ENSG00000138829	ENSG00000142871	ENSG00000136478	ENSG00000147576	ENSG00000163431	ENSG00000081320
ENSG00000163884	ENSG00000099998	ENSG00000137672	ENSG00000163110	ENSG00000139220	ENSG00000166979	ENSG00000196616	ENSG00000185813
ENSG00000143127	ENSG00000124440	ENSG00000185950	ENSG00000154856	ENSG00000131979	ENSG00000085117	ENSG00000170485	ENSG00000083223
ENSG00000173838	ENSG00000169271	ENSG00000149218	ENSG00000169031	ENSG00000137959	ENSG00000138166	ENSG00000077684	ENSG00000164104
ENSG00000163083	ENSG00000154734	ENSG00000119138	ENSG00000167191	ENSG00000198431	ENSG00000136546	ENSG00000168994	ENSG00000148120
ENSG00000128917	ENSG00000197301	ENSG00000169715	ENSG00000103175	ENSG00000162878	ENSG00000206538	ENSG00000138678	ENSG00000078053
ENSG00000096060	ENSG00000261468	ENSG00000095637	ENSG00000137880	ENSG00000137265	ENSG00000138356	ENSG00000008256	ENSG00000174306
ENSG00000101342	ENSG00000230018	ENSG00000116194	ENSG00000164647	ENSG00000161647	ENSG00000123358	ENSG00000155324	ENSG00000185112
ENSG00000128045	ENSG00000177283	ENSG00000157214	ENSG00000172403	ENSG00000172260	ENSG00000254254	ENSG00000107796	ENSG00000132970
ENSG00000211445	ENSG00000116285	ENSG00000237697	ENSG00000146122	ENSG00000163171	ENSG00000072163	ENSG00000120899	ENSG00000177674
ENSG00000101347	ENSG00000145390	ENSG00000211448	ENSG00000270885	ENSG00000174437	ENSG00000164105	ENSG00000185022	ENSG00000166825
ENSG00000170214	ENSG00000159167	ENSG00000169738	ENSG00000035664	ENSG00000171793	ENSG00000163803	ENSG00000183044	ENSG00000145386
ENSG00000152779	ENSG00000141401	ENSG00000145244	ENSG00000106617	ENSG00000156675	ENSG00000206190	ENSG00000138073	ENSG00000034152
ENSG00000187193	ENSG00000102554	ENSG00000196569	ENSG00000162998	ENSG00000072571	ENSG00000162426	ENSG00000168646	ENSG00000083290
ENSG00000157152	ENSG00000060718	ENSG00000154930	ENSG00000135362	ENSG00000261490	ENSG00000138615	ENSG00000196975	ENSG00000119711
ENSG00000233117	ENSG00000167549	ENSG00000116675	ENSG00000151726	ENSG00000164442	ENSG00000023909	ENSG00000137124	ENSG00000009413
ENSG00000157514	ENSG00000205364	ENSG00000162772	ENSG00000175471	ENSG00000131386	ENSG00000114098	ENSG00000020577	ENSG00000143878
ENSG00000165995	ENSG00000197381	ENSG00000162493	ENSG00000179820	ENSG00000149591	ENSG00000240859	ENSG00000114270	ENSG00000008311
ENSG00000112936	ENSG00000125148	ENSG00000123685	ENSG00000133142	ENSG00000213626	ENSG00000099204	ENSG00000134294	ENSG00000123562
ENSG00000174697	ENSG00000140807	ENSG00000068831	ENSG00000048540	ENSG00000221869	ENSG00000138074	ENSG00000167645	ENSG00000214944
ENSG00000189221	ENSG00000119508	ENSG00000170873	ENSG00000168621	ENSG00000068383	ENSG00000128311	ENSG00000226950	ENSG00000107968
ENSG00000136383	ENSG00000122035	ENSG00000174944	ENSG00000134121	ENSG00000158813	ENSG00000084090	ENSG00000133816	ENSG00000128923
ENSG00000179300	ENSG00000158246	ENSG00000137393	ENSG00000147119	ENSG00000245812	ENSG00000131459	ENSG00000159212	ENSG00000003402
ENSG00000179094	ENSG00000261685	ENSG00000108387	ENSG00000179862	ENSG00000180139	ENSG00000069535	ENSG00000154127	ENSG00000163697
ENSG00000198624	ENSG00000134243	ENSG00000127083	ENSG00000267669	ENSG00000108924	ENSG00000106123	ENSG00000169908	ENSG00000137962
ENSG00000102760	ENSG00000185432	ENSG00000119139	ENSG00000165507	ENSG00000148175	ENSG00000108604	ENSG00000151690	ENSG00000108950
ENSG00000157150	ENSG00000139132	ENSG00000154736	ENSG00000241399	ENSG00000187498	ENSG00000173918	ENSG00000196850	ENSG00000158079
ENSG00000221866	ENSG00000150907	ENSG00000137801	ENSG00000135678	ENSG00000163513	ENSG00000124151	ENSG00000147027	ENSG00000074590
ENSG00000004799	ENSG00000163017	ENSG00000070404	ENSG00000278621	ENSG00000169218	ENSG00000168556	ENSG00000272870	ENSG00000144362
ENSG00000120129	ENSG00000163251	ENSG00000081052	ENSG00000162616	ENSG00000243244	ENSG00000130066	ENSG00000165030	ENSG00000103064
ENSG00000099337	ENSG00000143869	ENSG00000259426	ENSG00000110756	ENSG00000153904	ENSG00000153207	ENSG00000269728	
ENSG00000136237	ENSG00000162614	ENSG00000184307	ENSG00000140511	ENSG00000099860	ENSG00000102996	ENSG00000136436	
ENSG00000135821	ENSG00000067798	ENSG00000246430	ENSG00000145675	ENSG00000235927	ENSG00000181061	ENSG00000122641	
ENSG00000167641	ENSG00000176928	ENSG00000100206	ENSG00000111859	ENSG00000186314	ENSG00000134686	ENSG00000140545	
ENSG00000046653	ENSG00000165899	ENSG00000128262	ENSG00000116962	ENSG00000115419	ENSG00000158716	ENSG00000141298	

The above gene ensembleID was converted to gene symbol  using BioMart, the resulting gene list in gene symbol is the file  /projects/oarc/Genomics_Workshop/RNA-Seq_analysis/GO_analysis/up_regulated_gene_FC2_FDR0.05_symbol
The ontology file "go-basic.obo" is downloaded from http://geneontology.org/page/download-ontology
			 wget http://purl.obolibrary.org/obo/go/go-basic.obo

We will first install apps (BINGO, yFiles layout)
Then test run GO enrichment analysis using the gene list in the above file "up_regulated_gene_FC2_FDR0.05_symbol"

#############################################################
Lab 5: ID mapping and conversion
#############################################################

  ----Learn about gene identifiers, g:profiler, Synergizer and BioMart
Use the above gene list:
1)	Convert Gene IDs to Entrez Gene, gene name: Use g:Profiler
	Explore more functions, what the site can do for you
2)	Get gene name, GO annotation + evidence codes Use Ensembl BioMart
3)	Do it again with your own gene list

