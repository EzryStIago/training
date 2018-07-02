# Lab 3  Running enrichment analysis using GSEA 

The command to start the gsea:
```
srun --x11 -p main --reservation=genomics_2 -N 1 -c 2 -n 1 -t 2:00:00 --pty /bin/bash -i  ##get onto a reserved compute node
        module load java
        java -jar /scratch/$USER/Genomics_Workshop/gsea-3.0.jar      

```
Prepare files required to run GSEA  
For detailed file format, see [here]( https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29)

>Expression data files   
>>Gene cluster text file format (.gct)

>Gene set files
>>gene matrix file format  (.gmx)
>>gene matrix transposed format (.gmt)

>Phenotype data files
>> Categorical class file (.cls)  (defining experimental group)

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
** Additional gene set database downloading source:**

  * <http://software.broadinstitute.org/gsea/msigdb/index.jsp>
  * <http://download.baderlab.org/EM_Genesets/>
  * <http://www.go2msig.org/cgi-bin/prebuilt.cgi>
  * or build your [own](http://www.go2msig.org/cgi-bin/go2msig.cgi)

                        
        


