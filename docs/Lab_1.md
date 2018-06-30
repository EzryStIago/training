# LAB 1: Visualization using IGV

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


 







