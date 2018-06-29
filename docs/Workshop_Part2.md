# Practice Lab Part II

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

 then in the R workspace do the following:
        > source("https://bioconductor.org/biocLite.R")
        > biocLite("MKmisc")
        Would you like to use a personal library instead?  (y/n)  y

        Would you like to create a personal library
        ~/R/x86_64-pc-linux-gnu-library/3.4
        to install packages into?  (y/n) y

 Wait till it finishes.

        > biocLite("Heatplus")
        > biocLite("affycoretools")
        > biocLite("flashClust")
        > biocLite("affy")
        > biocLite("GenomicFeatures")
        > quit()   ###quit R, no save
Save workspace image? [y/n/c]: n
```

## Getting the files ready:

By now, you are assumed to have completed the tophat mapping, the output alignment files shall be located at
 `/scratch/$USER/Genomics_Workshop/untreated/tophat_out/untreated_SRR1039508`(untreated_SRR1039512  untreated_SRR1039516)  
 `/scratch/$USER/Genomics_Workshop/dex_treated/tophat_out/dex_SRR1039509` (dex_SRR1039513, dex_SRR1039517);  
 You are also expected to have the read counts ready, located at:
`/scratch/$USER/Genomics_Workshop/read_counts`  

**Shortcut**  
To ensure the workshop progress, we have those files prepared already. 
 Please run the following code to have those files copied and linked into your folders:  
 `/projects/oarc/Genomics_Workshop/misc/lab_PartII.sh`  
 Start an interactive session on the compute node  
 ` srun --x11 -p main --reservation=genomics -N 1 -c 2 -n 1 -t 2:00:00 --pty /bin/bash -i`  
 The above command will land you on a reserved compute node

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

 then in the R workspace do the following:
        > source("https://bioconductor.org/biocLite.R")
        > biocLite("MKmisc")
        Would you like to use a personal library instead?  (y/n)  y

        Would you like to create a personal library
        ~/R/x86_64-pc-linux-gnu-library/3.4
        to install packages into?  (y/n) y

 Wait till it finishes.

        > biocLite("Heatplus")
        > biocLite("affycoretools")
        > biocLite("flashClust")
        > biocLite("affy")
        > biocLite("GenomicFeatures")
        > quit()   ###quit R, no save
Save workspace image? [y/n/c]: n
```

## Starting the Job

now,start a new interactive job on the compute node or switch to another terminal if you still have an interactive job running)
`srun --x11 -p main --reservation=genomics_2 -N 1 -c 2 -n 1 -t 02:00:00 --pty /bin/bash -i`
Go to your working directory
 `cd /scratch/$USER/Genomics_Workshop/DE_analysis/`
 Then, start R
 ```
 module load intel/17.0.4
 module load R-Project/3.4.1 
```
 `>getwd()` #check which directory you are in  
 You should be in the directory where the files merged_counts_clean2.txt and "Targets.txt", annotation gtf fileare.
If not there,  set your working directory :











