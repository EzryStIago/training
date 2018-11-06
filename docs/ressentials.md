# To start R shell on the cluster 

See workshop for explanation of the various options of `srun`
```bash
        module load intel/17.0.4
        module load R-Project/3.4.1
        ## run an interactive shell for 1 hr 40 min on 1 node, 2 cores; you will be placed on a compute node this way
        srun -p main -N 1 -c 2 -n 1 -t 01:40:00 --pty /bin/bash
        ##start R on compute node now
        R
```
# To submit R job using sbatch:
        
        #!/bin/bash
        #SBATCH -n 1
        #SBATCH -t 1:00
        #SBATCH -p main
        #SBATCH --export=ALL
        #SBATCH -o example-%j.out
        #SBATCH -e example-%j.err

        Rscript yourRprogram.r
        

# Packages used from BioConductor

If these packages are not installed, you can install them yourself. On login node, start R and inside R copy-paste the following commands: 

```r
        source("https://bioconductor.org/biocLite.R") 
        biocLite("ape")
        biocLite("MKmisc")
        biocLite("Heatplus")
        biocLite("affycoretools")
        biocLite("flashClust")
        biocLite("affy")
```

## Example: Calculate gene length  

Get some data from ENSEMBLE  
```wget ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz```

In R shell, you can execute these commands to compute gene lengths: 
```r

         library(GenomicFeatures)
         gtfdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.78.gtf",format="gtf")
         exons.list.per.gene <- exonsBy(gtfdb,by="gene")
         exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
         class(exonic.gene.sizes)

         Hg20_geneLength <-do.call(rbind, exonic.gene.sizes)
         colnames(Hg20_geneLength) <- paste('geneLength')    
```

# Some R essentials

## Arithmetic functions
```r
        2+2
        3*3
        3*8+2
        log10(1000)
        log2(8)
        abs(-10)
        sqrt(81)
```

## Creating objects

```r
        ls()  #see what objects are in the workspace
        x <- 4
        x
        x = 3  #a single = is an assignment operator
        x
        x == 5 #a double == asks "is the left side equivalent to the right side?"
        x + 2   #objects can be used in equations
        y <- "anyname"
        y
        class(x)
        class(y)
        ls()
```


## Vector and Matrix

```r
        x1 <- c(1,2,3,4,5)
        x1
        class(x1)
        length(x1)
        x <- cbind(x1, x1+1)    #1 will be added to all the numbers in x1
        x
        class(x)       #what kind of object is x?
        dim(x)         #the dimension of matrix
        x1[1:3]        #use [] to get subsets of a vector
        x[1,]          #use [,] to get subsets of a matrix (or dataframe)
        x[,1]
        x[,-1]
        x[c(1,2),]
        x[-c(1,3),]
        colnames(x)
        colnames(x) <-c("A","B")
        rownames(x) <-c("C","D","E","F","G")
        x
```

## Data Frames

```r
        z <- data.frame(A=x[,1], B=rownames(x), C=factor(rownames(x)), D=x[,1]==3, stringsAsFactors=F)
        class(z)
        names(z)
        dim(z)
        class(z$A)
        class(z$B)
        class(z$C)
        class(z$D)
        z$B
        z$C
```

## More ways to subset dataframes


```r
        z$B
        z[[2]]
        z[,2]   #these first 3 give equivalent results
        z[,1:2]
        z[,c(1,3)]
        z[c(1,3:5),]
```


## Lists

```r
        mylist <- list(first=z,second=x,third=c("W","X","Y","Z"))
        class(mylist)
        mylist
        names(mylist)
        class(mylist$first)
        class(mylist$second)
```


## Functions

```r
        my.add <- function(a, b) {a - b}
        class(my.add)
        my.add(4,99)
        my.add(99,4)
        my.add(b = 99, a = 4)
```

## Various directory/file/library manipulations

```r
        library(limma)  #load the limma package


        #### Make sure the working directory is set to your file on the computer;

        getwd()  #see what the current working directory is
        setwd("????????????????")  #change the working directory


        #### Output a single object as a comma separated value file

        write.csv(z, file="test.csv")
```


## Save all the objects you have created to your workspace
```r
        save.image()                #creates a default file named ".RData"
        save.image("intro.Rdata")   #creates a named file
```



## Remove objects from your workspace

```r
        ls()
        rm(x)          #remove a single object by name
        ls()
        rm(z,x1)       #remove multiple objects by name
        ls()
        load("intro.Rdata")
        ls()
        rm(list=ls())  #remove all objects
        ls()
```


## Save a history of all the commands entered

```r
        savehistory("introhistory.Rhistory")
```
