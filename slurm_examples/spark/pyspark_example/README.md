# Pyspark batch example

This is an example of how to run a pyspark job as a batch script. It is doubly a batch job:  
  - it is a **slurm** batch job
  - it is a **spark** batch job  (**inside** the wrapper slurm batch job)

## What the code does

The main logic is written in `run_pyspark.py`. Look at the source code for comments.  
  - It has an input file, `iris.csv` - one of the most famous multi-label datasets with 50 measurements each of the 3 types of iris flower.
  - It calculates the average sepal lenth per group (one mean for each of the 3 groups). 
  - It has an output file, `grouped_iris.csv`, recording the results of the computation.  
  - `slurm.*.out` and `slurm.*.err` are captured console outputs for this job. 

## Run the code

To run this example:  modify the `WORKING_DIRECTORY` variable to your path (where you place the example) and execute this command in your terminal: 
```
sbatch slurm_run_pyspark.sh
```

Note: `output_data` directory should be empty (although, there is an overwrite mode line, so it should work nevertheless - When I tried it, I didn't erase the hidden files of the type `.part-...`, so Spark didn't want to overwrite otherwise empty directory despite mode overwrite!). 
