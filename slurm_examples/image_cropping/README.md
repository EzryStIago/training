## Array job example


This directory contains an example of how to crop specified patches from an image (the famous `lena.png`). 

- It uses a program called ImageMagick for cropping. ImageMagick does much more than cropping. See [this link](https://www.imagemagick.org/) for details. 
- It uses a list of coordinates as input to the slurm job. This is the file `list_of_coordinates.txt`. An example is 50x50+100+0, and the interpretation is, this is a 50x50 pixel patch, whose upper left corner starts at x=100px and y=0px. An example of a command that would crop lena.png to the above patch is:  
```convert lena.png -crop 50x50+100+0 lena_cropped.png```
- It uses a slurm array job to distribute each patch crop to a different core and possibly different node, for processing in parallel on the cluster

## To run this job:

- produce your own `list_of_coordinates.txt` and modify the `slurm_image_cropping.sh` script's 
    + (make sure you have `results` directory, otherwise therere will be errors)
    + `WORKING_DIRECTORY` 
    + `lena.png`  
    + `#SBATCH --array=0-5` to your number of patches (I had 6 image patches in `list_of_coordinates.txt`). 
- execute the following command in a bash terminal: 
```sbatch slurm_image_cropping.sh```

## What is an array job?

We say that a task is embarassingly (or pleasantly) parallel when the computation performed in each task is the same, but the input data is different for each task. In that case it is possible to parcel out computation in parallel to speed up the entire computation. 
This is called an "array job" in Slurm. 
The way to do that is through the use of two things:  
  - an option `#SBATCH --array=0-100`  which specifies how many tasks are to be performed (in this case 101)
  - a slurm variable in the array yob `SLURM_ARRAY_TASK_ID` which keeps track of which task you are currently on. 

In this example, if the value of `SLURM_ARRAY_TASK_ID` is 10, then it will be processing 10th (actually, 11th since numbering starts with 0) image patch. 

Array jobs are just a convenient way of organizing your jobs, since the job number will be the same, let's say `10234`, but each task gets a subscript, so the individual tasks will be `10234_0`, `10234_1`, `10234_2` and so on. However, you are still launching the same number of jobs as you would if you launched them all separately, so please keep this number less than 5,000. 
