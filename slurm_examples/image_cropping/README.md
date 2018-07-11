## Array job example

This directory contains an example of how to crop specified patches from an image (the famous `lena.png`). 

- It uses a program called ImageMagick for cropping. ImageMagick does much more than cropping. See [this link](https://www.imagemagick.org/) for details. 
- It uses a list of coordinates as input to the slurm job. This is the file `list_of_coordinates.txt`. An example is 50x50+100+0, and the interpretation is, this is a 50x50 pixel patch, whose upper left corner starts at x=100px and y=0px. An example of a command that would crop lena.png to the above patch is:  
```convert lena.png -crop 50x50+100+0 lena_cropped.png```
- It uses a slurm array job to distribute each patch crop to a different core and possibly different node, for processing in parallel on the cluster

To run this job:

- produce your own `list_of_coordinates.txt` and modify the `slurm_image_cropping.sh` script's 
    + `WORKING_DIRECTORY` 
    + `lena.png`  
    + `#SBATCH --array=0-5` to your number of patches (I had 6 image patches in `list_of_coordinates.txt`). 
- execute the following command in a bash terminal: 
```sbatch slurm_image_cropping.sh```

