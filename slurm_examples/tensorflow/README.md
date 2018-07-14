# Tensorflow Example

This directory contains examples of how to train your neural networks with tensorflow using the MNIST dataset
  * For a more detailed explaination check out the tensorflow [website](https://www.tensorflow.org/versions/r1.0/get_started/mnist/beginners)
  * Make sure to clone the github [repo](https://github.com/tensorflow/models/tree/r1.8.0) for this model to a directory on the cluster and switch the branch to the version required, to ensure you have the correct packages or copy paste the commands below 
  ```
  git clone https://github.com/tensorflow/models.git

  git branch -a #This will list all branches that exist. Any that are prefixed with origin/ are on the server and you will need to fetch them.

  git checkout r1.8.0 #In this example we will be working with version 1.8.0
  git pull origin r1.8.0 #pulls any changes
  ```
  * For further reference, we have also included a benchmark that shows how the job would scale by changing the number of CPU cores and with or without GPU.
  * All training/jobs were ran for a standard of 40 epochs.

## Understanding the Code 

  * There are slight nuances between the ``mnist_1.8.0.py`` and ``mnist_1.7.0.py`` scripts.
  * Each code block is explained with comments detailing what is going on

## Understanding the Script

  * The ``tf_script_broadwell_4.sh`` trains the network on CPU. You can ask to be allocated more or less cores by specifying ``#SBATCH -c 8 # Number of cores``
  * The ``tf_script_p100.sh`` and ``tf_script_m6000`` trains the network on GPU. The difference being that one uses the P100, while the other utilizes the M6000 GPU.
  * You can change the number of GPUS allocated by specifying ``#SBATCH --gres=gpu:2       # Number of GPUs``





