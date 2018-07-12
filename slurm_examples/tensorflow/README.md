# Tensorflow Example
___
This directory contains examples of how to train your neural networks with tensorflow using the MNIST and CIFAR-10 dataset
  * For a more detailed explaination check out the tensorflow [website](https://www.tensorflow.org/versions/r1.0/get_started/mnist/beginners)
  * Make sure to clone the github [repo](https://github.com/tensorflow/models/tree/r1.8.0) for this model to a directory on the cluster and switch the branch to the version required. Copy paste the commands below 
  ```
  git clone https://github.com/tensorflow/models.git

  git branch -a #This will list all branches that exist. Any that are prefixed with origin/ are on the server and you will need to fetch them.

  git checkout r1.8.0 #In this example we will be working with version 1.8.0
  git pull origin r1.8.0 #pulls any changes
  ```
  * For further reference, we have also included a benchmark that shows how the job would scale by changing the number of CPU cores and with or without GPU.
  * All training/jobs were ran for a standard of 40 epochs.

## Understanding the Code and Job Script
___
  * There are slight nuances between the ``mnist_1.8.py`` and ``mnist-1.7.0.py`` scripts.
  * Each code block is explained with comments detailing what is going on

## Understanding the Script
___
  * The ``tf_script_broadwell_4.sh`` trains the network on CPU. You can ask to be allocated more or less cores by specifying ``#SBATCH -c 8 # Number of cores``
  * The ``tf_script_p100.sh`` trains the netwokr on GPU. You can chnage the number of GPUS allocated by specifying ``#SBATCH --gres=gpu:2       # Number of GPUs``





