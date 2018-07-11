# Using Spark

[Apache Spark](http://spark.apache.org) is a good way to process tabular or unstructured data at scale. Spark runs usually with one of the several schedulers:  
  - Spark Standalone
  - YARN
  - Mesos
  - Kubernetes

We will describe how to spin up a virtual Spark cluster inside the cluster than runs slurm. 

## Easy - use Spark Standalone local (single node)

1. load Spark module (see [environment modules](guides/cluster_guide.md)), which will put environment variables that Spark uses in your path
2. Allocate resources, using Slurm, that Spark cluster will use for its cluster. This means allocate cpu cores, memory, and maxtime for your job. Number of nodes will be 1 if you use  
3. Run `spark-submit` or `spark-shell` with --master local[4] (don't use `local[*]`, as it might confuse the Spark Standalone if it derives its information from node spec rather than what you have been allocated)

In code:   
```
srun -N 1 -n 1 -c 4 -t 01:00:00 --mem=10Gb --pty bash   # to ask for single node and 4 cores for 1 hour. This has the effect of placing you on a compute node
module use /projects/community/modulefiles
module load spark/2.3.0-kp807    # finds spark installed in community software
spark-shell --master local[4]    # starts spark-shell. Use the same number of cores (4 in this case) that you asked in srun
```

To look at Spark GUI, see how to do [tunnelling](howtos/jupyter.md). 

## Advanced - spin up multinode cluster

TODO

## Submit a batch script

Check out the `slurm_examples` directory in our [github repo](https://github.com/rutgers-oarc/training) for an example. 

## Start Spark notebook 

```
srun -N 1 -n 1 -c 4 -t 01:00:00 --mem=10Gb --pty bash   # to ask for single node and 4 cores for 1 hour. This has the effect of placing you on a compute node
module use /projects/community/modulefiles
module load spark-notebook/2.2.0-kp807
spark-notebook                # starts spark notebook on port 9000, now you have to tunnel to reach it
```

See [spark notebook demo](https://www.youtube.com/watch?v=tms_KpmTTkQ) for more information and options.
To tunnel, in another terminal establish tunnelling with this command (supposing your compute node is slepner022, and your netid is kp807:   
```ssh -L 9000:slepner022:9000 kp807@amarel.hpc.rutgers.edu```

