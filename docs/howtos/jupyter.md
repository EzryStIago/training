# Tunneling

This is a technique for connecting from your local machine to a remote web server, such as Jupyter notebook running on the compute node of a cluster. 

1. run jupyter notebook as a slurm job 
2. find out on which compute node jupyter notebook ended up
3. in another terminal establish "port forwarding" - from local port 9999 to port 8889 on the specific node (slepner009 here)
```
# This procedure works not only with jupyter notebook, but any web app running on a compute node

#command to run jupyter notebook 
srun -p main -c 1 -t 10:00 --error=slurm.%N_%j.err --output=slurm.%N_%j.out jupyter notebook --no-browser --ip=0.0.0.0 --port=8889

# to start jupyter succintly 
sbatch start_jupyter.sh      #start jupyter job on port 8889
squeue -u kp807              #find out on which node is the job running - say it's slepner009

# in another terminal establish "port forwarding" - from local port 9999 to port 8889 on the specific node (slepner009 here)
ssh -L 9999:slepner009:8889 kp807@amarel.hpc.rutgers.edu   # modify slepner009, the ports, the netID
```

Video expaining the steps above: 
[![Port Tunneling](http://img.youtube.com/vi/Nb2ERdr5KD0/0.jpg)](http://www.youtube.com/watch?v=Nb2ERdr5KD0)


# How to launch Jupyter notebook on the cluster

There is an already available installation of Anaconda 5.1.0 with Python 3.6.4 on the cluster, which contains both Jupyter notebook and Jupyter lab. To load the necessary module execute these commands: 
```
module use /projects/community/modulefiles
module load py-data-science-stack
```

Copy this into a script file like `start_jupyter.sh`

```
#!/bin/bash

#SBATCH --partition=main             # Partition (job queue)
#SBATCH --job-name=jupyter          # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=4000                 # Real memory (RAM) required (MB)
#SBATCH --time=01:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file
#SBATCH --error=slurm.%N.%j.err      # STDERR output file (optional) 

export XDG_RUNTIME_DIR=$HOME/tmp   ## needed for jupyter writting temporary files

module use /projects/community/modulefiles 
module load py-data-science-stack         # loads anaconda
source activate pytorch-0.4.0

#run system-installed jupyter notebook on port 8889. ip 0.0.0.0 means "any ip as interface"
srun jupyter notebook --no-browser --ip=0.0.0.0 --port=8889
```

Then run `sbatch start_jupyter.sh` on amarel login node. Now you need to find which node the jupyter notebook is running at. Do `squeue -u <your net id>` to see the slurm jobs you are running. 
You should then do the port tunneling described in the previous section and open local browser at that port. 

Youtube video that explains this: 
[![Graphical programs via FastX](http://img.youtube.com/vi/8UQHLO85rlQ/0.jpg)](https://www.youtube.com/watch?v=8UQHLO85rlQ)



