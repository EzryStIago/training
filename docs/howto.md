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


