Connecting to the remote Linux cluster makes running graphical programs remotely tricky, because the graphical program runs on the remote computer, yet, it must be displayed on the local machine such as laptop. In general, command line interaction with the cluster is normally preferred and even more efficient than a GUI. However, there are times when running a graphical program cannot be avoided, for example, running a debugger for a code running on the cluster. 

There is a convenient way to run a graphical program remotely, for example, using FastX. Here is the procedure: 

1. Go to htts://amarel.hpc.rutgers.edu:3443  and log in (you must be either on campus, or connected through VPN)
2. Click on Launch session
3. Click on xterm
4. Run this command that will ask the resource scheduler to put you on a compute node, rather than a login node (where you shouldn't be running intensive computations)  
```srun -p main -N 1 -c 2 --mem=4Gb -t 1:00:00 --pty /bin/bash```  
This command puts you on main partition, asks for 2 cores on 1 node, asks for 4 Gb or memory and time of 1 hour and runs interactive shell. After executing this, you should notice that the name of the node changed from `amarel` to `slepner036` or some such. You can request whatever resources you deem necessary for your work, but keep in mind that bigger requests are typically placed further down in the queue. 
5. Start the program with graphical interface from the terminal window.

Here is the video walking through these steps: 

[![Graphical programs via FastX](http://img.youtube.com/vi/2LaGZtOv8Lc/0.jpg)](https://www.youtube.com/watch?v=8UQHLO85rlQ)


