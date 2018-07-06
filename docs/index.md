Welcome to the Office of Advanced Research Computing at Rutgers! 

> OARC is a university-wide initiative that aims to develop and implement a strategic vision for centralizing the advanced research computing and data cyberinfrastructure (ACI) ecosystem at Rutgers. OARC has the goal of providing Rutgers researchers with essential computing and data handling capabilities, and students with necessary exposure and training, through centralized resources, services and training.

For more information on OARC, including how to [get access](https://oarc.rutgers.edu/access) or [become owners](https://oarc.rutgers.edu/services/condo-model/), please visit [our web page](https://oarc.rutgers.edu/)

These pages are a collection of resources to help you to utilize the cluster more effectively. Even if you are a very experienced Linux user, you will want to read [cluster user guide](guides/Cluster_User_Guide/) as it has slurm tips and examples. 

## WARNING - READ!

- Do not run large computational jobs on the **login** node. Use slurm to allocate resources on the compute node. Failure to respect the golden rule can get your account suspended.
- Number of jobs to submit should not exceed 5,000 at a time - queue is currently 10,000 and slurm will go to sleep if that is exceeded.
- Right-size your jobs: jobs may run sooner if the resources requested are smaller and your future priority is degraded if you have asked for (even if not used) a lot of resources lately.
- use `-n` option in sbatch only if you know what you are doing. Multithreaded but not multi-node jobs should use `-c` option. 
- TODO: add more complete list

## Learning paths

- If you understand material in [lmod cheatsheet](CheatSheet/#environment-modules-lmod-cheatsheet) and [slurm cheatsheet](CheatSheet/#slurm-cheatsheet) you are good to go.
- For users familiar with Linux but new to `slurm`, check out [environment modules](guides/Cluster_User_Guide/#loading-software-modules) and [intro to slurm](guides/Cluster_User_Guide/#running-slurm-jobs). 
- For users not familiar with Linux, please familiarize yourself with Linux through tutorials and cheatsheets, or use [OnDemand](guides/Cluster_User_Guide/#ondemand-gui-for-the-cluster) until you feel more comfortable with Linux.


