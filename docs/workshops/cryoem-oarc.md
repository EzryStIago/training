# Overview

## As a part of ongoing collaboration between Rutgers [**Office of Advanced Research Computing (OARC)**](http://oarc.rutgers.edu) and Rutgers [**Institute for Quantitative Biomedicine**](http://iqb.rutgers.edu/) we are proud to announce that OARC clusters were chosen to be a primary computational platform to support the 6th Annual Interdisciplinary Quantitative Biology Winter Boot Camp for  Single-Particle Cryo-Electron Microscopy January 7 - 11, 2019 https://iqb.rutgers.edu/bootcampwinter2019

### We are proud of being at the frontier of science and helping to spread the knowledge of Rutgers scientists with the  national and international scientific community.

# Setup

## Connect to the OARC cluster 
The boot camp participants are expected to bring their own laptops with a **Chrome web browser** installed.
The preferred method to connect to the cluster is through remote desktop (fastX) in a web browser.
Each participant will be granted an access to an individual compute node for the duration of the boot camp. 
The compute node will be accessible only on the Rutgers network or if you are connected via [VPN](https://ssl-vpn.rutgers.edu/dana-na/auth/url_default/welcome.cgi). The access port is unique for each participant and should not be shared.


 **Remote desktop via FastX**: 
 In your browser, go to https://cryoem-oarc.hpc.rutgers.edu:port <p>
 Enter your Rutgers netID and Password and click on the **Log in** button.
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/fastX_login.png" width="500" height="500" center> <p>
 
 In the new window, click on **Launch Session** on the left, select **XFCE desktop** and then click on **Launch** button on the right.<p>
<img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_sessions.png"> 
  
<img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_XFCE.png"> 
<p> 
 Ignore the warning<p> 
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_popup.png" width="500" height="155" center> 
<p>
 It will open a new desktop window on a compute node<br> 
 At the bottom of the page find a terminal icon, a black rectangular that looks like an old fashion TV box. Click on it. <br>
 In the terminal window type <i>scipion</i> to start Scipion <p>
 
 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/FX_desktop.png" width="600" height="349"> <p> 
 
 ### Please be patient, it may take up to 30sec to start Scipion for the first time <p>

 <img src="https://github.com/rutgers-oarc/training/blob/master/workshops/scipion_projects.png" width="500" height="212"> <p>
  
 
 **Reconnecting to the FastX remote desktop session**:
 - Open the web browser and go to https://cryoem-oarc.hpc.rutgers.edu:port <br>
 - Login into your account with your Rutgers netID and a password.<br>
 - Click on the **XFCE Desktop** session in the My Session list. <br>
 





## CRYO-EM Software

This is a list of software installed on OARC clusters for the workshop

|Software|Description|URL|
|------|-----------------------------|-----------------------|
|Scipion|An image processing framework for obtaining 3D models of macromolecular complexes using Electron Microscopy (3DEM)|[[link]](https://github.com/I2PC/scipion/wiki)|
|Scipion packages|Scipion integrated packages: Bsoft, cryoFF, CTFFIND, EMAN, ETHAN, Frealign, Gautomatch, gCTF, gEMpicker, IMAGIC-4D, Localrec, Motioncorr/dosefgpu, Motioncorr2, Relion, ResMap, SIMPLE, SPIDER, Xmipp, Unblur & summovie|[[link]](https://github.com/I2counting)|
|Chimera|A program for interactive visualization and analysis of molecular structures and related data, including density maps|[[link]](https://www.cgl.ucsf.edu/chimera/)| <p>
 
You can also check what Scipion packages are installed from a command line.<br>
Open a new terminal window and type
```
scipion install --help
```



### Understanding your data space

You have two main spaces on the OARC cluster. These are: 

- your home directory (100Gb) - `/home/netid/` 
- your scratch directory (unlimited, no-backup)- `/scratch/netid/` 

The difference between two is that scratch is not backed up and the read/write speed is higher for the scratch. <br> 
For the purpose of our workshop all the data and computational output will be held in `/scratch/netid/ScipionUserData`.<br>
**In the linux environment the names of files/directories are case sensative. Make sure that you properly type upper and lower case letters in the path when needed.**  <p>

When you start Scipion for the first time you may see some warning messages about security of python cashe directory. 
```
>>>>> python  /gpfs/home/projects/community/scipion/openmpi3/pyworkflow/apps/pw_manager.py 
/gpfs/home/projects/community/scipion/openmpi3/software/lib/python2.7/site-packages/setuptools-5.4.1-py2.7.egg/pkg_resources.py:1049: UserWarning: /home/db1102/.python-eggs is writable by group/others and vulnerable to attack when used with get_resource_filename. Consider a more secure location (set with .set_extraction_path or the PYTHON_EGG_CACHE environment variable).
```
You may ignore this message or if you want to get rid off its annoyance, in a new terminal window type the following command:
```
chmod g-w ~/.python-eggs
```





