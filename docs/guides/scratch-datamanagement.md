# Guidelines for managing data on scratch

### Table: Data management on OARC cluster


<table>
  <tr>
    <td> Fileset </td>
    <td> Users </td>
    <td> quota </td>
    <td> Time </td>
    <td> Back up </td>
    <td> External storage options </td>
    <td> Transfer Tools </td>
  </tr>
  <tr>
    <td> scratch
    <td> general
    <td> No  Quota
    <td> 90 days
    <td> no
   <td rowspan="2"> - Personal devices: Laptop, Desktop, USB drives, etc. <br>
         - Online drives: Box(unlimited for Rutgers), Google Drive, One Drive, Dropbox, etc. <br>
         -Public cloud storages: GCP, Azure, AWS, etc. <br>
       </td>
   <td rowspan="2">  -CLI tools: scp, rsync, rcloud <br>
         -GUI tools: Filezilla, WinSCP, Open OnDemand, rclone-browser (via sirius3) <br>
       </td>
  </tr>
  <tr>
    <td> home </td>
    <td> general </td>
    <td> 100GB </td>
    <td> Valid account </td>
    <td> yes </td>
    
  </tr>

  <tr>
    <td> project </td>
    <td> node owners </td>
    <td> varies </td>
    <td> varies </td>
    <td> yes </td>
   
  </tr>
  
</table>


Scratch offers unlimited temporary storage. You can processing large volumes of data in scratch, but the data cannot be older than 90 days. The files older than 90 days are automatically removed, and there is no backup. If you exceed the disk quota, your jobs would fail with a message that you cannot write anymore. So, it is a best practice to move the data from scratch as soon as possible.

## Internal and external storage options

For long term storage, you can utilize the internal storage at OARC's cluster or any of the external storage. OARC's internal storage includes home and project file sets.  If your data can be accommodated in your home directory, it is okay to keep the final results in the home. Node owners may want to use their project space. If you don't have enough space in the OARC cluster, you need to move the data to an external storage.

Here we point out a few options to consider for external storage. Personal devices such as laptop, desktop, workstations, or USB drives can serve as external storage.  Third part vendors offer several storage options that are either free or paid services. Box offers unlimited free storage for Rutgers community. Google offers unlimited  free storage for the academic community  by a service called Google Drive File Stream. One Drive and DropBox services are free to use up to a set quota.  Public cloud providers such as AWS, Azure, and GCP charge a fee to move the data in and out of the cloud.

## Data transfer tools

There are several tools and methods available to move the data that are based on command line interface (CLI) or Graphical User Interface (GUI). You can use command line tools such as scp, rsync, or rcloud to move the data from OARC cluster to any external storage. scp and rsync are suitable choices if you are moving data from OARC to your laptop or desktop that runs linux or MacOS. rcloud is an excellent choice to move your data to many storage services such as Box, DropBox, cloud, etc. GUI tools such as Filezilla, WinSCP, rclone-browser, and Open OnDemand are popular choices to move the data from variety of storage end points. Among GUI tools, FIlezilla, rclone-browser, and Open OnDemand tools should work on any operating systems (Linux, MacOS, Windows) since they interact through a web-browser.

Please refer our user guide for more details about how to use specific CLI tools such as scp, rsync, and rclone. Here is a simple command line example to move your data using scp. 

    scp your-netid@amarel.hpc.rutgers.edu:<file location on oarc cluster> <file location on your device> 

where ‘your-netid’ is your Rutger’s netid, which is same as your account name on OARC cluster. For example, to move data from scratch space located at /scratch/your-netid/your-data-file to a home directory on your laptop using scp command, type the following in your terminal shell prompt on your laptop.

    scp  your-netid@amarel.hpc.rutgers.edu:/scratch/your-net-id/your-data-file ~/your-data-file




