

https://software.rutgers.edu/product/3085  - Download tab

https://oit.rutgers.edu/help - Help Desk and Local IT Support | Office of Information Technology - Rutgers University


Have you connected the vpn through the instruction for Mac: https://sasit.rutgers.edu/documentation-help/accessing-sas-servers/304-configuring-the-vpn-client-for-mac? 


2.       If you connected the vpn, could you see this website: https://ask.oarc.rutgers.edu/question/157/i-cant-connect-to-ondemand/ ?



do you see that in the URL, it starts with https://ssl-vpn.rutgers.edu? That means, that VPN used is the web-based VPN, and not the client VPN. 


I get the same error as failed to connect to server

please ask OIT for help https://oit.rutgers.edu/help

# Didact for classes

Didact is a small cluster dedicated to teaching. It can be accessed at [this didact URL](https://didact.oarc.rutgers.edu/). It is equipped with R, Jupyter notebooks, and a Linux desktop which can be used to start other graphical applications like stata and SPSS. 

# Info for instructors

- we have an instructor packet for you - please contact us at help@oarc.rutgers.edu
- Please [fill in access request for Amarel](https://oarc.rutgers.edu/access) before the start of class and say you want to teach a class
- we will create accounts in one go for your students and notify you. We need a list of student netids and a list of student email addresses in a pasteable form (csv or xlsx, not pdf)
- you need to be able to troubleshoot students connection problems to some extent - please read VERY carefully VPN access problems section - these is the most common pitfall students face
- you should communicate to students when their accounts are ready, expectations on OARC support, and information about expected account removal


# Info for students

- OARC can't provide help for issues outside the control of our office - most frequently VPN access. For this, you need to contact [OIT helpdesk](https://oit.rutgers.edu/help). 
- You should submit access requests through [course registration form](https://oarc.rutgers.edu/course-registration-form/) even if your account already exists. This is for our accounting purposes. 
- Every didact account requires also an Amarel account, which will be automatically provided to you. You will be able to continue using the software for your research after the course ends.

# Steps to run stata/spss or other programs

- Go the the form to launch an interactive desktop [on didact](https://didact.oarc.rutgers.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop/didact/session_contexts/new)  
or [on Amarel](https://ondemand.hpc.rutgers.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop_local/session_contexts/new) and click blue button "Launch"
- A new tab will open, where the session is being created. After about 30 seconds when the session is ready, a blue button will appear "Launch noVNC in New Tab" - click on it, and your desktop is ready to be used
- icons should appear on your desktop. If you see no icons, please see [this question](https://ask.oarc.rutgers.edu/question/55/how-can-i-get-icons-on-amarel-desktop/) - you need to copy the files from /projects/community/desktops to Desktop folder in your home folder. You can do this eitherby using a terminal linux command, or by using the file explorer. 

# User forum

Please find answers to some common questions and help each other by asking and answering questions at [our user forum](https://ask.oarc.rutgers.edu/questions/)

# VPN access

To access any OARC services from off-campus, you need a VPN. This is for Rutgers security reasons, and is a campus-wide service administered outside our office. Please direct any questions to [OIT helpdesk](https://oit.rutgers.edu/help). For your convenience, here is a brief collection of links and directions: 

- you can download VPN client at https://software.rutgers.edu/product/3085, from the "Download" tab. 
- Please DO NOT USE [web-based VPN access](https://ssl-vpn.rutgers.edu/), as it interferes with the OpenOnDemand application (see Known Issues below)
- instructions to install VPN client on a Mac from [School of Arts and Sciences](https://sasit.rutgers.edu/documentation-help/accessing-sas-servers/304-configuring-the-vpn-client-for-mac)
- if you have a newer Windows laptop, you need 64bit version, and not 32bit version

# Known issues

- `this page does not exist` - this means you connected from off campus and you are not connected through VPN. If you can't see [this page](https://ask.oarc.rutgers.edu), you are not on campus and not connected through VPN. Please download the VPN client. 
- `failed to connect to VNC server` - one of two things is possible: 
  * you connected through web-based VPN, which doesn't work well with ondemand. You will recognize this problem by looking at the URL and it starts with `https://ssl-vpn.rutgers.edu/pun/sys/dashboard`and your screenshot [looks like this](cant_connect_to_vnc_server.png)
  * you are using Safari as your browser. Please try Firefox instead. 

# Debugging questions

If you have trouble connecting to didact, please try these steps to isolate the problem: 

- please read about VPN access - this trips users most often. Either you didn't connect through VPN, or you didn't download the VPN client, and instead used [this webpage](https://ssl-vpn.rutgers.edu/) to login to VPN. 
- try a different browser (e.g. Safari doesn't work for some people)
- try a different computer
- try a different operating system
- ask your peers if they are successful in connecting
- raise the problem with your instructor


