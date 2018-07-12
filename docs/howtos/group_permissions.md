
Suppose I have a directory that I want to share with a collaborator, `/home/kp807/myproject`, which I want to make readable for the collaborator and a subdirectory `/home/kp807/myproject/code` where I will also give the collaborator the permission to write. I do not want the collaborator to see everything in my home directory `/home/kp807`. 

## Using granular FACL controls

FACL (File Access Control List) has two important commands, `setfacl` and `getfacl`. `getfacl` works just like `ls`. Example: 

```
[kp807@amarel1 ~]$ getfacl /home/kp807
getfacl: Removing leading '/' from absolute path names
# file: home/kp807
# owner: kp807
# group: kp807
user::rwx
user:dm1084:--x
group::---
other::---
```

So to setup the correct permissions, we'd need to do: 

```
setfacl -m u:kholodvl:x /home/kp807/  #modify permissions (-m) for user kholodvl and make /home/kp807/ executable for his commands
setfacl -m u:kholodvl:rx /home/kp807/myproject
setfacl -m u:kholodvl:rwx /home/kp807/myproject  
```

**Important Note**: Check defaults with which your files and directories are created. 

## Using group membership

Use commands `chmod` (or `chown` as appropriate) to modify permissions and ownership, respectively, of directories or files: 

```
chmod g+x /home/kp807/                  #make the directory executable so other person's command can tranverse the file tree to myproject
chmod g+rx /home/kp807/myproject/       #make the directory readable and executable for the group I belong to
chmod g+rwx /home/kp807/myproject/code  #make the directory readable, writable and executable for the group I belong to
```

[kp807@amarel1 ~]$ mkdir myproject
[kp807@amarel1 ~]$ mkdir myproject/code
[kp807@amarel1 ~]$ setfacl -m u:kholodvl:x /home/kp807/  #modify permissions (-m) to user kholodvl and make /home/kp807/ executable for his commands
[kp807@amarel1 ~]$ setfacl -m u:kholodvl:rx /home/kp807/myproject
[kp807@amarel1 ~]$ setfacl -m u:kholodvl:rwx /home/kp807/myproject
[kp807@amarel1 ~]$ touch /home/kp807/myproject/readme.txt
[kp807@amarel1 ~]$ touch /home/kp807/myproject/code/writeme.txt
[kp807@amarel1 ~]$ ls /home/kp807/myproject/code/

rm -rf /home/kp807/myproject
mkdir /home/kp807/myproject
mkdir /home/kp807/myproject/code
ls -la /home/kp807/myproject
getfacl /home/kp807/myproject
ls -la /home/kp807/myproject/code
getfacl /home/kp807/myproject/code

