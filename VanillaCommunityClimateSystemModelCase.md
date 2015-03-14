# Introduction #

## Creating a case ##

First attempt:
```
-bash-2.05b$ ./create_newcase -case vanilla -mach jazz -rez T31_gx3v5 -compset B
unknown input, invoke create_newcase with no arguments for usage
```
Changed _rez_ to _res_

Second attempt:
```
-bash-2.05b$ ./create_newcase -case vanilla -mach jazz -compset B -res T31_gx3v5
Successfully added env_mach.jazz to /home/steder/ccsm3_0/scripts/vanilla
Successfully created new case root directory /home/steder/ccsm3_0/scripts/vanilla
```

Okay, time to configure and build:
```
-bash-2.05b$ ls
ccsm_utils/      create_test*  dead/       deadtest3/  pyccsm_dead/  vanilla/
create_newcase*  data/         deadtest2/  pyccsm/     README*       xdeadtest/
-bash-2.05b$ 
-bash-2.05b$ cd vanilla
-bash-2.05b$ ls
configure*  env_conf  env_mach.jazz*  env.readme*  env_run  SourceMods/
-bash-2.05b$ 
```

## Configure ##

```
-bash-2.05b$ ls
configure*  env_conf  env_mach.jazz*  env.readme*  env_run  SourceMods/
-bash-2.05b$ ./configure -mach jazz
-------------------------------------------------------------------------
 Generating resolved setup namelist-prestage and build scripts
 and installing build scripts for libraries  
    See directory /home/steder/ccsm3_0/scripts/vanilla/Buildexe/ 
    See directory /home/steder/ccsm3_0/scripts/vanilla/Buildnml_Prestage/ 
    See directory /home/steder/ccsm3_0/scripts/vanilla/Buildlib/ 
-------------------------------------------------------------------------
-------------------------------------------------------------------------
Generating resolved build script, batch run script 
   See /home/steder/ccsm3_0/scripts/vanilla/vanilla.jazz.build
   See /home/steder/ccsm3_0/scripts/vanilla/vanilla.jazz.run
-------------------------------------------------------------------------
no long term archiving script (vanilla.jazz.l_archive) will be created for this machine
*************************************************************************
    Successfully configured the model for machine jazz
*************************************************************************
```

## Build ##

Command:
```
-bash-2.05b$ ./vanilla.jazz.build 
```

Output:
```
-------------------------------------------------------------------------
 Preparing T31_gx3v5 component models for execution 
-------------------------------------------------------------------------
 - Create execution directories for atm,cpl,lnd,ice,ocn
 - If a restart run then copy restart files into executable directory 
ccsm_getrestart: get /home/steder/exe/vanilla restarts from /home/steder/archive/vanilla/restart
 - Check validity of configuration
 - Determine if build must happen (env variable BLDTYPE)
 - Build flag (BLDTYPE) is TRUE
 - Build Libraries: esmf, mph, mct
Mon Jan 29 15:58:40 CST 2007 esmf.buildlib.070129-155838
Mon Jan 29 15:59:15 CST 2007 mph.buildlib.070129-155838
Mon Jan 29 15:59:18 CST 2007 mct.buildlib.070129-155838
 - Create model directories for each platform
 - Determine if models must be rebuilt
 - Build model executables, create namelist files, prestage input data
Mon Jan 29 15:59:28 CST 2007 /home/steder/exe/vanilla/cpl/cpl.log.070129-155838
Mon Jan 29 15:59:28 CST 2007 /home/steder/exe/vanilla/cpl/cpl.buildexe.070129-155838
Mon Jan 29 15:59:42 CST 2007 /home/steder/exe/vanilla/ice/ice.log.070129-155838
Mon Jan 29 15:59:42 CST 2007 /home/steder/exe/vanilla/ice/ice.buildexe.070129-155838
Mon Jan 29 16:00:00 CST 2007 /home/steder/exe/vanilla/lnd/lnd.log.070129-155838
Mon Jan 29 16:00:01 CST 2007 /home/steder/exe/vanilla/lnd/lnd.buildexe.070129-155838
Mon Jan 29 16:00:22 CST 2007 /home/steder/exe/vanilla/ocn/ocn.log.070129-155838
Mon Jan 29 16:00:22 CST 2007 /home/steder/exe/vanilla/ocn/ocn.buildexe.070129-155838
Mon Jan 29 16:06:01 CST 2007 /home/steder/exe/vanilla/atm/atm.log.070129-155838
Mon Jan 29 16:06:01 CST 2007 /home/steder/exe/vanilla/atm/atm.buildexe.070129-155838
 - Create MPH input file and link into all model dirs
 - Create shr_msg_stdio chdir/stdin/stdout data file
-------------------------------------------------------------------------
 - CCSM BUILD HAS FINISHED SUCCESSFULLY                                  
-------------------------------------------------------------------------
```