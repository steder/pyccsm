| [PyCCSM](PyCCSM.md) |
|:--------------------|

# Installation #

## PyMCT ##
[InstallPyMCT](InstallPyMCT.md)

## CCSM ##

Assuming you have a CCSM built on your platform of choice, it's relatively straightforward to modify it to be compatible with our Python coupler.

The main difference between PyCCSM and CCSM is that MPH is not used for initialization.  Other then that the 2 codebases are identical.

Replacing MPH without own initialization routine is a 3 step process.  First, we add our initialization routine to the CCSM shared source.  Next, we copy the main program of each component model into our CCSM Case's _SourceMods_ directory and slightly modify the source to call the new initialization routine.  Finally, we recompile CCSM using the standard CCSM build system.

### Step 1: Adding our alternate start up ###

**Download** the PyCCSM initialization code from this site:
[cpl\_comm\_mod.F90.withoutMPH](http://pyccsm.googlecode.com/svn/trunk/cpl_comm_mod.F90.withoutMPH)

**Locate** your CCSM source.  It is likely in a subdirectory like _ccsm3\_0_on your system.
Find it and change directories to it.
```
$ pwd
/home/steder/
$ ls | grep -i ccsm
ccsm3_0/
CCSM.txt
$ cd ccsm3_0/
$ ls
Copyright*  models/  README*  scripts/
$ cd models/
$ ls
atm/  bld/  cpl/  csm_share/  dead/  ice/  lnd/  ocn/  utils/
$ cd csm_share/
$ ls
cpl/  shr/
$ cd cpl/
$ ls
cpl_bundle_mod.F90*    cpl_control_mod.F90*  cpl_interface_mod.F90*  cpl_map_mod.F90*
cpl_comm_mod.F90*      cpl_domain_mod.F90*   cpl_iobin_mod.F90*      cpl_mct_mod.F90*
cpl_const_mod.F90*     cpl_fields_mod.F90*   cpl_iocdf_mod.F90*      netcdf.inc*
cpl_contract_mod.F90*  cpl_infobuf_mod.F90*  cpl_kind_mod.F90*
```_

Now we need to take this file and put it in place of the original.  So first we'll move the original file someplace so we can retrieve it if we ever need it.

```
$ cp cpl_comm_mod.F90 cpl_comm_mod.F90.backup
```

Then we'll go ahead and download the new version to this directory
```
$ wget http://pyccsm.googlecode.com/svn/trunk/cpl_comm_mod.F90.withoutMPH 
--17:44:39--  http://pyccsm.googlecode.com/svn/trunk/cpl_comm_mod.F90.withoutMPH
           => `cpl_comm_mod.F90.withoutMPH'
Resolving pyccsm.googlecode.com... done.
Connecting to pyccsm.googlecode.com[66.102.1.82]:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 17,137 [text/plain]

100%[============================================================>] 17,137        74.05K/s    ETA 00:00

17:44:40 (74.05 KB/s) - `cpl_comm_mod.F90.withoutMPH' saved [17137/17137]

$ ls
cpl_bundle_mod.F90*          cpl_contract_mod.F90*  cpl_interface_mod.F90*  cpl_mct_mod.F90*
cpl_comm_mod.F90*            cpl_control_mod.F90*   cpl_iobin_mod.F90*      netcdf.inc*
cpl_comm_mod.F90.backup*     cpl_domain_mod.F90*    cpl_iocdf_mod.F90*
cpl_comm_mod.F90.withoutMPH  cpl_fields_mod.F90*    cpl_kind_mod.F90*
cpl_const_mod.F90*           cpl_infobuf_mod.F90*   cpl_map_mod.F90*
```

And finally, we'll copy the new initialization code over the original.  Remember, we should also have a backup of the original initialization code somewhere.
```
$ cp cpl_comm_mod.F90.withoutMPH cpl_comm_mod.F90
cp: overwrite `cpl_comm_mod.F90'? y
```

### Step 2: Setting up our SourceMods ###

Now we need to modify the main of each component program to use our replacement MPH:

We'll most likely want to create a brand new CCSM case for our PyCCSM runs.  So we'll start by locating our CCSM scripts.  They're usually located in a subdirectory with the CCSM source distribution.

Assuming our CCSM 3 source is downloaded to our home directory in a directory called _ccsm3\**0** our scripts will be:
```
steder@joxaren:~/ccsm3_0/scripts$ pwd
/home/steder/ccsm3_0/scripts
steder@joxaren:~/ccsm3_0/scripts$ cd  
steder@joxaren:~$ cd ccsm3_0/
steder@joxaren:~/ccsm3_0$ pwd
/home/steder/ccsm3_0
steder@joxaren:~/ccsm3_0$ ls
Copyright  README  models/  scripts/
steder@joxaren:~/ccsm3_0$ cd scripts/
steder@joxaren:~/ccsm3_0/scripts$ ls
CCSM.txt  ccsm_utils/      create_test*  deadtest/  jazz_cases/
README    create_newcase*  data/         fulltest/  pyccsmtest/
```_

### Dead Model Case ###

First you will need to create a dead model case in CCSM.  You can do this by running the CCSM script _create\**newcase** with the following flags:_

```
steder@joxaren:~/ccsm3_0/scripts$ ./create_newcase -case pyccsm-dead-models -mach jazz -res T31_gx3v5 -compset X 
Successfully added env_mach.jazz to /home/steder/ccsm3_0/scripts/pyccsm-dead-models
Successfully created new case root directory /home/steder/ccsm3_0/scripts/pyccsm-dead-models
```

Now let's change into our newly created test case directory:

```
steder@joxaren:~/ccsm3_0/scripts$ cd pyccsm-dead-models/
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ ls
SourceMods/  configure*  env.readme  env_conf  env_mach.jazz*  env_run
```

Assuming that the CCSM model source code is in the default location:

```
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ ls ../../models/                    
atm/  bld/  cpl/  csm_share/  dead/  ice/  lnd/  ocn/  utils/
```

We'll now copy the dead model main program into our case's directory:

```
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ cp ../../models/dead/dead6/dead.F90 .
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ ls
SourceMods/  configure*  dead.F90  env.readme  env_conf  env_mach.jazz*  env_run
```

Now we'll edit dead.F90 to use our new initialization routine.  We're going to make just a few small modifications to this file.  Here's the first.

**Add a use statement to include our new code**
```
subroutine dead_init()

! !USES:                                                                                        

   use cpl_interface_mod
   use cpl_fields_mod
   use cpl_contract_mod
   use cpl_control_mod
   use shr_const_mod
   use shr_timer_mod
   use shr_sys_mod
   use shr_msg_mod
   use shr_kind_mod
   use data_mod

   implicit none
```

The above needs to be modified by adding the following use statement (these are lines 59-73):
```
subroutine dead_init()

! !USES:                                                                                        

   use cpl_interface_mod
   use cpl_fields_mod
   use cpl_contract_mod
   use cpl_control_mod
   use shr_const_mod
   use shr_timer_mod
   use shr_sys_mod
   use shr_msg_mod
   use shr_kind_mod
   use data_mod
   use cpl_comm_mod 

   implicit none
```



Next, look for a little later on in the file(lines 204-206):
```
   !----------------------------------------------------------------------------   
   ! initialize local MPH & MPI                                                                 
   !----------------------------------------------------------------------------                
   call cpl_interface_init(myModelName,local_comm)
```

We'll comment out the old call to cpl\_interface\_init and replace it with the following:
```
   !----------------------------------------------------------------------------   
   ! initialize local MPH & MPI                                                                 
   !----------------------------------------------------------------------------                
   !call cpl_interface_init(myModelName,local_comm)
   call cpl_comm_init(myModelName, local_comm)
```

Now we have a modified main for all 4 component dead models.  Let's place this in _SourceMods_ so our modifications will be included the next time we build this case.

```
$ ln -s ../../dead.F90 SourceMods/src.xice/dead.F90
$ ln -s ../../dead.F90 SourceMods/src.xatm/dead.F90
$ ln -s ../../dead.F90 SourceMods/src.xlnd/dead.F90
$ ln -s ../../dead.F90 SourceMods/src.xocn/dead.F90
$ ls
SourceMods/  configure*  dead.F90  env.readme  env_conf  env_mach.jazz*  env_run
$ ls -l SourceMods/src.x*/
SourceMods/src.xatm/:
total 8
-rw-r--r--  1 steder users 59 Feb  8 23:10 README
lrwxrwxrwx  1 steder users 14 Feb  8 23:30 dead.F90 -> ../../dead.F90

SourceMods/src.xice/:
total 8
-rw-r--r--  1 steder users 59 Feb  8 23:10 README
lrwxrwxrwx  1 steder users 14 Feb  8 23:30 dead.F90 -> ../../dead.F90

SourceMods/src.xlnd/:
total 8
-rw-r--r--  1 steder users 59 Feb  8 23:10 README
lrwxrwxrwx  1 steder users 14 Feb  8 23:30 dead.F90 -> ../../dead.F90

SourceMods/src.xocn/:
total 8
-rw-r--r--  1 steder users 59 Feb  8 23:10 README
lrwxrwxrwx  1 steder users 14 Feb  8 23:30 dead.F90 -> ../../dead.F90
```

### Live Model Case ###

The live model case is the most complicated because each live model has a different initialization code.  However, it really isn't much more complicated then the dead model case.  I'll simply go over the important files and their modifications briefly.

```
steder@joxaren:~/ccsm3_0/scripts/pyccsmtest$ ls SourceMods/src.cam SourceMods/src.clm/ SourceMods/src.csim/ SourceMods/src.clm/
SourceMods/src.cam:
README  cam.F90  

SourceMods/src.clm/:
README  clm_csmMod.F90

SourceMods/src.csim/:
README  ice_coupling.F

SourceMods/src.pop/:
README  forcing_coupled.F
```

**Editing cam.F90**
First, get the source file:
```
cp ../../models/atm/cam/src/control/cam.F90 .
ln -s ../../cam.F90 SourceMods/src.cam/cam.F90
```

Now update the use statements (starts at line 34):
```
#if (defined COUP_CSM)
   use shr_msg_mod
   use cpl_interface_mod
   use cpl_fields_mod
   ! Add this use statement:
   use cpl_comm_mod 
#endif
```

Now search for "cpl\_interface\_init" and replace it with "cpl\_comm\_init" (starts at line 100)
```
!                                                                                               
! Initialize internal/external MPI if appropriate                                               
!                                                                                               
#if ( defined COUP_CSM )
   call shr_msg_stdio('atm')
   !call cpl_interface_init(cpl_fields_atmname,mpicom)
   call cpl_comm_init( cpl_fields_atmname, mpicomm )
#endif
```

**Editing clm\_csmMod.F90**
First, get the source file:
```
$ cp ../../models/lnd/clm2/src/main/clm_csmMod.F90 .
$ ln -s ../../clm_csmMod.F90 SourceMods/src.clm/clm_csmMod.F90
```

Now update the use statements (starts at line 23):
```
! !USES:                                                                                        
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clm_varpar
#if (defined SPMD)
  use spmdMod        , only : masterproc, mpicom
  use spmdGathScatMod, only : gather_data_to_master
#else
  use spmdMod        , only : masterproc
#endif
  use mpiinc
  use cpl_fields_mod
  use cpl_contract_mod
  use cpl_interface_mod
  ! add this use statement!
  use cpl_comm_mod
```

Now search for "cpl\_interface\_init" and replace it with "cpl\_comm\_init" (Starting at line 206)
```
! !LOCAL VARIABLES:                                                                             
!------------------------------------------------------------------------                       

   !call cpl_interface_init(cpl_fields_lndname,mpicom)
   call cpl_comm_init(cpl_fields_lndname,mpicomm)
```

**Editing ice\_coupling.F**

First, get the source file:
```
$ cp ../../models/ice/csim4/src/source/ice_coupling.F .
$ ln -s ../../ice_coupling.F SourceMods/src.csim/ice_coupling.F
```

Now update the use statements (Starting at line 23):
```
! !USES:                                                                                        
!                                                                                               
      use ice_kinds_mod
      use ice_model_size
      use ice_constants
      use ice_calendar
      use ice_grid
      use ice_state
      use ice_flux
      use ice_albedo
      use ice_mpi_internal
      use ice_timers
      use ice_fileunits
      use ice_work, only: worka, work_l1
#ifdef coupled
      use shr_sys_mod, only : shr_sys_flush
      use ice_history, only : runtype
      use cpl_contract_mod
      use cpl_interface_mod
      use cpl_fields_mod
      ! Add this line:
      use cpl_comm_mod
```

Now search for "cpl\_interface\_init" and replace it with "cpl\_comm\_init" (Starting at line 102)
```
     write(nu_diag,*) 'calling cpl_interface_init for model: ',
     &     in_model_name,' ', trim(cpl_fields_icename)

!      call cpl_interface_init(cpl_fields_icename,model_comm)
     call cpl_comm_init(cpl_fields_icename, model_comm)
```


**Editing forcing\_coupled.F**

First, get the source file:
```
$ cp ../../models/ocn/pop/source/forcing_coupled.F .
$ ln -s ../../forcing_coupled.F SourceMods/src.pop/forcing_coupled.F
```

Now update the use statements (Starting at line 13):
```
      use kinds_mod
      use model_size
      use domain
      use communicate
      use constants
      use io
      use stencils
      use time_management
      use grid
      use prognostic
      use exit_mod
      use ice
      use forcing_shf
      use forcing_sfwf
      use qflux_mod
      use ms_balance
      use timers
      use cpl_contract_mod
      use cpl_interface_mod
      ! Add this line!
      use cpl_init_mod
      use cpl_fields_mod
      use shr_sys_mod
      use registry
```

Now search for "cpl\_interface\_init" and replace it with "cpl\_comm\_init" (Starting at line 1509)
```
!      call cpl_interface_init(cpl_fields_ocnname,MPI_COMM_OCN)
      call cpl_comm_init(cpl_fields_ocnname, MPI_COMM_OCN)
      end subroutine ocn_coupling_setup
```


### Modifying the Fortran Coupler ###

If you would like to run CCSM with the Fortran coupler (perhaps to simply verify that CCSM runs properly with the new initialization routine, here's a quick description of the modifications necessary to force the CCSM coupler to use the new initialization.

**Get main.F90**
```
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ cp ../../models/cpl/cpl6/main.F90 .
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ ln -s ../../main.F90 SourceMods/src. 
src.cam   src.cpl   src.datm  src.dlnd  src.latm  src.xatm  src.xlnd  
src.clm   src.csim  src.dice  src.docn  src.pop   src.xice  src.xocn  
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ ln -s ../../main.F90 SourceMods/src.
src.cam   src.cpl   src.datm  src.dlnd  src.latm  src.xatm  src.xlnd  
src.clm   src.csim  src.dice  src.docn  src.pop   src.xice  src.xocn  
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ ln -s ../../main.F90 SourceMods/src.cpl/main.F90
steder@joxaren:~/ccsm3_0/scripts/pyccsm-dead-models$ ls -l SourceMods/src.cpl/
total 8
-rw-r--r--  1 steder users 152 Feb  8 23:10 README
lrwxrwxrwx  1 steder users  14 Feb  8 23:35 main.F90 -> ../../main.F90
```
**Add the use statement**
Lines 22-41:
```
program cpl

! !USES:                                                                                        

   use shr_date_mod       ! shr code: date/time data type & methods                             
   use shr_msg_mod        ! shr code: i/o redirect, chdir                                       
   use shr_sys_mod        ! shr code: wrappers to system calls                                  
   use shr_timer_mod      ! shr code: timing utilities                                          
   use shr_mpi_mod        ! shr code: mpi layer                                                 

   use cpl_kind_mod       ! kinds for strong typing                                             
   use cpl_mct_mod        ! wrapper/access to mct lib                                           
   use cpl_comm_mod       ! MPI process & communicator group info                               
   use cpl_fields_mod     ! fields lists common to cpl & components                             
   use cpl_contract_mod   ! contract (message-passing) data type & methods                      
   use cpl_interface_mod  ! contract/message-passing wrapper routines                           
   use cpl_infobuf_mod    ! information buffer module                                           
   use cpl_control_mod    ! control flags/logic subsystem module                                
   use cpl_map_mod        ! mapping subsystem module                                            
!  use cpl_iocdf_mod      ! netCDF file i/o routines (for debug data dumps)   
```

Simply add our updated module to this list:
```
program cpl

! !USES:                                                                                        

   use shr_date_mod       ! shr code: date/time data type & methods                             
   use shr_msg_mod        ! shr code: i/o redirect, chdir                                       
   use shr_sys_mod        ! shr code: wrappers to system calls                                  
   use shr_timer_mod      ! shr code: timing utilities                                          
   use shr_mpi_mod        ! shr code: mpi layer                                                 

   use cpl_kind_mod       ! kinds for strong typing                                             
   use cpl_mct_mod        ! wrapper/access to mct lib                                           
   use cpl_comm_mod       ! MPI process & communicator group info                               
   use cpl_fields_mod     ! fields lists common to cpl & components                             
   use cpl_contract_mod   ! contract (message-passing) data type & methods                      
   use cpl_interface_mod  ! contract/message-passing wrapper routines                           
   use cpl_infobuf_mod    ! information buffer module                                           
   use cpl_control_mod    ! control flags/logic subsystem module                                
   use cpl_map_mod        ! mapping subsystem module                                            
!  use cpl_iocdf_mod      ! netCDF file i/o routines (for debug data dumps)   
   use cpl_comm_mod
```

**Change the init call**
Lines 168-172:
```
   !----------------------------------------------------------------------------                
   ! determine cpl6's MPI communicator group, redirect stdin/out as necessary                   
   !----------------------------------------------------------------------------                
   call cpl_interface_init(cpl_fields_cplname,local_comm)
```
Changes to:
```
   !----------------------------------------------------------------------------                
   ! determine cpl6's MPI communicator group, redirect stdin/out as necessary                   
   !----------------------------------------------------------------------------                
   !call cpl_interface_init(cpl_fields_cplname,local_comm)
   call cpl_comm_init(cpl_fields_cplname, local_comm)
```

That's it!  Your Fortran Coupler should now be using the new initialization code.  All you should have to do is rebuild the case.

## PyCPL ##

### Dependencies ###

First we'll need to install [MaroonMPI](MaroonMPI.md) and [PyCDF](PyCDF.md).

Instructions for MaroonMPI can be found:
[On MaroonMPI's Google Code page](http://code.google.com/p/maroonmpi/)

Instructions for [PyCDF](PyCDF.md) can be found on our [PyCDF](PyCDF.md) page.

### Installation ###

Finally, we have to
  * download the PyCPL mapping data files for our specific CCSM resolution
  * Check, modify, and confirm the settings in cpl.nml
  * Update the logpath inside the main.py script.

**Download the PyCCSM Source**

You can download the source either through anonymous SVN or from the download page on this site.

Place the source in a convenient place on your system and unpack it (i.e.):
```
$ pwd
/home/steder/Sandbox
$ ls
pyccsm_r112_2-8-07.tar.gz
$ tar -zxf pyccsm_r112_2-8-07.tar.gz
$ cd pyccsm_r112_2-8-07
```

**Download the PyCPL Mapping Data from Earth Systems Grid**

Visit [http://www.earthsystemgrid.org/](http://www.earthsystemgrid.org/)

Look for a link to "Browse" the CCSM "Dataset Catalog" and follow the link.

Then select
CCSM 3.0 Release( CCSM 3.0 source code and input data )

And then:
CCSM 3.0 input data

Finally you'll see a listing of data files available for download.
If you know the specific resolution you're interested in running you can browse for and download only those files.  However, there's a simple option for users who may be interested in doing more general runs on all sorts of resolutions and that is the "coupler data files for all resolutions" link on this page.

Download and unpack this file and you should see:
```
$ tar -xf ccsm3.0.inputdata.cpl.tar
$ ls
ccsm3.0.inputdata.cpl.tar  inputdata/                 inputdata_user/
$ cd inputdata/cpl/cpl6/
$ ls
map_2x25d_to_gx1v3_aave_da_030226.nc    map_T85_to_gx1v3_aave_da_020405.nc
map_2x25d_to_gx1v3_bilin_da_030226.nc   map_T85_to_gx1v3_bilin_da_020405.nc
map_T31_to_gx3v5_aave_da_040122.nc      map_gx1v3_to_2x25d_aave_da_030226.nc
map_T31_to_gx3v5_bilin_da_040122.nc     map_gx1v3_to_T42_aave_da_010709.nc
map_T42_to_gx1v3_aave_da_010709.nc      map_gx1v3_to_T62_aave_da_010806.nc
map_T42_to_gx1v3_bilin_da_010710.nc     map_gx1v3_to_T85_aave_da_020405.nc
map_T42_to_gx3v5_aave_da_040122.nc      map_gx3v5_to_T31_aave_da_040122.nc
map_T42_to_gx3v5_bilin_da_040122.nc     map_gx3v5_to_T42_aave_da_040122.nc
map_T62_to_gx1v3_aave_da_010806.nc      map_gx3v5_to_T62_aave_da_040506.nc
map_T62_to_gx1v3_bilin_da_010806.nc     map_r05_to_gx1v3_roff_smooth_010718.nc
map_T62_to_gx3v5_aave_da_040506.nc      map_r05_to_gx3v5_e2000r500_040209.nc
map_T62_to_gx3v5_bilin_da_040506.nc
```

You'll need to copy the appropriate data files (or symbolic link them) to the PyCCSM data subdirectory.  Your data subdirectory should look something like:
```
-bash-2.05b$ ls -l
total 126524
-rwxr-xr-x    1 steder   mcsz          811 Oct 17 09:48 cpl.nml*
-rwxr-xr-x    1 steder   mcsz          817 Oct 17 09:54 cpl.stdin*
-rwxr-xr-x    1 steder   mcsz          118 Oct 17 09:49 cpl_stdio.nml*
-rwxr-xr-x    1 steder   mcsz      2343760 Oct 17 09:49 map_gx3v5_to_T31_aave_da_040122.nc*
-rwxr-xr-x    1 steder   mcsz     122436572 Oct 17 09:54 map_r05_to_gx3v5_e2000r500_040209.nc*
-rwxr-xr-x    1 steder   mcsz      2343760 Oct 17 09:49 map_T31_to_gx3v5_aave_da_040122.nc*
-rwxr-xr-x    1 steder   mcsz      2128656 Oct 17 09:49 map_T31_to_gx3v5_bilin_da_040122.nc*
```

**Update cpl.nml**

The most important variables to change here are:
  * map\_a2of\_fn: filename of the mapping between the atmosphere and ocean fluxes
  * map\_a2os\_fn: filename of the mapping between the atmosphere and ocean surfaces
  * map\_o2af\_fn: filename of the mapping datafile between the ocean and atmosphere fluxes
  * map\_r2o\_fn: filename of the mapping datafile between river and ocean data
```
  &inparm                                                                                       
  case_name     = 'deadtest '                                                                   
  case_desc     = 'deadtest deadtest '                                                          
  start_type    = "initial"                                                                     
  start_date    =  00010101                                                                     
  start_pfile   = 'rpointer.cpl '                                                               
  start_bfile   = "null"                                                                        
  stop_option   = 'ndays'                                                                       
  stop_n        =  2                                                                            
  rest_option   = 'ndays'                                                                       
  rest_n        =  5                                                                            
  diag_option   = 'ndays'                                                                       
  diag_n        =  10                                                                           
  hist_option   = 'never'                                                                       
  hist_n        =  -999                                                                         
  hist_64bit    =  .false.                                                                      
  avhist_option = "never"                                                                       
  avhist_n      =  -999                                                                         
  map_a2of_fn   = "map_T31_to_gx3v5_aave_da_040122.nc"                                          
  map_a2os_fn   = "map_T31_to_gx3v5_bilin_da_040122.nc"                                         
  map_o2af_fn   = "map_gx3v5_to_T31_aave_da_040122.nc"                                          
  map_r2o_fn    = "map_r05_to_gx3v5_e2000r500_040209.nc"                                        
  orb_year      =  1990                                                                         
  flx_epbal     = 'off'                                                                         
  flx_albav     = .false.                                                                       
  info_dbug     =  1                                                                            
  bfbflag       = .false.                                                                       
  /                                                                                             
```

The above variables have to be edited to point to the appropriate data files for the resolution of the models you've chosen.

**Update main.py script**

Finally you can modify the coupler script to change a number of behavior variables, choose the location for NetCDF log output, and decide the number of days you'd like the model to run.

Here's what a typical main.py will look like:
```
#!/usr/bin/env python2.4                                                                        
import os.path                                                                                  
import time                                                                                     
import traceback                                                                                
abspath = os.path.abspath(__file__)                                                             
basedir = os.path.dirname(abspath)                                                              
ROOT=basedir                                                                                    
                                                                                                
import coupler                                                                                  
                                                                                                
def main():                                                                                     
    print "pyCPL: Start!"                                                                       
    # Coupler Path, doComputations, doMapping, initializeRiverMapping,Logging)                  
    coupler.TIMER.start("total")                                                                
    coupler.TIMER.start("init")                                                                 
    print "pyCPL: Coupler root is",ROOT                                                         
    LOGROOT= os.path.join("/home/steder/exe/pyccsm/nc",str(time.ctime()).replace(" ",""))       
    print "pyCPL: NC Output path is",LOGROOT                                                    
    c = coupler.Coupler(os.path.join(ROOT,"cpl","cpl.nml"),True,True,True,True,LOGROOT)         
                                                                                                
    coupler.TIMER.stop("init")                                                                  
    c.run(5)                                                                                    
    coupler.TIMER.stop("total")                                                                 
    #print "StartDate:",cpl.date                                                                
    coupler.cpl.comm.finalize()                                                                 
    print "pyCPL: Done!"                                                                        
    print coupler.TIMER                                                                         
                                                                                                
if __name__=="__main__":                                                                        
    main()                                                                                      
```

This single line instantiates a coupler and sets a number of values for the length of the run:
```
    c = coupler.Coupler(os.path.join(ROOT,"cpl","cpl.nml"),True,True,True,True,LOGROOT)
```

The arguments to Coupler are as follows:
Coupler( "/path/to/cpl.nml", DoComputations, DoMapping, DoRiverInitialization, DoLogging, "/path/to/output/logfiles/")

The boolean values:
  * DoComputations: turn off calculations before/after send/recv operations
  * DoMapping: turn off mapping calculations
  * DoRiverInitialization: the river initialization for some resolutions is extremely long.  This lets you turn it off to speed debugging/testing.
  * DoLogging:  Enables a bunch of NetCDF output of data at various points in the coupler.  Useful for debugging.

The other interesting thing about this file is the _run_ method which allows us to set the number of model days we'd like to simulate:
```
    c.run(5)                              
```

Simply change that 5 to the number of days you'd like to run.  At this time, the _run_ method doesn't support the ability to run for a fractional number of days (i.e. run a variable number of steps rather then days) but modifying _run_ to allow this behavior is pretty trivial.

# Running PyCCSM #

## Modifying CCSM Run Scripts ##

### Make a copy of one of your CCSM run scripts ###

Imagine we're working with that pyccsm-dead-models case we were using earlier:
```
$ pwd
/home/steder/ccsm3_0/scripts/pyccsm-dead-models
$ ls
Buildexe/           configure*  env_mach.jazz*                  pyccsm-dead-models.jazz.run*
Buildlib/           dead.F90    env_run
Buildnml_Prestage/  env.readme  main.F90
SourceMods/         env_conf    pyccsm-dead-models.jazz.build*
$ cp pyccsm-dead-models.jazz.run custom.run
$ ls
Buildexe/           SourceMods/  dead.F90    env_mach.jazz*  pyccsm-dead-models.jazz.build*
Buildlib/           configure*   env.readme  env_run         pyccsm-dead-models.jazz.run*
Buildnml_Prestage/  custom.run*  env_conf    main.F90
```

### Editing the copy ###

Let's take a look inside custom.run.  We're looking specifically at lines 82-94:
```
# -------------------------------------------------------------------------                     
# Run the model                                                                                 
# -------------------------------------------------------------------------                     

env | egrep '(MP_|LOADL|XLS|FPE|DSM|OMP|MPC)' # document env vars                               

cd $EXEROOT/all
paste ${PBS_NODEFILE} mpirun.pgfile1 > mpirun.pgfile
echo "`date` -- CSM EXECUTION BEGINS HERE"
mpirun -pg mpirun.pgfile ./$COMPONENTS[1]
wait
echo "`date` -- CSM EXECUTION HAS FINISHED"
```

**Assuming we're using MPICH2**

First we'll comment out the above lines between "cd $EXEROOT/all" and "wait"
```
# -------------------------------------------------------------------------                     
# Run the model                                                                                 
# -------------------------------------------------------------------------                     

env | egrep '(MP_|LOADL|XLS|FPE|DSM|OMP|MPC)' # document env vars                               

cd $EXEROOT/all
##paste ${PBS_NODEFILE} mpirun.pgfile1 > mpirun.pgfile
##echo "`date` -- CSM EXECUTION BEGINS HERE"
##mpirun -pg mpirun.pgfile ./$COMPONENTS[1]
wait
echo "`date` -- CSM EXECUTION HAS FINISHED"
```

Now we'll add our new MPICH2 run code:
  * We need to start MPD on the nodes in _${PBS\_NODEFILE}_
```
# -------------------------------------------------------------------------                     
# Run the model                                                                                 
# -------------------------------------------------------------------------                     

env | egrep '(MP_|LOADL|XLS|FPE|DSM|OMP|MPC)' # document env vars                               

cd $EXEROOT/all
##paste ${PBS_NODEFILE} mpirun.pgfile1 > mpirun.pgfile
##echo "`date` -- CSM EXECUTION BEGINS HERE"
##mpirun -pg mpirun.pgfile ./$COMPONENTS[1]
mpdboot -n `wc -l < $PBS_NODEFILE` -f $PBS_NODEFILE
wait
echo "`date` -- CSM EXECUTION HAS FINISHED"
```

  * We can now go ahead and add the MPICH2 equivalent of mpirun (This example runs all the Fortran components, including the Fortran coupler using MPICH2):
```
# -------------------------------------------------------------------------                     
# Run the model                                                                                 
# -------------------------------------------------------------------------                     

env | egrep '(MP_|LOADL|XLS|FPE|DSM|OMP|MPC)' # document env vars                               

cd $EXEROOT/all
##paste ${PBS_NODEFILE} mpirun.pgfile1 > mpirun.pgfile
##echo "`date` -- CSM EXECUTION BEGINS HERE"
##mpirun -pg mpirun.pgfile ./$COMPONENTS[1]
mpdboot -n `wc -l < $PBS_NODEFILE` -f $PBS_NODEFILE
mpiexec -l -n $NTASKS[1] $COMPONENTS[1] : -n $NTASKS[2] $COMPONENTS[2] : -n $NTASKS[3] $COMPONE
NTS[3] : -n $NTASKS[4] $COMPONENTS[4] : -n $NTASKS[5] $COMPONENTS[5]
wait
echo "`date` -- CSM EXECUTION HAS FINISHED"
```

  * This example runs the Python Coupler with the Fortran Components using MPICH2:
( $PYTHON must be defined as a variable that points to your Python interpreter and $COUPLER\_PY must be defined as the PyCCSM/main.py script )
```
# -------------------------------------------------------------------------                     
# Run the model                                                                                 
# -------------------------------------------------------------------------                     

env | egrep '(MP_|LOADL|XLS|FPE|DSM|OMP|MPC)' # document env vars                               

cd $EXEROOT/all
##paste ${PBS_NODEFILE} mpirun.pgfile1 > mpirun.pgfile
##echo "`date` -- CSM EXECUTION BEGINS HERE"
##mpirun -pg mpirun.pgfile ./$COMPONENTS[1]
mpdboot -n `wc -l < $PBS_NODEFILE` -f $PBS_NODEFILE
mpiexec -l -n 1 $PYTHON $COUPLER_PY : -n $NTASKS[2] $COMPONENTS[2] : -n $NTASKS[3] 
$COMPONENTS[3] : -n $NTASKS[4] $COMPONENTS[4] : -n $NTASKS[5] $COMPONENTS[5]
wait
echo "`date` -- CSM EXECUTION HAS FINISHED"
```
### NOTE ###
Very dependent on your specific computing environment.

| PyCCSM |
|:-------|
