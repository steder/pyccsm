| [StageZeroInstallPyMCT](StageZeroInstallPyMCT.md) |  [InstallPyMCT](InstallPyMCT.md) | [ModelCouplingToolkitInstall](ModelCouplingToolkitInstall.md) |
|:--------------------------------------------------|:---------------------------------|:--------------------------------------------------------------|

# Stage 1 #

Here we're concerned with simply setting up environment variables to make the installation of PyMCT straightfoward.  Take some time and research what these values should be for your system.

## Setting Environment Variables ##

If you're unfamiliar with setting shell environment variables we can provide a little information on how to make these modifications yourself.  However, first we need to identify which shell you are using.

### Which shell am I running? ###
A simple way to check which shell you're running is to simply type:
```
echo $SHELL
```
at the prompt.

The result will likely be either of these:
```
/bin/bash
/bin/csh
```

### BASH ###
To setup environment variables on a system using BASH you will want to edit a hidden file like:
```
/home/user/.bash_profile
/home/user/.bashrc
```

In Bash setting an environment variable looks like:
```
VARIABLE = value
export VARIABLE
```

For instance
```
PATH=/usr/bin:/bin
export PATH
```

### CSHELL ###
To setup environment variables on a system using CSHELL you will want to edit a hidden file like:
```
/home/user/.cshrc
```

In CShell setting an environment variable looks like:
```
setenv VARIABLE "value"
```

For instance:
```
setenv PATH "/bin:/sbin:/usr/bin:/usr/sbin"
```

### Other Systems ###

There are many different shells available and it's very difficult to cover every possible system here.  We strongly recommend that you contact a system administrator for your system and take a look at the documentation for your shell.  You only really need to know how to set environment variables.  Once you have that information we can proceed.

## Variables to set ##

```
PYMCT_ROOT
```
```
F90_COMPILER_ROOT
MPI_COMPILER_ROOT
MPIHEADER
MPILIBS
NON_MPI_F77
NON_MPI_F90
MPI_HOME
```

```
FFLAGS
FCFLAGS
F90FLAGS
FC
F77
F90
```

```
BABEL
PATH
LD_LIBRARY_PATH
LD_LOAD_LIBRARY
INCLUDE_PATH
SIDL_DLL_PATH
SIDL_DEBUG_DLOPEN
```

```
PYTHONPATH
PYTHON
CHASM_FORTRAN_VENDOR
```

## Example Values ##

### Installation directory ###

PYMCT\_ROOT is used to determine where all the necessary PyMCT packages and modules will be installed.  Possible directories might be "/usr/local" for a system-wide installation or a "/home/user" to install PyMCT for use by a single person.  Let's use the following to install PyMCT for a single user.
```
PYMCT_ROOT=$HOME/pymct/
```

### Compiler Flags ###
These values are talked about in more detail on the StageZeroInstallPyMCT page.  Hopefully they are already set for you.  You may be interested in adding to them to specify a certain optimization level, but assuming they are set you can leave them alone.

If they are not set, you need to look into the documentation for your compiler, or speak with your system admin, to determine what linking conventions are in use for your Fortran compiler.  Some compilers will require you to make changes, many will simply work _out of the box_ without setting these flags.

The Absoft compiler is a notable compiler that has very different linking conventions then most other Fortran compilers.  These flags change the default behavior of Absoft Fortran to match that of a Gnu, Intel, and Portland Group compilers.
```
FFLAGS=-f -N15
FCFLAGS=-YEXT_NAMES=LCS -YEXT_SFX=_
F90FLAGS=-YEXT_NAMES=LCS -YEXT_SFX=_
```

### Default compilers ###

These are shortcuts to the compilers we'll be using most.

```
FC=mpif90
F77=mpif77
F90=mpif90
```

### Compilers  and Paths for MCT ###

This next group of variables are used to setup the compiler paths needed by MCT.
F90\_COMPILER\_ROOT can be found by running the following command:
```
-bash-2.05b$ mpif90 -show
f90 -YEXT_NAMES=LCS -YEXT_SFX=_ -YEXT_NAMES=LCS -YEXT_SFX=_ -I/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/include -p/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/include -L/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/lib -lmpichf90 -X -rpath -X /soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/lib -lmpichf90 -lmpich -lpthread -lpvfs -lrt
```

You'll notice that f90 is the compiler itself.  You might see _gfortran_, _ifcbin_, _pgf90_ or some other compiler at the start of that line and if that's the case you'll need to modify the next command to match.
```
-bash-2.05b$ which f90
/soft/com/packages/absoft-9.0/opt/absoft/bin/f90
```

Now we can set F90\_COMPILER\_ROOT by dropping /bin/f90 from the path:
```
F90_COMPILER_ROOT=/soft/com/packages/absoft-9.0/opt/absoft/
```

MPI\_COMPILER\_ROOT can be determined a little more simply:
```
-bash-2.05b$ which mpif90
/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/bin/mpif90
```

Now we can set MPI\_COMPILER\_ROOT by dropping the /bin/mpif90 from the path:
```
MPI_COMPILER_ROOT=/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/
```

MPIHEADER and MPILIBS are set from MPI\_COMPILER\_ROOT.

So, these variables in our example are set to:
```
F90_COMPILER_ROOT=/soft/com/packages/absoft-9.0/opt/absoft/
MPI_COMPILER_ROOT=/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/
MPIHEADER=$MPI_COMPILER_ROOT/include/mpif.h
MPILIBS=$MPI_COMPILER_ROOT/lib
MPI_HOME=$MPI_COMPILER_ROOT
```

### Compiler paths for Babel ###
These variables are used to specify the Fortran compilers that should be used by Babel.  We don't need MPI libraries and includes at this point so we can use the original compilers directly.

You can find out the names of these compilers by typing
```
-bash-2.05b$ mpif77 -show
f77 -f -N15 -f -B108 -I/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/include -L/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/lib -Wl,-rpath -Wl,/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/lib -lmpich -lpthread -lpvfs -lrt
-bash-2.05b$ which f77
/soft/com/packages/absoft-9.0/opt/absoft/bin/f77
```

And:
```
-bash-2.05b$ mpif90 -show
f90 -YEXT_NAMES=LCS -YEXT_SFX=_ -YEXT_NAMES=LCS -YEXT_SFX=_ -I/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/include -p/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/include -L/soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/lib -lmpichf90 -X -rpath -X /soft/apps/packages/mpich2-sock-1.0.4p1-absoft-9.0-1-shared/lib -lmpichf90 -lmpich -lpthread -lpvfs -lrt
-bash-2.05b$ which f90
/soft/com/packages/absoft-9.0/opt/absoft/bin/f90
```

The output of _which_ gives you the values for these two variables.  Go ahead and set them in your environment.
```
NON_MPI_F77=/soft/com/packages/absoft-9.0/opt/absoft/bin/f77
NON_MPI_F90=/soft/com/packages/absoft-9.0/opt/absoft/bin/f90
```

### Babel and PyMCT library settings ###
First we setup BABEL to "babel"+"-"+"version number".  This should be the same as the directory name for the babel distribution.  This variable is used to set a number of other variables and is necessary for properly locating and running Babel.

PATH is set to include the executables installed in later stages and the babel executables.  We want to simply prepend our $PYMCT\_ROOT/bin and $PYMCT\_ROOT/$BABEL/bin directories to the $PATH variable.  Depending on your shell there may be better ways to do this.  See AppendingToEnvironmentVariables for more details.

LD\_LIBRARY\_PATH and LD\_LOAD\_LIBRARY are system paths referenced whenever you start an application and are used to look up dynamic libraries at runtime.  If these are not set properly you'll see a number of linking errors when attempting to use PyMCT.  Basically the directories containing the PyMCT client and server libraries need to be included in both LD\_LIBRARY\_PATH and LD\_LOAD\_LIBRARY.

INCLUDE\_PATH is modified primarily to include Python and Babel header files.

SIDL\_DLL\_PATH must point to the directory where the PyMCT server library is installed.  In this case we'll install it to $PYMCT\_ROOT/lib/.  The reasoning here is similar to that for LD\_LIBRARY\_PATH and LD\_LOAD\_LIBRARY.

SIDL\_DEBUG\_DLOPEN can be set equal to either 0 or 1.  If it's set to 1 we enable some extra debugging statements whenever PyMCT is imported.  Once everything is properly installed and working this could be set to 0, but until then leave it 1.

```
BABEL=babel-0.11.2
PATH=$PYMCT_ROOT/bin:$PYMCT_ROOT/$BABEL/bin:$PATH
LD_LIBRARY_PATH=$PYMCT_ROOT/lib:$PYMCT_ROOT/$BABEL/lib:$MPILIBS:$LD_LIBRARY_PATH
LD_LOAD_LIBRARY=$LD_LIBRARY_PATH:$LD_LOAD_LIBRARY
INCLUDE_PATH=$PYMCT_ROOT/include/python2.4/:$PYMCT_ROOT/include/:$PYMCT_ROOT
/$BABEL/include:$MPIHEADER:$INCLUDE_PATH
SIDL_DLL_PATH=$PYMCT_ROOT/lib/
SIDL_DEBUG_DLOPEN=1
```

### Python Settings ###
PYTHON is simply a pointer to the Python executable you wish to use for all installations.

PYTHONPATH needs to be updated to include the site-packages that will be installed in $PYMCT\_ROOT and in $PYMCT\_ROOT/$BABEL.

```
PYTHON=$PYMCT_ROOT/bin/python
PYTHONPATH=$PYMCT_ROOT/lib/python2.4/site-packages/:$PYMCT_ROOT/$BABEL/lib/python2.4/site-packages:./
```

### Chasm Settings ###

CHASM\_FORTRAN\_VENDOR needs to be set to one of the following options:
| Absoft | Alpha | Cray | GNU | IBMXL | Intel | Intel\_7 | Lahey | MIPSpro | NAG | SUNWspro |
|:-------|:------|:-----|:----|:------|:------|:---------|:------|:--------|:----|:---------|

Note that Portland Group Fortran is not on that list.  This is because the vendor has refused to make their array interface public.
```
CHASM_FORTRAN_VENDOR=Absoft
```

### Java Settings ###

The Babel installation requires that the CLASSPATH environment variable be set to ".".

```
CLASSPATH=.
```

# Lonestar (TACC) Specific Settings #

Lonestar uses CSH as its default shell.

You aren't supposed to modify _.cshrc_, instead you're supposed to create and modify a file named _.cshrc\_user_.  SO:

```
% pwd
/home/utexas/ig/steder
% touch .cshrc_user
% emacs .cshrc_user
```

Contents of _.cshrc\_user_
```
# PYMCT_ROOT

setenv PYMCT_ROOT "$HOME/pymct"

# F90_COMPILER_ROOT

setenv F90_COMPILER_ROOT "/opt/intel/compiler9.1/fc/"

# MPI_COMPILER_ROOT

setenv MPI_COMPILER_ROOT "/opt/MPI/intel9/mvapich-gen2/0.9.8/"

# MPIHEADER

setenv MPIHEADER "$MPI_COMPILER_ROOT/include/mpif.h"

# MPILIBS

setenv MPILIBS "$MPI_COMPILER_ROOT/lib"

# MPI_HOME
setenv MPI_HOME "$MPI_COMPILER_ROOT"

# NON_MPI_F77
setenv NON_MPI_F77 "/opt/intel/compiler9.1//fc/bin/ifort"
# NON_MPI_F90
setenv NON_MPI_F90 "/opt/intel/compiler9.1//fc/bin/ifort"

# FFLAGS
# FCFLAGS
# F90FLAGS

# FC
# F77
# F90
setenv FC "mpif90"
setenv F77 "mpif77"
setenv F90 "mpif90"

# PYVERSION
setenv PYVERSION "python2.5"
# PYTHON
setenv PYTHON "$PYMCT_ROOT/bin/python"

# BABEL
setenv BABEL "babel-1.0.2"

# PYTHONPATH
setenv PYTHONPATH "$PYMCT_ROOT/lib/$PYVERSION/site-packages/:$PYMCT_ROOT/$BABEL/lib
/$PYVERSION/site-packages/:./"

# CHASM_FORTRAN_VENDOR
setenv CHASM_FORTRAN_VENDOR "Intel"

# JAVA - CLASSPATH
setenv CLASSPATH "."

# PATH
setenv PATH "$PYMCT_ROOT/bin:$PYMCT_ROOT/$BABEL/bin:${PATH}"

# LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH ${PYMCT_ROOT}/lib:${PYMCT_ROOT}/${BABEL}/lib:${MPILIBS}:${LD
_LIBRARY_PATH}

# LD_LOAD_LIBRARY
setenv LD_LOAD_LIBRARY "${LD_LIBRARY_PATH}"

# INCLUDE_PATH
setenv INCLUDE_PATH "${PYMCT_ROOT}/include/${PYVERSION}/:${PYMCT_ROOT}/include/:${P
YMCT_ROOT}/${BABEL}/include:${MPIHEADER}"

# SIDL_DLL_PATH
setenv SIDL_DLL_PATH "${PYMCT_ROOT}/lib/"

# SIDL_DEBUG_DLOPEN
setenv SIDL_DEBUG_DLOPEN "1"

```

| [StageZeroInstallPyMCT](StageZeroInstallPyMCT.md) |  [InstallPyMCT](InstallPyMCT.md) | [ModelCouplingToolkitInstall](ModelCouplingToolkitInstall.md) |
|:--------------------------------------------------|:---------------------------------|:--------------------------------------------------------------|
