| [InstallPyMCT](InstallPyMCT.md) | [StageOneInstallPyMCT](StageOneInstallPyMCT.md) |
|:--------------------------------|:------------------------------------------------|

# Stage 0 #

## MPICH ##
The only configuration options that have to be set are those flags to enable Fortran 90 support for your compiler
and another single flag to enable support for shared (dynamic) MPI libraries.

### Fortran ###

For compatibility reasons with other Fortran compilers it is strongly encouraged that when you build MPI for your Fortran compiler that you specify the linkage conventions for symbols.

| Compiler | Convention | Example |
|:---------|:-----------|:--------|
| G77 | All Symbols are lowercase with an appended underscore | subroutine ComputeFlux -> `computeflux_`|
| Absoft | Defaults to Mixed Case, No Underscore | subroutine ComputeFlux -> `ComputeFlux` |
| PGI | Same conventions as G77 | subroutine ComputeFlux -> `computeflux_` |
| Intel | Same conventions as G77 | subroutine ComputeFlux -> `computeflux_` |

So if you perchance have to use the Absoft F90 compiler, or some other compiler with differing linkage conventions it's important that you look into compiler options to modify those linkage conventions to make them compatible with other compilers.

The consequences of not doing this are building specialized versions of nearly every system library for each compiler which uses a different convention.  Generally this is far more work then simply changing the compilers default conventions to match the other compilers on your system.

**Flags for the Absoft Fortran Compilers**
```
FFLAGS=-f -N15                                                                    
FCFLAGS=-YEXT_NAMES=LCS -YEXT_SFX=_                                               
F90FLAGS=-YEXT_NAMES=LCS -YEXT_SFX=_    
```
Set these environment variables in your shell before continuing with the configuration if you can.

Alternatives:
1) you can set these flags everywhere you invoke Absoft
2) you can build alternative libraries for Absoft and other compilers.

We think most installations with more than one Fortran compiler will prefer to set the flags globally so as to avoid these problems.

### Configuration ###

Unpack the mpich source into a convenient directory and move to where the MPICH _configure_ script is located.  For example:
```
-bash-2.05b$ cd mpich2-1.0.2p1/
-bash-2.05b$ pwd
/home/steder/mpich2-1.0.2p1
-bash-2.05b$ ls
bin/            COPYRIGHT.rtf*    mpich2.sln*           README.romio*
CHANGES*        doc/              mpich2s.vcproj*       README.testing*
conf19599.sh*   examples/         mpich2s.vs05.vcproj*  README.winbin.rtf*
confdb/         lib/              mpich2.vcproj*        README.windeveloper*
config.log*     maint/            mpich2.vs05.sln*      README.windows*
config.status*  Makefile*         mpich2.vs05.vcproj*   RELEASE_NOTES*
config.system*  Makefile.in*      mpi.def*              src/
configure*      Makefile.sm*      mpi.vcproj*           test/
configure.in*   makewindist.bat*  mpi.vs05.vcproj*      winconfigure.wsf*
CONFTEST.mod*   man/              README*               www/
COPYRIGHT*      mpich2.def*       README.developer*
```

The first  issue we want to make sure we address is enabling dynamic/shared library support.  For MPICH1/MPICH2 this can be accomplished by adding the following flags to the configure line.

For example:
```
-bash-2.05b$ ./configure --enable-sharedlibs=gcc
```
OR
```
-bash-2.05b$ ./configure --enable-sharedlibs=osx-gcc
```

The second issue is regarding the Fortran 90 flags we picked out to change the Fortran compilers linking conventions.  Whether are environment is updated or not we can simply specify the flags right on the configure line:
```
-bash-2.05b$ ./configure --enable-sharedlibs=osx-gcc  --enable-f90 FCFLAGS="-f -N15" F90FLAGS="-YEXT_NAMES=LCS -YEXT_SFX=_"                 
```

Assuming FFLAGS, FCFLAGS, and F90FLAGS are set in our environment our final configure line can be just the following.  Additional options, for instance device or debugging support are perfectly acceptable to add if desired.
```
-bash-2.05b$ ./configure --enable-sharedlibs=gcc --enable-f90
```

| [InstallPyMCT](InstallPyMCT.md) | [StageOneInstallPyMCT](StageOneInstallPyMCT.md) |
|:--------------------------------|:------------------------------------------------|