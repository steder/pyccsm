| [PyMCT](PyMCT.md)  |
|:-------------------|
| [StageZeroInstallPyMCT](StageZeroInstallPyMCT.md) | [StageOneInstallPyMCT](StageOneInstallPyMCT.md) |

# Introduction #
This goal of this document is to walk through the steps to configure, build, and install PyMCT and setup an environment for writing and running MCT programs in Python.

# Requirements #
### Compilers and Languages ###
We expect you to have:
  * C Compiler
  * Fortran 90 Compiler (Portland Group Fortran is incompatible 

&lt;link&gt;

!)
  * Java JDK
  * Python 2.4 or later

### Libraries ###
  * MPI (We prefer MPICH2, but others should work)
  * MPI with Dynamic Library support enabled (--enable-sharedlibs configuration option)
  * NetCDF
If you do not know if you have any of the above software installed and properly configured on the system you intend to use for PyMCT you will need to install or reconfigure the missing or misconfigured packages before you can begin the installation of PyMCT.

# Installation Overview #
## Stage 0 ##
Stage 0 is useful for administrators looking to setup their machine to use PyMCT from a fresh installation.  If you do not have  necessary compilers or basic libraries you should start here. [StageZeroInstallPyMCT](StageZeroInstallPyMCT.md)
## Stage 1 ##
Initial Configuration of the PyMCT environment [StageOneInstallPyMCT](StageOneInstallPyMCT.md)

## Stage 2 ##

### MCT ###
This stage builds and installs the Fortran version of MCT.
[ModelCouplingToolkitInstall](ModelCouplingToolkitInstall.md)

## Stage 3 ##
### Babel Prerequisites ###
This stage adds the basic modules required by Babel and PyMCT.  These are straightforward installations.
  * Chasm
  * Python Array Modules
[ChasmInstall](ChasmInstall.md)
[NumpyInstall](NumpyInstall.md)
[LibXML](LibXML.md)

### Babel ###
This stage builds the Babel Compiler.  This is a large and sophisticated application.
[BabelInstall](BabelInstall.md)

## Stage 4 ##
### PyMCT ###
Here we use Babel to update the bindings.  We then modify makefiles to fit our system and compile and install the MCT Server and Client libraries.
[ServerInstall](ServerInstall.md)
[ClientInstall](ClientInstall.md)

### Utilities ###
This is our final stage which consists of just installing some useful utility packages.
[UtilitiesInstall](UtilitiesInstall.md)
[MaroonInstall](MaroonInstall.md)

| [PyMCT](PyMCT.md)  |
|:-------------------|
| [StageZeroInstallPyMCT](StageZeroInstallPyMCT.md) | [StageOneInstallPyMCT](StageOneInstallPyMCT.md) |

