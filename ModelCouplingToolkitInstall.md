| [StageOneInstallPyMCT](StageOneInstallPyMCT.md) | [InstallPyMCT](InstallPyMCT.md) |
|:------------------------------------------------|:--------------------------------|

# Stage 2 #
# Installing the Model Coupling Toolkit #

## Step 1:  Obtaining MCT ##

### Downloading the source ###

You can download the most recent version of the Model Coupling Toolkit directly from
http://www-unix.mcs.anl.gov/mct/

You only need the source so the following link [MCT source only (.tar.gz). Only 0.3MB](ftp://ftp.mcs.anl.gov/pub/acpi/MCT/MCTnodata.tar.gz) is all you need.

### Unpacking ###
Let's create a directory to hold our work.  In this case, we're using a subdirectory named Sandbox in our user directory.

```
$ pwd
/home/steder/Sandbox
$ ls
MCTnodata.tar.gz
$ tar -zxf MCTnodata.tar.gz 
$ ls
MCT/  MCTnodata.tar.gz
$ cd MCT
$ ls
aclocal.m4    COPYRIGHT  install-sh*       mct/            mpi-serial/
configure*    doc/       Makefile          mkinstalldirs*  protex/
configure.ac  examples/  Makefile.conf.in  mpeu/           README
```

## Step 2:  Configuration and Installation ##

Now, assuming we've followed the instructions in [StageOneInstallPyMCT](StageOneInstallPyMCT.md) and properly configured our environment variables we should be all set to install MCT.

### ./configure ###

This is where we need to be to execute the following commands:
```
$ pwd
/home/steder/Sandbox/MCT
$ ls
aclocal.m4    COPYRIGHT  install-sh*       mct/            mpi-serial/
configure*    doc/       Makefile          mkinstalldirs*  protex/
configure.ac  examples/  Makefile.conf.in  mpeu/           README
```

Now we're going to run configure with options and parameters we setup in Stage 1.

```
$ ./configure --prefix=$PYMCT_ROOT --enable-babel   FC=$FC     F90=$F90    MPIF90=$F90 
    BABELROOT=$PYMCT_ROOT/$BABEL/     PYTHON=$PYTHON       COMPILER_ROOT=$F90_COMPILER_ROOT
    FCFLAGS=$F90FLAGS      F90FLAGS=$F90FLAGS
```

Note that the above is a single line.

--prefix determines where we'll install what we're building here when we type "make install".

--enable-babel turns on various MCT build options required for compatibility with Babel types.

FC, F90, MPIF90 just force MCT to use our _mpif90_ compiler everywhere.  This is mostly useful for consistencies sake.

BABELROOT is used later on in Stage 4.  For now we're just setting it.

PYTHON is simply the path to the python interpreter.  Luckily, we set it to a variable of the same name in our environment.

COMPILER\_ROOT is used by MCT to get access to some other files that are usually located relative to the compilers root.

FCFLAGS and F90FLAGS are there to enforce consistent linking conventions.  You can read more about this and why it's important on [StageZeroInstallPyMCT](StageZeroInstallPyMCT.md).

### make ###

Immediately after running configure we should be able to run:

```
$ make 
```

The last few lines of output of the make command should look something like this:
```
rm -f libmct.a
ar cq libmct.a m_MCTWorld.o m_AttrVect.o m_GlobalMap.o m_GlobalSegMap.o m_GlobalSegMapComms.o m_Accumulator.o m_SparseMatrix.o m_Navigator.o m_AttrVectComms.o m_AttrVectReduce.o m_AccumulatorComms.o m_GeneralGrid.o m_GeneralGridComms.o m_SpatialIntegral.o m_SpatialIntegralV.o m_MatAttrVectMul.o m_Merge.o m_GlobalToLocal.o m_ExchangeMaps.o m_ConvertMaps.o m_SparseMatrixDecomp.o m_SparseMatrixToMaps.o m_SparseMatrixComms.o m_SparseMatrixPlus.o m_Router.o m_Rearranger.o m_Transfer.o          
make[2]: Leaving directory `/home/steder/Sandbox/MCT/mct'
make[1]: Leaving directory `/home/steder/Sandbox/MCT/'
```

### make install ###

```
$ make install 
```

After this you should see quite a few lines of output as all the files we just built are copied into  $PYMCT\_ROOT.

## Lonestar ##

Configuring MCT on lonestar is a little simpler as you shouldn't need to pass any flags to the Intel compiler.  Configure looks like this:

```
$ ./configure --prefix=$PYMCT_ROOT --enable-babel   FC=$FC     F90=$F90    MPIF90=$F90 
    BABELROOT=$PYMCT_ROOT/$BABEL/     PYTHON=$PYTHON       COMPILER_ROOT=$F90_COMPILER_ROOT
```

| [StageOneInstallPyMCT](StageOneInstallPyMCT.md) | [InstallPyMCT](InstallPyMCT.md) |
|:------------------------------------------------|:--------------------------------|