| [InstallPyMCT](InstallPyMCT.md) | [ClientInstall](ClientInstall.md) |
|:--------------------------------|:----------------------------------|

# PyMCT Server Installation #

## Obtaining the PyMCT Babel Bindings ##

You will have to download the Babel cross-language bindings for MCT from the [MCT Website](http://www-unix.mcs.anl.gov/mct/).  Here's a direct link to the [bindings tarball](ftp://ftp.mcs.anl.gov/pub/acpi/MCT/babel.tar.gz).

Once you have the sources downloaded, place them in your _Sandbox_ directory OR _wherever you placed your MCT **source**_.  For example:

```
$ pwd
/home/steder/Sandbox
$ ls
babel.tar.gz  MCT/  MCTnodata.tar.gz
```

Now go ahead and decompress the babel.tar.gz file and move it into the MCT subdirectory:

```
$ ls
babel.tar.gz  MCT/  MCTnodata.tar.gz
$ tar -zxf babel.tar.gz 
$ ls 
babel/  babel.tar.gz  MCT/  MCTnodata.tar.gz
$ mv babel MCT/
$ ls
babel.tar.gz  MCT/  MCTnodata.tar.gz
$ cd MCT
$ ls
aclocal.m4  configure.ac  examples/    Makefile.conf.in  mpeu/        README
babel/      COPYRIGHT     install-sh*  mct/              mpi-serial/
configure*  doc/          Makefile     mkinstalldirs*    protex/
$ cd babel/
$ ls
example-c++/     genscl.py  Makefile        mctserver/  README
example-python/  impls/     mctclient-c++/  MCT.sidl
```

## Configuring ##

First we'll copy a few implementation files to the mctserver directory.
```
$ cp impls/*.F90 mctserver/
```

From this directory you can regenerate the babel bindings with your specific version of Babel by running the following command:

```
$ $PYMCT_ROOT/$BABEL/bin/babel -sf90 MCT.sidl -o mctserver/
```

This should take a few minutes to finish processing.

## Editing the Makefiles ##

Space reserved.

## Building the PyMCT Server Library ##

```
$ pwd 
/home/steder/Sandbox/MCT/babel
$ ls
example-c++/     genscl.py  Makefile        mctserver/  README
example-python/  impls/     mctclient-c++/  MCT.sidl
$ cd mctserver/
$ make
```

```
$ make install
```
OR all that really has to be done to complete the installation is copy this library to the appropriate place:
```
$ cp libmctserver.so $PYMCT_ROOT/lib/
```

| [InstallPyMCT](InstallPyMCT.md) | [ClientInstall](ClientInstall.md) |
|:--------------------------------|:----------------------------------|


