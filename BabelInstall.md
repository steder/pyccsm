| [InstallPyMCT](InstallPyMCT.md) | [ServerInstall](ServerInstall.md) |
|:--------------------------------|:----------------------------------|

# Babel Installation #

## About Babel ##
Stolen blatently from http://www.llnl.gov/casc/components/babel.html
Babel (pronounced babble) addresses the language interoperability problem using Interface Definition Language (IDL) techniques. An IDL describes the calling interface (but not the implementation) of a particular software library. IDL tools such as Babel use this interface description to generate glue code that allows a software library implemented in one supported language to be called from any other supported language. We have designed a Scientific Interface Definition Language (SIDL) that addresses the unique needs of parallel scientific computing. SIDL supports complex numbers and dynamic multi-dimensional arrays as well as parallel communication directives that are required for parallel distributed components. SIDL also provides other common features that are generally useful for software engineering, such as enumerated types, symbol versioning, name space management, and an object-oriented inheritance model similar to Java.

## Step 1:  Obtaining the Babel source ##

Babel can be downloaded directly from:  http://www.llnl.gov/casc/components/software.html

You will want to download the complete tarball and not the "runtime only" tarball.
[Complete link to v1.0.2 complete](http://www.llnl.gov/casc/components/docs/babel-1.0.2.tar.gz)

Any babel version later then 0.11.x should work.

```
$ pwd
/home/steder/Sandbox
$ ls
babel-1.0.2.tar.gz
$ tar -zxf babel-1.0.2.tar.gz 
$ ls
babel-1.0.2/  babel-1.0.2.tar.gz
$ cd babel-1.0.2
$ ls
acinclude.m4  BUGS        configure.ac  examples/  Makefile.am  runtime/
aclocal.m4    CHANGES     contrib/      INSTALL    Makefile.in  share/
ANNOUNCE      compiler/   COPYRIGHT     lib/       README       THANKS
bin/          configure*  doc/          LICENSE    regression/
```

## Step 2:  Configuring Babel ##

### Lonestar Specific ###
Java is required to build Babel.  On lonestar, the you can make java available by doing the following:
```
$ module add java
```

Additionally, on lone star you should be able to omit the FFLAGS and FCFLAGS from the configure:

The following is a single line.
```
$ ./configure --prefix=$PYMCT_ROOT/$BABEL    PYTHON=$PYTHON      CHASMPREFIX=$PYMCT_ROOT
   --with-F90-vendor=Intel      F77=$NON_MPI_F77    FC=$NON_MPI_F90 
```

### General ###

The following is a single line.
```
$ ./configure --prefix=$PYMCT_ROOT/$BABEL    PYTHON=$PYTHON      CHASMPREFIX=$PYMCT_ROOT
   --with-F90-vendor=Absoft      F77=$NON_MPI_F77    FC=$NON_MPI_F90      FFLAGS=$FFLAGS 
    FCFLAGS=$FCFLAGS
```

Configuration takes a while.  When it completes successfully you'll be greeted with the following:
```
configure: creating ./config.status
config.status: creating Makefile
config.status: creating m4/Makefile
config.status: creating bin/Makefile
config.status: creating config/Makefile
config.status: creating java/Makefile
config.status: creating sidl/Makefile
config.status: creating sidlx/Makefile
config.status: creating sidl/babel_internal.h
config.status: sidl/babel_internal.h is unchanged

          Fortran90 enabled.
          Python enabled.
          Fortran77 enabled.
          C++ enabled.
          Java support disabled against request
            (no jvm.dll/libjvm.so/libjvm.a found!)

When reporting issues to babel-bugs@cca-forum.org, please include the
output from babel-config --version-full.
*** Configured BABEL
touch timestamps/babelconfig
```

If Java is disabled you should still be fine.  It is not necessary to use the Java bindings for PyMCT (of course).

## Step 3:  Building and Installing ##

```
$ make
$ make install
```

| [InstallPyMCT](InstallPyMCT.md) | [ServerInstall](ServerInstall.md) |
|:--------------------------------|:----------------------------------|





