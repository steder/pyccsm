| [InstallPyMCT](InstallPyMCT.md) | [NumpyInstall](NumpyInstall.md) |
|:--------------------------------|:--------------------------------|

# Chasm Installation #

## About Chasm ##
Blatently stolen from http://chasm-interop.sourceforge.net/:

Chasm is a language transformation system created to improve Fortran and C/C++ language interoperability. Chasm (with the aid of associated tools) parses Fortran and C/C++ source code and automatically generates bridging code that can be used to make calls to routines in the foreign language. Chasm also supplies an array-descriptor library that provides an interface between C and F90 assumed-shape arrays. This allows arrays to be created in one language and then passed to and used by the other (foreign) language.

## Step 1:  Obtaining the Chasm source ##

Chasm can be downloaded directly from http://sourceforge.net/projects/chasm-interop/

Download version 1.3 or later.

Once you have Chasm downloaded to a directory on your system (let's assume it's  /home/user/Sandbox/ again) you can extra the source like this:

```
$ pwd
/home/steder/Sandbox
$ ls
chasm_1.4.RC1.tar.gz
$ tar -zxf chasm_1.4.RC1.tar.gz 
$ ls     
chasm/  chasm_1.4.RC1.tar.gz
$ cd chasm
$ ls
aclocal.m4      configure*    doc/      INSTALL   Makefile.in  test/
bin/            configure.in  example/  interop/  README       xform/
configuration/  contrib/      include/  LICENSE   src/         xml/
```

## Step 2:  Configuring Chasm ##

### General ###
```
$ ./configure --prefix=$PYMCT_ROOT --with-F90=$F90 --with-F90-vendor=$CHASM_FORTRAN_VENDOR
```

### Lonestar ###
```
$ ./configure --prefix=$PYMCT_ROOT --with-F90=$F90 --with-F90-vendor=$CHASM_FORTRAN_VENDOR CFLAGS="-fPIC" F90FLAGS="-fPIC"
```

I think that StageOneInstallPyMCT may need to be revamped to include -fPIC (for creating position independent code, useful for dynamic linking) in the FCFLAGS and F90FLAGS environment variables, at least on lonestar.

## Step 3:  Building and Installing ##

```
$ make
$ make install
```

| [InstallPyMCT](InstallPyMCT.md) | [NumpyInstall](NumpyInstall.md) |
|:--------------------------------|:--------------------------------|



