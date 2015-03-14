| [InstallPyMCT](InstallPyMCT.md) | [BabelInstall](BabelInstall.md) |
|:--------------------------------|:--------------------------------|

# Numpy Installation #

## About Numpy ##
Blatently stolen from http://numpy.scipy.org/

NumPy

The fundamental package needed for scientific computing with Python is called NumPy.  This package contains:
  * a powerful N-dimensional array object
  * sophisticated (broadcasting) functions
  * basic linear algebra functions
  * basic Fourier transforms
  * sophisticated random number capabilities
  * tools for integrating Fortran code.
  * NumPy derives from the old Numeric code base and can be used as a replacement for Numeric.   It also adds the features introduced by Numarray and can also be used to replace Numarray.

Numeric users should find the transition relatively easy (although not without some effort).  There is a module (numpy.oldnumeric.alter\_code1) that can make most of the necessary changes to your Python code that used Numeric to work with NumPy's Numeric compatibility module.

## Step 1a:  Obtaining the Numpy source ##

You can download Numpy from the above site.  Here's a direct link to their SourceForge project [page](http://sourceforge.net/projects/numpy/)  Please download the most recent version.

Assuming you download the package to our Sandbox directory:
```
$ pwd
/home/steder/Sandbox
$ ls
numpy-1.0.1.tar.gz
$ tar -zxf numpy-1.0.1.tar.gz 
$ ls
numpy-1.0.1/  numpy-1.0.1.tar.gz
$ cd numpy-1.0.1
$ ls
COMPATIBILITY   MANIFEST.in  README.txt           site.cfg
DEV_README.txt  numpy/       scipy_compatibility  THANKS.txt
LICENSE.txt     PKG-INFO     setup.py*
```

## Step 1b:  Obtaining the Numeric source ##

You can also download the Numeric source from [the Numpy Site](http://numpy.sf.net/).  You should download the most recent version, 24.X+.

Assuming you download the package to our Sandbox directory:
```
lslogin1% pwd
/san/hpc/A-ig2/steder/Sandbox
lslogin1% ls
babel-1.0.2         chasm                 MCTnodata.tar.gz        numpy-1.0.1         Python-2.5.tgz
babel-1.0.2.tar.gz  chasm_1.4.RC1.tar.gz  numpy-1.0.1.tar.gz
babel.tar.gz        MCT                   Numeric-24.2.tar.gz     Python-2.5
lslogin1% tar -zxf Numeric-24.2.tar.gz 
cd Numeric-2lslogin1% 
lslogin1% cd Numeric-24.2
lslogin1% ls
changes.txt   DEVELOPERS  Legal.htm     makedist.bat  Misc      README      setup.py
customize.py  Include     Lib           makedist.sh   Packages  RPM.README  Src
Demo          INSTALL     makeclean.sh  MANIFEST      PKG-INFO  setup.cfg   Test
```

## Step 2a:  Building  and Installing Numpy ##

```
$ $PYTHON setup.py install
```

## Step 2a: Building and Installing Numpy on Lonestar ##

Edit site.cfg in the numpy directory to point to BLAS and LAPACK on your system.

```
lslogin1% more site.cfg
[atlas]
library_dirs = /opt/apps/gotoblas/gotoblas-1.02
atlas_libs = goto
```

You should then be able to proceed with the build normally:
```
$ $PYTHON setup.py install
```

## Step 2b: Building and Installing Numeric ##

```
$ $PYTHON setup.py install
```

| [InstallPyMCT](InstallPyMCT.md) | [BabelInstall](BabelInstall.md) |
|:--------------------------------|:--------------------------------|




