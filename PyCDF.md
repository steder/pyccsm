# Installing PyCDF #

## Obtaining PyCDF ##

First you'll need to visit [http://pysclint.sourceforge.net/pycdf/](http://pysclint.sourceforge.net/pycdf/) and download the PyCDF source.

Assuming you've downloaded the PyCDF source to a convenient directory:
```
$ pwd      
/home/steder/Sandbox
$ ls
pycdf-0.6-2-rc1.tar.gz
$ tar -zxf pycdf-0.6-2-rc1.tar.gz 
$ ls
pycdf-0.6-2-rc1/  pycdf-0.6-2-rc1.tar.gz
$ cd pycdf-0.6-2-rc1
$ ls
CHANGES  INSTALL  PKG-INFO  README  doc/  examples/  pycdf/  setup.cfg  setup.py
```

## Configuring PyCDF ##

You will most likely need to modify the setup.py script in the pycdf directory to point to your specific installation of Numeric/Numpy and NetCDF.

Look for the variable _USE=_ and set it to NUMERIC.
```
USE = NUMERIC         # use the Numeric package (default)                                       
#USE = NUMARRAY         # use the numarray package       
```

Then look for:
```
if USE == NUMERIC:
    _pycdf_ext = Extension('pycdf._pycdfext',
                           sources   = ["pycdf/pycdfext/numeric/pycdfext_wrap.c"],
                           #library_dirs=["non standard path where libs live"],                 
                           libraries = ["netcdf"])
    shutil.copy("pycdf/pycdfext/numeric/pycdfext.py", "pycdf")
    shutil.copy("pycdf/pycdfext/numeric/pycdfext_array.py", "pycdf")
```

Modify the above "library\_dirs" variable to point at your NetCDF libraries, wherever they happen to be installed.

For instance:

```
if USE == NUMERIC:
    _pycdf_ext = Extension('pycdf._pycdfext',
                           sources   = ["pycdf/pycdfext/numeric/pycdfext_wrap.c"],
                           library_dirs=["/home/steder/software/lib/","/opt/netcdf/lib/",],                 
                           libraries = ["netcdf"])
    shutil.copy("pycdf/pycdfext/numeric/pycdfext.py", "pycdf")
    shutil.copy("pycdf/pycdfext/numeric/pycdfext_array.py", "pycdf")
```

## Building PyCDF ##
Building PyCDF is quite simple, just run the setup.py script like this:
```
$ $PYTHON setup.py build
```

## Installing PyCDF ##

```
$ $PYTHON setup.py install
```
