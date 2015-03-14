| [InstallPyMCT](InstallPyMCT.md) | [MaroonInstall](MaroonInstall.md) |
|:--------------------------------|:----------------------------------|

# PyMCT Python Client Installation #

## Configuring ##

First of all you need to be in the _babel_ directory we created in the PyMCT Server step.

```
$ pwd 
/home/steder/Sandbox/MCT/babel
$ ls
example-c++/     genscl.py  Makefile        mctserver/  README
example-python/  impls/     mctclient-c++/  MCT.sidl
```

From this directory you can regenerate the babel bindings with your specific version of Babel by running the following command:

```
$PYMCT_ROOT/$BABEL/bin/babel -cpython MCT.sidl -o mctclient-python/
```

This should take a few minutes to finish processing.

## Editing the Makefiles ##

Space reserved.

## Building the PyMCT Client Library(ies) ##

There's really a single command to run to build the client libraries.  However, it is a fairly long and complex command.  The only thing that you have to be careful off in the following commands is that you will need to replace python2.4 with the appropriate pythonX.Y major and minor version numbers for the version you are using.  This means that Python 2.5 users will substitude python2.5 for python2.4 in the following.

It may be helpful to define a simple variable in your environment, perhaps called PY\_VERSION that you can use in place of _python2.4_ or _python2.5_ in the commands below.

```
$ pwd 
/home/steder/Sandbox/MCT/babel
$ ls
example-c++/     genscl.py  Makefile        mctserver/  README
example-python/  impls/     mctclient-c++/  mctclient-python/  MCT.sidl
$ cd mctclient-python/
$ CC=mpicc LDSHARED=mpicc LDFLAGS="-shared" $PYTHON setup.py build_ext \
        --include-dirs=$PYMCT_ROOT/include \
    --include-dirs=$PYMCT_ROOT/$BABEL/include \
    --include-dirs=$PYMCT_ROOT/$BABEL/include/$PYVERSION \
    --include-dirs=$PYMCT_ROOT/$BABEL/include/$PYVERSION/babel \
    --include-dirs=$PYMCT_ROOT/include/$PYVERSION \
    --include-dirs=$PYMCT_ROOT/include/$PYVERSION/Numeric \
    --include-dirs=$PYMCT_ROOT/include/$PYVERSION/babel \
    --library-dirs=$PYMCT_ROOT/lib \
    --library-dirs=$PYMCT_ROOT/$BABEL/lib \
    --libraries=sidlstub_f90 build --force
```

And to install the module:
```
$ $PYTHON setup.py install 
```

## Generate the SCL File ##

```
$ pwd
/home/steder/Sandbox/MCT/babel/mctclient-python
$ cd ../
$ pwd
/home/steder/Sandbox/MCT/babel
$ ls
CVS             genscl.py  mctclient-c++           mctserver        README
example-c++     impls      mctclient-python        mctserver-1.0.2
example-python  Makefile   mctclient-python-1.0.2  MCT.sidl
$ $PYTHON genscl.py $PYMCT_ROOT/lib/libmctserver.so
$ cat $PYMCT_ROOT/lib/libmctserver.so.scl 
<?xml version="1.0" ?>
        <scl>
                <library uri="libmctserver.so"
                scope="local" resolution="lazy" >
                        <class name="MCT.Accumulator" desc="ior/impl"/>
                        <class name="MCT.AttrVect" desc="ior/impl"/>
                        <class name="MCT.GeneralGrid" desc="ior/impl"/>
                        <class name="MCT.GlobalSegMap" desc="ior/impl"/>
                        <class name="MCT.Rearranger" desc="ior/impl"/>
                        <class name="MCT.Router" desc="ior/impl"/>
                        <class name="MCT.SparseMatrixPlus" desc="ior/impl"/>
                        <class name="MCT.SparseMatrix" desc="ior/impl"/>
                        <class name="MCT.World" desc="ior/impl"/>
                </library>
        </scl>
$ 
```

## Quick Test ##

First, start python.
```
$ $PYTHON
```

Now import the MCT module and try instantiating a basic object.
```
>>> import MCT
>>> from MCT import AttrVect
Loading main: ok
Searching for class MCT.AttrVect, target ior/impl, file /home/steder/pyMCT//lib/
>>> av = AttrVect.AttrVect()
>>> av.init("hello:world","from:mct",10)
```

| [InstallPyMCT](InstallPyMCT.md) | [MaroonInstall](MaroonInstall.md) |
|:--------------------------------|:----------------------------------|


