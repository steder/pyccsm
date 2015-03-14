# Introduction #

First download 2.4+ series python either tarball or bzip2 source distributions and download it to your Sandbox:

For example:
```
% ls
babel-1.0.2.tar.gz  babel.tar.gz  MCT  MCTnodata.tar.gz  Python-2.5.tar.bz2
% tar -jxf Python-2.5.tar.bz2 
```

# Configuration #

```
  $ ./configure --prefix=$PYMCT_ROOT --enable-shared 
     --enable-unicode=ucs4 --with-threads --with-universal-newlines --with-pymalloc 
     --with-cxx=g++
```

# Build and Install #

```
$ make
$ make install
```