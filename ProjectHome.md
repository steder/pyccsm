# PyMCT Start Page #

## Introduction ##

PyMCT consists of a suite of software packages necessary to build and run a Python Coupler like PyCCSM.
MCT is a high performance regridding and parallel communication package designed to address issues of coupling multiple scientific models on different scales and grids to one another.

# Architecture #

## System Overview ##

PyMCT wraps the Fortran 90 MCT library and makes it available to Python Programmers.  This wrapping is not normally possible due to complex issues with the many different interpretations of the Fortran 90 standard.  However, a powrful tool, the CCA Babel project allowed us to address these issues and completely wrap the MCT functionality into Python.

## Components ##

The following components are necessary to build and use MCT and PyMCT

### MCT ###
  * Fortran 90 Compiler
  * MPI (we use MPICH2)

### Babel ###
  * C, Java, and Fortran compilers ( Portland Group Fortran is NOT supported! <email them!> )
  * Python Interpreter
  * Chasm
  * Numpy or Numeric Python array packages

### Python MPI ###

There are a number of Python bindings to MPI out there, but the only one that completely meets the requirements of our package is our own solution, MaroonMPI.  Depending on your application you may have reasonable success with other bindings.

### NetCDF ###

# Configuration and Installation #

Please see our [InstallPyMCT](InstallPyMCT.md) page.