| [PyCCSM](PyCCSM.md) |
|:--------------------|
| [InstallPyCCSM](InstallPyCCSM.md) |

# Introduction #

# Overview #

## Python files ##

```
main.py
coupler.py
flux.py
frac.py
cpl/
cpl/__init__.py
cpl/attributevector.py  - Python wrapper class for MCT.AttrVect object exported by Babel.
cpl/bundle.py - Defunct "bundle" object.
cpl/comm.py - Sets up communication at run time automatically.  Essentially, a simple MPH replacement tailored to our appplication.
cpl/const.py - Defines a series of physical constants, many simply refer to the values defined in cpl/shr/const.py
cpl/contract.py - Defines the main coupler setup routines and intercommunication
cpl/control.py - Defines a series of coupler control flags and a function for reading a namelist file
cpl/debug.py - Defines a simple replacement for "print <msg>". i.e.:  "debugPrint(<msg>)".
cpl/decomposition.py - Defines the primary decomposition scheme used by domain objects along with placeholder functions for several other possible schemes.
cpl/domain.py - A python wrapper class for the MCT.GlobalSegMap type exposed by Babel.
cpl/error.py - Defines a number of coupler exceptions for the various classes. (AttributeVectorError, BundleError, etc).  They provide a base to create customized exceptions for each class.
cpl/fields.py - Defines a number of lists and dictionaries defining all the fields into and out of each component.  This is just a simple translation of a fortran file.
cpl/infobuffer.py - A complete reimplementation of the infobuffer class, which is not an MCT type, but a CPL6 type.
cpl/interface.py - Placeholder for simple interface functions to be used by coupler.py.  This file corresponds not to a MCT file but to a CPL6 file.
cpl/map.py - A conversion of the CPL6 mapping routines.  Uses MCT.SparseMatrix and Rearranger.py.
cpl/rearranger.py - Wrapper class around MCT.Rearranger which is exported by Babel.
cpl/timer.py - Pure python timer class.  Used to do basic timings in the python coupler.
cpl/shr/__init__.py
cpl/shr/cal.py - Reimplementation of simple F90 file cpl_shr_cal.F90.  This is a simplified julian calendar.
cpl/shr/const.py - Simply defines a number of constants
cpl/shr/date.py - A basic date class taht corresponds to cpl_shr_date.F90 (reimplemenation of a CPL6 file)
cpl/shr/orb.py - Basic orbital physics functions (Reimplemented from CPL6)
```


# In-depth #

## Reimplemented from CPL6 ##

### Basic Dependencies ###

  * cpl/shr/cal.py
  * cpl/shr/date.py
  * cpl/shr/orb.py
  * cpl/shr/const.py

### Coupler ###

  * cpl/contract.py
  * cpl/comm.py
  * cpl/domain.py
  * cpl/map.py
  * cpl/fields.py
  * cpl/const.py
  * cpl/infobuffer.py
  * cpl/control.py

  * main.py
  * frac.py
  * flux.py

## New Code ##

  * cpl/attributevector.py
  * cpl/decomposition.py
  * cpl/timer.py
  * cpl/debug.py
  * cpl/error.py
  * coupler.py

NewCouplerCode

## Dropped Code ##

  * cpl/bundle.py -> Replaced by an attributevector
  * merge.py (merge\_mod.F90) -> Rolled into coupler.py

WhyRemoveMergeMod

## Wrapped Babel+MCT types ##

  * cpl/attributevector.py -> AttrVect
  * cpl/rearranger.py -> Rearranger

For some justification of why these layers were deemed necessary take a look at WhyWrapBabelTypes