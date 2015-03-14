# Introduction #

Babel exports the following MCT Types:
  * Accumulator
  * AttrVect
  * GeneralGrid
  * GlobalSegMap
  * Rearranger
  * Router
  * SparseMatrix
  * SparseMatrixPlus
  * World

Of these types, all but
  * Accumulator
  * GeneralGrid
  * SparseMatrixPlus
are used directly in the python version of the Coupler.  (I say "directly" because they may still be used internally by the other MCT routines)

The most commonly used object in this Python Coupler is:
  * AttrVect

Attribute Vectors are also one of the few MCT types that do not require any MPI environment to be initialized, so they're useful for testing that PyMCT is actually working.

# Why wrap them again? #

The main reason that drove me to wrap the AttrVect type was the idea of an **Interactive** Python Coupler or other MCT application.

Without my wrapping:
  * Any mistyped field names would cause _Abort_ to be called and kill the process.
  * There was no good way to query the names of the fields in an AttrVect
  * There was no way to loop over the fields of an AttrVect (this is necessary in my opinion to enable us to really take advantage of Python).
  * We're choosing to use AttributeVectors **in place of** another CPL6 type called "Bundle", and Bundle's need to support a few coupler specific functions.  By overloading the AttrVect class from babel we can provide those methods.

Because of these problems I have wrapped the AttrVect class exported by Babel in a "shadow class" named AttributeVector.  (Shadow class:  It "shadows" (provides it's own version of) each of the methods of AttrVect).
Here's a quick display of the contents of the two classes:

| AttrVect (babel) | AttributeVector (shadow) |
|:-----------------|:-------------------------|
| av.init( "ifields", "rfields", lsize )                  |  AttributeVector( [ifields](ifields.md), [rfields](rfields.md), lsize ) OR av.initialize( ifields, rfields, lsize )     |
| av.initv( in\_av, size ) | av.initv( in\_av, size ) |
| av.nIAttr() / av.nRAttr() | av.nIAttr() / av.nRAttr() |
| av.lsize() | av.lsize() |
|  | av.countFields() |
|  | av.getFields() |
| av.importIAttr(...) / av.exportIAttr(...) | av.importIAttr(...) / av.exportIAttr(...) |
| av.send(...) | av.send(...) |
|  | av.writeNC( path, ni, nj, n ) |
|  | av.nanCheck( name ) |
|  | av.extremes( name ) |
|  | av.mult( in\_av, field, bfields ) |
|  | av.divide( in\_av, field, bfields ) |

All the above operations have error checking and can throw a simple exception from the error.py module instead of calling Abort.

This allows someone running interactively to notice their mistake early -at the Python level, before any Fortran code is called- and recover from the mistake.  In production code exceptions can be used to display more detailed errors before shutting down.

AttributeVectors supply a number of application specific methods that normal AttrVect objects do not (at least in the current version of MCT).
  * writeNC( ... ):  Simply writes the contents of the real data fields to a NetCDF file
  * nanCheck( ... ):  Debugging / diagnostic function, performs a simple test that should determine if ANY fields (real or integer) contain NAN.  The exceptions (if they occur) are caught and used to simply report **all** the malformed fields.
  * extremes( ... ): Similar to nanCheck, this method finds the min and max of each field in an AttributeVector and displays it.  Useful for a quick sanity check: i.e.:  "OH!  This is probably wrong as the temperature range is between -9999 and 0!"
  * mult( ... ):  Allows you to multiply an entire AttributeVector (or a subset) by a field from another AttributeVector.  Only affects real fields.
  * divide( ... ):  Allows you to divide an entire AttributeVector (or a subset) by a field from another AttributeVector.   Only affects real fields.  Performs division with multiplication by the reciprocal so as to avoid DivisionByZero type errors.

