"""
attributevector.py

Defines baseclase AttrVect, which is inherited by Bundle and Domain
to define how those classes access their internal data.  

This internal data is stored in MCT Fortran 90 derived types known
as Attribute Vectors which define a fairly simple interface for 
getting and setting fields, and for sending and receiving data 
to and from the coupler.

The following interface is defined in this module for all
inherited types:

---
void    self.initialize( ifields="", rfields="", lsize = 10 )# actually
initializes the fields and datasize of the AttrVect after creating the object.
(In a lot of cases we don't know this information until a while after we 
actually create the object.  This is a short term solution to make it
easier to translate code. )

void    self.send( Router r )
Array1D self.recv( Router r )
void    self.__setitem__( String f, Array1D data )
Array1D self.__getitem__( String f )
Integer self.__len__( void )
void    self.addIntegerField( String f )
Boolean self.safe( void ) # returns true if initialize has been called

Listing of all methods in AttrVect that have NOT been implemented:
['AllReduce',
 'CopyI',
 'CopyR',
 'GlobalReduce',
 'GlobalWeightedSumRAttr',
 'LocalReduce',
 'LocalReduceRAttr',
 'LocalWeightedSumRAttr',
 'MaskedSpatAverage',
 'MaskedSpatAverageV',
 'MaskedSpatIntegral',
 'MaskedSpatIntegralV',
 'Permute',
 'SharedAttrIndex',
 'Sort',
 'SortPermute',
 'SpatAverage',
 'SpatAverageV',
 'SpatIntegral',
 'SpatIntegralV',
 'appendIAttr',
 'appendRAttr',
 'bcast',
 'clean', -> __del__
 'deleteRef',
 'initIList',
 'initRList',
 'irecv',
 'isend',
 'rearrange',
 'sMatAvMult_sMPlus',
 'waitrecv',
 'waitsend',
 'zero']

---
"""

# Builtin Python Modules (sys, os, string, math, etc)
import os

#
import pycdf

# External Modules (Numeric, mpi, etc)
import Numeric
from MCT import AttrVect

# Internal Modules
import error
AttributeVectorError = error.AttributeVectorError

from debug import newDebugPrint
debugPrint = newDebugPrint(False)


class AttributeVector:
    """
    Input Arguments:
    ifields - List of integer fields defined in this vector
    rfields - List of real fields defined in this vector
    lsize   - Integer describing the size of all the defined fields

    returns None
    """
    def __init__( self, ifields=[], rfields=[], lsize=10 ):
        self._av = AttrVect.AttrVect()
        self._ifields = {}
        self._rfields = {}
        self._lsize = -1
        self._count = 0
        self.initialize( ifields, rfields, lsize )
        return

    def count(self):
        return self._count

    def increment(self,inc=1):
        self._count += inc
        return self._count

    def reset(self,s=0):
        self._count = s
        return self._count
    
    def __repr__(self):
        s = "AttributeVector - Field Size: %d\n"%(self._lsize)
        i = 0
        for f in self._ifields.keys():
            s += "  Integer Field %d: %s\n" %(i, f)
            i += 1
        i = 0
        for f in self._rfields.keys():
            s += "  Real Field %d: %s\n" %(i,f)
            i += 1
        return s
    
    def __str__(self):
        return self.__repr__()
    
    def __del__(self):
        """
        Deallocates the memory of the internal MCT datatype AttrVect.
        """
        self._av.clean()

    def initialize( self, ifields, rfields, lsize ):
        """
        Used to initialize this attribute vector after the contructor
        is called.

        Required Arguments:
        ifields - List of integers field names (a list of strings)
        rfields - List of real field names( a list of strings )
        lsize   - Integer size to be allocated for all fields

        Returns no values
        """
        if not ifields:
            ifields = []
        else:
            i = 0
            for field in ifields:
                self._ifields[field] = i
                i+=1
        
        if not rfields:
            rfields = []
        else:
            i = 0
            for field in rfields:
                self._rfields[field] = i
                i+=1

        if lsize < 0:
            raise AttributeVectorError,"Attribute vectors cannot have negative sizes"
        self._lsize = lsize
        debugPrint("Before av.init:")
        self._av.init( ":".join(ifields), ":".join(rfields), lsize )
        debugPrint("After av.init:")
        return

    def initv( self, in_av, size ):
        if( isinstance( in_av, AttributeVector ) ):
            self._ifields,self._rfields = in_av._ifields, in_av._rfields
            retval = self._av.initv( in_av._av, size )
        else:
            retval = self._av.initv( in_av, size )
        self._lsize = size
        return retval 

    def nIAttr(self):
        return self._av.nIAttr()

    def nRAttr(self):
        return self._av.nRAttr()
    
    def __len__( self ):
        """
        This function returns the total number of fields inside
        this attribute vector.
        """
        a,b = self.countFields()
        return (a+b)
        
    def size( self ):
        """
        This method returns the length of a single field 
        of this vector.
        """
        return self._av.lsize()

    def lsize(self):
        return self._av.lsize()

    def copy( self, inav ):
        """
        Copies and overwrites the values of this AttributeVector with
        the AttributeVector 'av' that is passed in.

        newcopy.copy( oldAV ) # newcopy is now a duplicate of oldAV.
        """
        if( isinstance(inav, AttributeVector) ):
            self._av.Copy( inav._av )
        else:
            self._av.Copy( inav )
        return

    def indexIA( self, fieldName ):
        print "Warning:  indexIA returns Fortran (1-based) indices!"
        return ( self._av.indexIA( fieldName ) )

    def indexRA( self, fieldName ):
        print "Warning:  indexRA returns Fortran (1-based) indices!"
        return ( self._av.indexRA( fieldName ) )
        
    def countFields(self ):
        """
        This method returns a tuple containing the 
        number of Integer and Real fields contained
        in the attribute vector or bundle.
        
        num_ifields, num_rfields = mybundle.countFields()
        """
        return ( self.nIAttr(), self.nRAttr() )
    
    def getFields( self ):
        """
        Returns a tuple ( integer_fields_list, real_fields_list )
        """
        #raise NotImplementedError
        return self._ifields, self._rfields
    
    def importIAttr( self, field, data ):
        """
        """
        if field in self._ifields.keys():
            if ( len(data) == self.size() ):
                self._av.importIAttr( field, data )
            else:
                raise AttributeVectorError,"Fields should be size: %d, Data passed in was %d"%(self.size(), len(data))
        else:
            raise AttributeVectorError,"Tried to set a value in a non-existant field:%s"%(field)
    
    def exportIAttr( self, field ):
        """
        """
        if field in self._ifields.keys():
            data,length = self._av.exportIAttr(field)
            return data,length
        else:
            raise AttributeVectorError,"Integer Attribute %s is not defined"%(field)
    
    def importRAttr( self, field, data ):
        """
        """
        debugPrint("field:",field,"\nrfields:",self._rfields)
        if field in self._rfields.keys():
            if ( len(data) == self.size() ):
                self._av.importRAttr( field, data )
            else:
                raise AttributeVectorError,"Fields should be size: %d, Data passed in was %d"%(self.size(), len(data))
        else:
            raise AttributeVectorError,"Tried to set a value in a non-existant field:%s"%(field)
            
    def exportRAttr( self, field ):
        """
        """
        if field in self._rfields.keys():
            data,length = self._av.exportRAttr(field) 
            return data,length
        else:
            raise AttributeVectorError,"Real Attribute %s is not defined" % (field)
    
    def send( self, router, tag=600 ):
        """
        This method is a thin wrapper around AttrVect.send
        so that anything that inherits from this class( bundles and
        domains ) will define a meaningful send method

        Required Arguments:
        router - MCT Router Object
        tag    - Integer similar to an MPI tag value

        Returns:
        router
        """
        
        debugPrint("Beginning AV Send...")
        retval = self._av.send( router, tag )
        debugPrint("Finished AV Send...")
        return retval

    def recv( self, router, tag=600, sum=False ):
        """
        This method is a thin wrapper around AttrVect.recv
        so that anything that inherits from this class (bundles nad
        domains ) will have a meaningful receive method.

        Required Arguments:
        router - MCT Router Object
        tag    - Integer similar to an MPI Tag value

        Optional Arguments:
        sum    - Boolean, not sure what it defines at the moment

        Returns:
        router
        """
        debugPrint("Beginning AV Recv...")
        self._av.recv( router,tag,sum )
        debugPrint("Finished AV Recv")
        return router

    def scatter( self, avin, gsmap, root, myid, comm ):
        """
        scatter( avin, gsmap, root, myid, comm )

        Scatter attribute vector 'avin' on communicator 'comm'.
        'gsmap' is a GlobalSegMap that defines
        how to divide up the attribute vector.
        Root is the root processor managing the scatter,
        and myid is the rank of the local processor taking
        part in the scatter.

    
        """
        if( isinstance( avin, AttributeVector ) ):
            self._av.scatter( avin._av, gsmap, root, myid, comm )
        else:
            self._av.scatter( avin, gsmap, root, myid, comm )
        return
    
    def gather( self, avin, gsmap, root, myid, comm ):
        """
        scatter( avin, gsmap, root, myid, comm )
        
        Gather attribute vector avin on communicator 'comm'.
        'gsmap' is a GlobalSegMap that defines
        how to divide up the attribute vector.
        Root is the root processor managing the scatter,
        and myid is the rank of the local processor taking
        part in the scatter.
        """
        debugPrint( "Entering gather..." )
        debugPrint("avin._av, gsmap, root, myid, comm=",avin._av, gsmap, root, myid, comm)
        if( isinstance(avin, AttributeVector) ):
            retval = self._av.gather( avin._av, gsmap, root, myid, comm )
        else:
            retval = self._av.gather( avin, gsmap, root, myid, comm )
        debugPrint("Leaving gather...")
        return retval

    def bcast( self, root, myid, comm ):
        retval = self._av.bcast( root, myid, comm )
        return retval
    
    def sMatAvMult_DataLocal( self, inav, sMat ):
        debugPrint("Entering sMatAvMult_DataLocal...")
        if(isinstance(inav, AttributeVector)):
            val = self._av.sMatAvMult_DataLocal( inav._av, sMat )
        else:
            val = self._av.sMatAvMult_DataLocal( inav, sMat )
        debugPrint("Leaving sMatAvMult_DataLocal...")
        return val

    def zero( self ):
        """
        Set all fields to 0s.
        """
        self._av.zero()
        return

    def writeNC( self, path, ni, nj, n ):
        """
        Writes all the fields in the Attribute Vector
        to a NetCDF file.
        """
        print path
        if(os.path.exists(path)):
            os.remove(path)
        ncfile = pycdf.CDF(path,(pycdf.NC.WRITE|pycdf.NC.CREATE))
        ncfile.definemode()
        ncfile.def_dim("time",pycdf.NC.UNLIMITED)
        ncfile.def_dim("nj", nj)
        ncfile.def_dim("ni", ni)
        ncfile.def_dim("n", n)
        ncfile.datamode()
        ncfile.sync()
        ifields,rfields = self.getFields()
        for r in rfields:
            debugPrint( "* outputing %s field to %s"%(r,path) )
            ncfile.definemode()
            ncfile.def_var(r,pycdf.NC.DOUBLE,("time","nj","ni"))
            ncfile.datamode()
            ncfile.sync()
            tmpdata,tmpsize = self.exportRAttr(r)
            
            # ncfile.var("ice_temp").put(ice_temperature)#[,start=(time,0)])
            data = Numeric.resize(tmpdata,(nj,ni))
            # statusPrint( "TIMESTEP:",TIMESTEP )
            ncfile.var(r).put( data, start =(0,0,0) )
            ncfile.sync()
        ncfile.close()
        # statusPrint( "NetCDF file written!" )
        # Finished outputting fields to NetCDF
        return

    def fakewriteNC( self, path, ni, nj, n ):
        print "fakewriteNC: %s,%s,%s,%s"%(path,ni,nj,n)
    
    def nanCheck( self, name="Unknown" ):
        tmpi = Numeric.ones((self.size()),Numeric.Int32)
        for i in self._ifields.keys():
            try:
                ifield,isize = self.exportIAttr( i )
                ifield *= tmpi
            except OverflowError:
                print "NAN(s) present in bundle %s field %s"%(name,i)
        tmpr = Numeric.ones((self.size()),Numeric.Float64)
        for r in self._rfields.keys():
            try:
                rfield,rsize = self.exportRAttr( r )
                rfield *= tmpr
            except OverflowError:
                print "NAN(s) present in bundle %s field %s"%(name,r)
        return

    def extremes( self, name="Unknown" ):
        for i in self._ifields.keys():
            ifield,isize = self.exportIAttr(i)
            mini,maxi = 0,ifield[0]
            for n,e in enumerate(ifield):
                if ( e > maxi ):
                    maxi = e
                if ( e < mini ):
                    mini = e
            print "%s bundle - integer field %s : min(%s) - max(%s)"%(name,i,mini,maxi)
        for r in self._rfields.keys():
            rfield,isize = self.exportRAttr(r)
            minr,maxr = 0,rfield[0]
            for n,e in enumerate(rfield):
                if ( e > maxr ):
                    maxr = e
                if ( e < minr ):
                    minr = e
            print "%s bundle - real field %s : min(%s) - max(%s)"%(name,r,minr,maxr)
            
    def mult( self, inav, field, bfields=None ):
        """
        """  
        if( bfields ):
            bunfields = bfields
        else:
            bunfields = self._rfields.keys()
        try:
            factor,size = inav.exportRAttr(field)
            for rfield in bunfields:
                tmp,size = self.exportRAttr(rfield)
                tmp *= factor
                self.importRAttr( rfield, tmp )
        except OverflowError:
            print "Overflow exception occured in %s(size of field:%s)"%(rfield,size)
            #for i in xrange(0,int(size*0.10)):
            #    print "%s[%s]*factor[%s] = %s*%s"%(rfield,i,i,tmp[i],factor[i])
            #print recip
            raise
        except:
            raise
        return

    def divide(self, inav, field, bfields=None ):
        """
        Calculates the product of self * the reciprocal of field from inav 
        """
        
        if( bfields ):
            bunfields = bfields
        else:
            bunfields = self._rfields.keys()
        try:
            recip,size = inav.exportRAttr(field)
            for i,v in enumerate(recip):
                if ( recip[i] != 0 ):
                    recip[i] = 1.0 / v
                else:
                    recip[i] = 0.0
            for rfield in bunfields:
                tmp,size = self.exportRAttr(rfield)
                tmp *= recip
                self.importRAttr( rfield, tmp )
        except OverflowError:
            print "Overflow exception occured in %s(size of field:%s)"%(rfield,size)
            #for i in xrange(0,int(size*0.10)):
            #    print "numerator[%s] / denominator[%s] = %s/%s"%(i,i,tmp[i],recip[i])
            #print recip
            raise
        except:
            raise
        return
    
if __name__=="__main__":
    # Testing creation and deletion
    print "Creating an attribute vector and providing the following values:"
    print "ifields = ['i1','i2'], rfields=['r1','r2'], and lsize = 50"
    av = AttributeVector( ["i1","i2"],["r1","r2"], 50 )
    print "Real AV Size: %d" %(av.size())
    del av
    
