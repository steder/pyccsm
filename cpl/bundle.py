"""
Bundles are simple objects that store data and manage the translation of 
data between the fortran and python layers.

Ultimately bundles are just simple wrappers around the fortran Attribute 
Vector structure.  

bundle.py:
---

Bundle Object now built on AttributeVector objects: 2/28/05 - Mike Steder
Bundle Skeleton:  1/5/05 - Mike Steder
"""

# Builtin Python Modules
import traceback

# External Modules
from MCT import AttrVect

# Internal Modules
import attributevector
import domain
import error

BundleError = error.BundleError

class Bundle(attributevector.AttributeVector):
    def __init__(self, ifields, rfields, lsize, dom, name="unnamed", count=0 ):
        """
        idata and rdata parameters should be 2D numeric arrays ar[i][j], where
        i indices correspond to specific data fields and j indices correspond 
        to space.
        """        
        attributevector.AttributeVector.__init__(self, ifields, rfields, lsize)
        self.name = name
        self.domain  = dom
        self.count   = count
        return

    def __del__(self):
        # MCT Types need to be cleaned up
        self.av.init("","",10)
        self.av.clean()
    
    def info(self):
        """
        Bundle.info()

        Prints information about the bundle object
        """
        print "bundle name = %s" % (self.name),
        print "domain name = %s" % (self.domain.name),
        print "accumulation count = %s" % (self.count)
        try:
            self.av.info(1)
        except:
            traceback.print_exc()
        return

    def fill(self):
        """
        Fill the bundle with a test case:
        
        mybundle = Bundle()
        mybundle.fill()

        Takes no arguments, modifies its bundle in place.
        """
        nflds = self.av.nRAttr()
        npts = self.av.lsize()
        # these two loops are reversed in order:
        # because the ordering of arrays in Numeric/C 
        # is the opposite of fortrans...
        # RIGHT?
        for i in xrange(npts):
            for j in xrange(nflds):
                self.idata[i][j] = ( (j * 10.0) + sin( i / 500.0 ) )
                self.rdata[i][j] = ( (j * 10.0) + sin( i / 500.0 ) )
        self.count = 1
        return

    def dump(self, file=None):
        """
        ** Not Implemented **
        
        Write bundle to file, standard out, whatever.

        For debugging purposes, prints contents of a bundle to an output file
        
        # Takes one argument:
        iun : base unit number
        """
        if file == None:
            pass
        else:
            raise NotImplementedError,"Bundle.dump Not Implemented"
        avnum = self.av.nRAttr()
        avsize = self.av.lsize()
        id1 = self.domain.lgrid.indexRA("index")

        print "Bundle.dump: info:",avnum,avsize,id1
        
        return

    # Copy and Fcopy should be just overloaded =?
    def copy(self,other,bunrlist=None, bunilist=None, buntrlist=None, buntilist=None):
        """
        Copy data from one bundle to another
        
        mybundle.copy( other, [bunrlist, bunilist[, buntrlist, buntilist]] )
        
        This routine copies from from the input argument other into 
        mybundle the data of all the attributes shared between the two.

        If only a subset of the shared attributes should be copied you
        can specify them in the optional arguments 'bunrlist' and 'bunilist'.

        If any attributes of 'mybundle' have different names for the same
        quantity in 'other' you can provide the optional translation
        list as 'buntrlist' and 'buntilist' which describes which field
        names correspond to which.

        bunilist -> bundle integer field list
        buntilist -> bundle translation list for integer fields
        bunrlist && buntrlist -> simply correspond to the real data.
        
        """
        # If either bundle has no attributes - return:
        if ( ( not self.hasAttr() ) or ( not other.hasAttr() ) ):
            return
        # 
        if ( other.count != 1 ):
            print "WARNING: bundle",other.name,"has accum count =",other.count

        # Copy real attributes if specified:
        if( bunrlist != None ):
            if( buntrlist != none ):
                self.av.CopyR( other.av, bunrlist, buntrlist )
            else:
                self.av.CopyR( other.av, bunrlist )
        # Copy integer attributes if specified:
        if( bunilist != None ):
            if( buntilist != None ):
                self.av.CopyI( other.av, bunilist, buntilist )
            else:
                self.av.CopyI( other.av, bunilist )
        # Otherwise copy everything:
        if( bunilist == None and bunrlist == None ):
            self.av.Copy( other.av )
        self.count = other.count
        
        # Hack:  For now just bruteforce update the python 
        # copy of the AttrVect object's data
        for f in self.ifields.keys():
            self.idata[self.ifields[f]] = self.av.exportIAttr(f)
        for f in self.rfields.keys():
            self.rdata[self.rfields[f]] = self.av.exportRAttr(f)
            
        return

    def fcopy(self, other, bunrlist=None, bunriist=None, buntrlist=None, buntilist=None):
        self.copy( other, bunrlist, bunilist, buntrlist, buntilist )
        """
        Alias to Bundle.copy, see Bundle.copy for usage.
        """
        return

    def split(self):
        """
        Bundle.split(

        This is silly...  This method simply makes copies
        of a bundle.  We can just make them directly using the
        copy method above.
        """
        raise NotImplemented
        return

    def gather(self, *bundles):
        """
        Bundle.gather( bun1[, bun2[, bun3, ...]] )
        
        mybundle.gather( other1, other2 )
        
        This method takes a variable number of bundle objects and
        copys all the data they have in common with the mybundle.
        """
        for bun in bundles:
            self.data.Copy( bun )
        return

    def hasAttr(self):
        """
        if ( mybundle.hasAttr() ):
            print "My bundle has some real or integer attributes!"
        else:
            print "My bundle is empty!"
        
        Return true if input bundle has any real of integer attributes
        """
        if( self.av.nIAttr() > 0 or self.av.nRAttr() > 0 ):
            return True
        else:
            return False

    # These functions, from 'zero' to 'gsum' should probably be
    # implemented inside of MCT.  They only affect the attribute
    # vectors and they will be much slower to operate on in Python..
    def zero(self):
        """
        """
        return

    def accum(self):
        """
        """
        return

    def avg(self):
        """
        """
        return

    def add(self):
        """
        """
        return

    def mult(self, other, field, fieldlist=None, bundle = None):
        """
        product = mybundle.mult( other, field[,fieldlist[,bundle]] )

        This method does a product of a single field in a bundle with
        all fields in another bundle

        Takes the following parameters:
        other # The bundle to multiply by
        field # The field in 'other' to multiply all fields in self by
        fieldlist # Optional sublist of fields to multiply
               # Rather then multiplying all fields you can pick 
               # and choose by specifying a list of fields  here.
        bundle # Optional initialization bundle:
        """
        return

    def divide(self):
        """
        """
        return

    def gsum(self):
        """
        """
        return
