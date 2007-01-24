"""
domain.py

Python version of cpl_domain_mod.F90

Defines the 'Domain' class

Domain objects are a fundamental coupler data type.  A domain
has both a grid and a decomposition.  A decomposition is
described by a global seg map.
"""
# Standard Imports
import sys

# External dependencies:
import Numeric
import mpi

# Coupler Imports
from MCT import GlobalSegMap
import attributevector
import comm
from control import infodbug
from error import DomainError
import fields
debugLevel = infodbug

import debug
debugPrint = debug.newDebugPrint(True)

import Numeric

class Domain:
    """
    Domain objects are a fundamental coupler data type.  A domain
    has both a grid(self.lgrid) and a decomposition.  The decomposition is
    described by a global seg map(gsmap).
    """
    def __init__(self,name=None,suffix=None,n=None,ni=None,nj=None):
        """
        * Initializes parameters:
          name:   Name of domain (eg. "ocean")
          suffix: netCDF domain suffic (eg. "o")
          n:      n = ni * nj - total number of grid points (global)
          ni:     Number of 2d array i indices (global)
          nj:     Number of 2d array j indices (global)

        * Note that self.lgrid and gsmap are not initialized,
        just declared.  (Where are they allocated?)
        """
        self.lgrid = attributevector.AttributeVector()
        self.gsMap = GlobalSegMap.GlobalSegMap()
        if not name:
            self.name = "null"
        else:
            self.name = name
        if not suffix:
            self.suffix = "null"
        else:
            self.suffix = suffix
        if not n:
            self.n = 0
        else:
            self.n = n
        if  not ni:
            self.ni = 0
        else:
            self.ni = ni
        if not nj:
            self.nj = 0
        else:
            self.nj = nj
        self._indices = None
        return

    def __repr__(self):
        s = "<"
        s += self.name
        s += ","+self.suffix
        s += ","+str(self.n)
        s += ","+str(self.ni)
        s += ","+str(self.nj)
        s += ">"
        return s

    def __str__(self):
        return self.__repr__()

    def lsize(self):
        return self.gsMap.lsize( comm.local_comm )

    def info(self):
        print "\tdomain info for",self.name
        print "\tsuffix =",self.suffix
        print "\tn,ni,nj =",self.n, self.ni, self.nj
        print "\tdebugLevel =",debugLevel
        print "\tgsMap comp_id =",self.gsMap.compid()
        print "\tngSeg =",self.gsMap.ngsegs()
        print "\tgSize =",self.gsMap.globalsize()
        print "\tlSize =",self.gsMap.lsize( comm.local_comm )
        return
        
    def __del__(self):
        """
        Replaces Domain.Clean in the fortran interface.
        Cleans up MCT types and deallocates memory
        """
        return

    def initGrid( self, name, suffix, n, ni, nj, ifields=[],rfields=[], lsize=10):
        """
        Takes the following arguments and initializes Domain.av, an 
        attribute vector.

        name    - simply the name of this domain
        
        suffix  - a string to use a suffix on printouts
        n       - ni * nj, total number of points in this domain
        ni      
        nj      
        ifields - list of names of the integer fields to be defined for this domain, defaults to an empty list
        rfields - list of names of the real fields to be defined for this domain, defaults to an empty list
        lsize   - the size of all fields in this domain(integer)
        """
        #print "initGrid: initializing %s"%(name)
        self.name = name
        self.suffix = suffix
        self.n  = n
        self.ni = ni
        self.nj = nj
        self.ifields = ifields
        self.rfields = rfields
        #print "ifields=%s\nrfields=%s\nlsize=%d"%(temp_ifields, temp_rfields, lsize)
        self.lgrid = attributevector.AttributeVector( ifields, rfields, lsize )
        #print "allocating a temp array..."
        temp = Numeric.zeros( lsize, Numeric.Float64 )
        #temp = -9999.0
        #print "Filling real fields with default values..."
        for f in rfields:
            #print "\tFilling field",f,":",
            self.lgrid.importRAttr( f, temp )
            #print "... OK!"
        print "initGrid: Initialized Grid!"
        # setup av complete
        return

    def indices(self):
        if( self._indices == None ):
            print "*** warning:  indices undefined until after call to initGsMap! ***"
        return self._indices

    def initGsmap( self, lsize, indices ):
        # print "CPL:Entering initGsmap: lsize=%d"%(lsize)
        # Initialize contract's gsMap (based on index data)
        # print "... OK!"
        self._indices = indices
        if ( lsize == 0 ):
            nseg = 0
            start = Numeric.zeros(nseg,Numeric.Int32)
            count = Numeric.zeros(nseg,Numeric.Int32)
            #raise ValueError,"CPL:This process has a null/empty segment"
            print "This process has a null/empty segment"
        else:
            try:
                # print "CPL:gsMap_init, Compute segment's start indicies and length counts"
                nseg = 1 
                for n in xrange( 1, lsize ): # changed from 2 to 1
                    i = indices[n-1]
                    j = indices[n]
                    if ( j - i != 1 ):
                        nseg = nseg + 1
                # print "CPL:",time.ctime(),"Computed indices and counts."

                # print "CPL:(initGsmap) allocating Numeric arrays:"
                start = Numeric.zeros(nseg, Numeric.Int32)
                count = Numeric.zeros(nseg, Numeric.Int32)
                
                nseg = 1
                # nseg - 1 is an adjustment from Fortran to C array indices
                start[nseg-1] = indices[0] 
                count[nseg-1] = 1
                # print "CPL:(initGsmap) setting up start and count arrays:"
                for n in xrange(1, lsize):
                    i = indices[n-1]
                    j = indices[n]
                    if( (j - i) != 1 ):
                        nseg = nseg + 1
                        start[nseg-1] = indices[n] 
                        count[nseg-1] = 1
                    else:
                        count[nseg-1] = count[nseg-1] + 1
                # print "CPL:gsMap start[0],count[0]=",start[0],count[0],"start[nseg],count[nseg]=",start[nseg],count[nseg]
            except IndexError:
                print "i = ",i
                print "j = ",j
                print "nseg = ",nseg
                print "start = ",start
                print "count = ",count
                raise
            # start = Numeric.array([0,0],Numeric.Int32)
            # count = Numeric.array([0,lsize],Numeric.Int32)
            # print "CPL:Calling:",
        sys.stdout.flush()
        local_rank = mpi.comm_rank( comm.local_comm )
        print "CPL:self.gsMap.initd0(%s,%s,%s,%s,%s)"%(start,count, 0, comm.local_comm, comm.mph_component_id)
        # gsMap.initd0( Array start, Array count, Int root, Comm comm, Int component_id )
        self.gsMap.initd0( start, count, 0, comm.local_comm, comm.mph_component_id )
        # print "CPL:self.gsMap.ngsegs()=%s"%(self.gsMap.ngsegs())
        # print "CPL:self.gsMap.nlseg()=%s"%(self.gsMap.nlseg(self.pid))
        # print "CPL:Returning from gsMap.initd0!"
        sys.stdout.flush()
        self.__initialized_gsmap=True
        return


    def compare(self, other, enforce_mask=False, enforce_grid=False,
                enforce_area=False, enforce_aream=False, enforce_all=False):
        """
        Compares two domains(self and other), optionally
        aborting if the domains are too different.
        
        Optional Abort Conditions:
        enforce_mask: Abort if masks differ with respect to zero/nonzero
        enforce_grid: Abort if grids differ by eps_grid
        enforce_area: Abort if area differs by eps_area
        enforce_aream: Abort if area differs by eps_area
        enforce_all: abort for all of the above conditions
        
        NOTE:
          eps_grid and eps_area are defined in this function. 
          Not terribly modular or configurable.
        """
        eps_mask = 1.0e-6
        eps_grid = 1.0e-2
        eps_area = 1.0e-1

        # Do a global gather to create a non-distributed attribute vector
        debugPrint( "self.lgrid:\n",self.lgrid )
        debugPrint( "other.lgrid:\n",other.lgrid )
        gGrid1 = attributevector.AttributeVector(self.ifields, self.rfields, self.lsize())
        gGrid1.initv(self.lgrid, self.lgrid.lsize())
        gGrid1.gather(self.lgrid, self.gsMap, comm.world_pe0, comm.component_pid, comm.local_comm) 
        gGrid2 = attributevector.AttributeVector(other.ifields, other.rfields, other.lsize())
        gGrid2.initv( other.lgrid, other.lgrid.lsize() )
        gGrid2.gather(other.lgrid, self.gsMap,comm.world_pe0, comm.component_pid, comm.local_comm)

        # From here on, everything is done by the root pe
        if( comm.component_pid != comm.world_pe0 ):
            return

        # Compare size of domain
        npts1 = gGrid1.lsize()
        npts2 = gGrid2.lsize()
        npts = npts1

        if ( npts1 == npts2 ):
            debugPrint( "the domain size is ",npts )
        else:
            debugPrint( "domain size #1 = ", npts1 )
            debugPrint( "domain size #2 = ", npts2 )
            debugPrint( "ERROR: domain size mis-match" )
            # call shr_sys_abort(subName // "ERROR: domain size mis-match")
            # Exceptions?

        # If there was no problem, continue:
        # Compare Domain masks:
        debugPrint("gData1:\n",gGrid1)
        debugPrint("gData2:\n",gGrid2)
        data1,data1_size = gGrid1.exportRAttr("mask")#rcode)?
        data2,data2_size = gGrid2.exportRAttr("mask")#rcode)?
        
        ndiff = 0
        debugPrint( "npts:",npts )
        debugPrint( "length of data1:",data1_size )
        for n in xrange(0,npts-1):
            if ( (( (abs(data1[n])) > eps_mask ) and (abs(data1[n]) < eps_mask )) or 
                  ( (( abs(data1[n])) < eps_mask ) and (( abs(data1[n])) > eps_mask) ) ):
                ndiff = ndiff + 1

        # Enforce consistency: 
        # Nested function declaration
        def enforce_consistency(msg,exception=None):
            if (enforce_mask or enforce_all):
                if (ndiff > 0):
                    debugPrint( msg )
                    # Raise Exception
                    
        enforce_consistency("ERROR: incompatible domain masks")
                    
        # Compute Maximum Latitude and Longitude Differences
        mask = data1
        ndiff = 0
        data1,data1_size = gGrid1.exportRAttr("lat")#,rcode))
        data2,data2_size = gGrid2.exportRAttr("lat")#,rcode))
        diff = 0
        max_diff = 0.0
        for n in xrange(npts):
            if( abs( mask[n] ) > eps_mask ):
                diff = abs( data1[n] - data2[n] )
            max_diff = max(max_diff, diff)
            if( diff > eps_grid ):
                ndiff = ndiff + 1
        debugPrint( "Maximum latitude difference = ",max_diff )

        data1,data1_size = gGrid1.exportRAttr("lon")#,rcode))
        data2,data2_size = gGrid2.exportRAttr("lon")#,rcode))
        max_diff = 0.0

        for n in xrange(npts):
            if( abs( mask[n] ) > eps_mask ):
                x1 = data1[n]
                x2 = data2[n]
                if( x1 > x2 ): #make sure x1 < x2
                    # swap(x1,x2)
                    x1 = data2[n]
                    x2 = data1[n]
                while( (x1+360.0) < (x2+180.0) ):#longitude is periodic
                    x1 = x1 + 360.0
                diff = abs( x2 - x1 )
                max_diff = max(max_diff,diff)
    
                if (diff > eps_grid):
                    ndiff = ndiff + 1
        debugPrint( "Maximum longitude difference = ",max_diff )

        enforce_consistency("ERROR: incompatible domain grid coordinates!")

        # Compare Area:
        data1,data1_size = gGrid1.exportRAttr( "area" )#, rcode )
        data2,data2_size = gGrid2.exportRAttr( "area" )#, rcode )

        ndiff = 0
        max_diff = 0.0

        for n in xrange(npts):
            if( abs( mask[n] ) > eps_mask ):
                if( data2[n] != 0.0 ):
                    diff = abs( (data2[n] - data1[n]) / data2[n] )
                max_diff = max(max_diff,diff)
                if( diff > eps_area ):
                    ndiff = ndiff + 1
        debugPrint( "Maxium relative error of area (model) = ", max_diff )

        enforce_consistency("ERROR: icompatible domain area(model)")

        # Compare aream
        data1,data1_size = gGrid1.exportRAttr("aream")#,rcode))
        data2,data2_size = gGrid2.exportRAttr("aream")#,rcode))

        ndiff = 0
        max_diff = 0.0
        for n in xrange(npts):
            if ( abs( mask[n] ) > eps_mask ):
                if( data2[n] != 0.0 ):
                    diff = abs((data2[n] - data1[n])/data2[n])
                max_diff = max(max_diff,diff)
                if( diff > eps_area ):
                    ndiff = ndiff + 1
        debugPrint( "maximum relative error of area(map) = ",max_diff )

        enforce_consistency("ERROR: incompatible domain area (map)")

        # Clean up, we're finished!
        return

#     def decompose( self, decomp, gi,gj, myid, npes ):
#         """
#         lsize = recvcontract.decompose(decomp, gi, gj, myid, npes, indx)
    
#         This function creats a decomposition given a global grid
#         and a decomposition type.  This is used for contract
#         initialization on the recv side.
        
#         This method takes the following arguments:
        
#         decomp   # the decomposition type(integer)
#         gi,gj    # global i and j sizes(integer)
#         myid     # this processes pe number(integer)
#         npes     # total number of pes(integer)
       
#         And returns a tuple containing: 
#         lsize    # local size of decomposition(integer)
#         indx     # global index of decomposition(integer array)
        
#         """
#         gsize = gi*gj
#         indx = None
#         lsize = 0
#         if( decomp == 1 ): # Assume 1D Decomposition in i direction
#             n = 0
#             nx = gi
#             ny = gj / npes
#             lsize = ny * nx
#             ng = 0
#             nglast = 0
#             print "decompose: nx = %s\tny = %s\tlsize = %s"%(nx,ny,lsize)
#             indx = Numeric.zeros( lsize, Numeric.Float32 )
#             for j in xrange(1,ny+1):
#                 for i in xrange(1,nx+1):
#                     ig = i
#                     jg = j + (myid*ny)
#                     nglast = ng
#                     ng = ((jg - 1)*gi) + ig
#                     if False:
#                         if ( ng > (nglast+1) ):
#                             print "nglast = %s - ng = %s" % (nglast, ng)
#                     indx[n] = ng 
#                     n += 1
#         elif( decomp == 2 ): # Assume 1D decomposition in the j direction
#             raise ContractError,"ReceivingContract.decompose: Invalid or Unavailable Decomposition!"
#         elif( decomp == 901 ): # test decomposition, do not use
#             raise ContractError,"ReceivingContract.decompose: Invalid or Unavailable Decomposition!"
#         else:
#             raise ContractError,"ReceivingContract.decompose: Invalid or Unavailable Decomposition!"
        
#         # Sort 'indx' for mapping performance?
#         # Skipping this sort for now!
#         return lsize,indx
    
