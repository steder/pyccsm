"""
file: frac.py -- handles surface fractions

   Defines, declares, initializes, and updates surface fractions.

REMARKS:
    These fractions are used for merging fields onto various domains.
    This particular implementation of this module makes certain assumptions 
    about which domains exist and the relationships between them.  These 
    assumptions are hard-coded into this software implementation.

 ASSUMPTIONS: 
    o atm & lnd grid cells and domain decompostion are identical
    o ice & ocn domains are identical (includes decomposition)
    o all atm cells are fully active
    o all ocn cells are either fully active or fully inactive
    o lnd cells can be partially active -- the fraction of a lnd cell that is 
      active is the fraction that is not occupied by ocn
    o ice cells can be partially active -- the fraction that is active is
      determined by the ice component itself
     
    For each domain (atm,lnd,ice,ocn) there are four fractions: fa,fi,fl,fo,
    three that could be used for merging, and one which indicatates the 
    fraction of the cell which is active.
    o merging on atm domain: Fa = fi*Fi + fl*Fl + fo*Fo   (fi + fl + fo = 1)
    o merging on ice domain: Fi = fa*Fa + fo*Fo = Fa + Fo (fa = fo = 1, fl=0)
    o merging on lnd domain: Fl = fa*Fa = Fa              (fa = 1, fo = fi = 0)
    o merging on ocn domain: Fo = fa*Fa + fo*Fo           (fa + fi = 1, fl = 0)
    o on the atm domain: fa = 1 (atm cells are fully active)
    o on the ice domain: fi = is time-variant and determined by the ice model
    o on the lnd domain: fl = 1 - fo and is time-invariant
    o on the ocn domain: fo = 1 (ocn cells are fully active)
"""
## Standard Library imports
import sys, os, time

##
import mpi
from MCT import AttrVect
##
import Numeric
##
import cpl
import cpl.comm
import cpl.attributevector
import cpl.domain
import cpl.fields
import cpl.const
import cpl.control
##
from cpl import debug
debugPrint = debug.newDebugPrint(False)
statusPrint = debug.newDebugPrint(False)

frac_fields = ["afrac",
               "ifrac",
               "lfrac",
               "ofrac",
               ]

frac_a = None
frac_i = None
frac_l = None
frac_o = None

def nint( num ):
    return int(round(num))

def init( map_o2a, domain_a, domain_i, domain_l, domain_o ):
    """
    ROUTINE: frac_init - initialize the surface fraction bundles

    DESCRIPTION:
    Initialize the fraction bundles.  All fractions are derived from the
    (time-invariant) ice/ocn domain masks plus the (time-variant) ice fraction.
    This initialization routine sets the time-invariant values.
    """
    ## Initialize fraction "bundles" with spvals.
    # lsize in this case is referring to the size
    # of the domain's GlobalSegMap.
    print "domain_a.lsize()=",domain_a.lsize()
    print "domain_i.lsize()=",domain_i.lsize()
    print "domain_l.lsize()=",domain_l.lsize()
    print "domain_o.lsize()=",domain_o.lsize()
    global frac_a, frac_i, frac_l, frac_o
    frac_a = cpl.attributevector.AttributeVector([],
                                                 frac_fields,
                                                 domain_a.lsize())
    frac_i = cpl.attributevector.AttributeVector([],
                                                 frac_fields,
                                                 domain_i.lsize())
    frac_l = cpl.attributevector.AttributeVector([],
                                                 frac_fields,
                                                 domain_l.lsize())
    frac_o = cpl.attributevector.AttributeVector([],
                                                 frac_fields,
                                                 domain_o.lsize())
    
    # Initialize frac_* to cpl.const.spval:
    # below is written in fortran as
    # frac_a.rAttr = cpl.const.spval
    default_a = Numeric.zeros( (domain_a.lsize()), Numeric.Float64 )
    default_i = Numeric.zeros( (domain_i.lsize()), Numeric.Float64 )
    default_o = Numeric.zeros( (domain_o.lsize()), Numeric.Float64 )
    default_l = Numeric.zeros( (domain_l.lsize()), Numeric.Float64 )
    default_a[:] = cpl.const.spval
    #default_i[:] = cpl.const.spval # ice is initialized to zeros
    default_o[:] = cpl.const.spval
    default_l[:] = cpl.const.spval
    for field in frac_fields:
        frac_a.importRAttr( field, default_a ) #initialized to all cpl.const.spval
        frac_i.importRAttr( field, default_i ) #initialized to all 0's
        frac_o.importRAttr( field, default_o ) #initialized to all cpl.const.spval
        frac_l.importRAttr( field, default_l ) #initialized to all cpl.const.spval
    # Here they set bundle counters(for accumulation)
    # we won't be doing it here.

    ## initialize values on ice grid (based on zero ice fraction)
    # already initialized -- frac_i.rdata[:] = 0.0
    mask,mask_size = domain_i.lgrid.exportRAttr("mask")
    debugPrint("ICE DOMAIN Mask slice[0:400]:\n",mask[0:400])
    # where( nint( domain_i.rdata[km][:] ) != 0 ):
    #     frac_i.rdata[ka][:] = 1.0
    #     frac_i.rdata[ko][:] = 1.0
    #frac_i.rdata[ka] = Numeric.where( (domain_i.lgrid.rdata[km] != 0),
    #                                  1.0, frac_i.rdata[ka] )
    #frac_i.rdata[ko] = Numeric.where( (domain_i.lgrid.rdata[km] != 0),
    #                                   1.0, frac_i.rdata[ko] )
    ice_afrac,size = frac_i.exportRAttr("afrac")
    ice_ofrac,size = frac_i.exportRAttr("ofrac")
    for i,v in enumerate(mask):
        if ( nint(v) != 0 ):
            ice_afrac[i] = 1.0
            ice_ofrac[i] = 1.0
    frac_i.importRAttr("afrac",ice_afrac)
    frac_i.importRAttr("ofrac",ice_ofrac)
    
    ## initialize values on ocean grid( same as for ice grid)
    #frac_o.rdata = frac_i.rdata
    for field in frac_fields:
        tmp,tmpsize = frac_i.exportRAttr(field)
        frac_o.importRAttr(field,tmp)

    ## initialize values on atm grid (needs/assumes zero ice)
    
    # map all fractions from ocn to atm grid:
    #call cpl_map_bun(bun_frac_o,bun_frac_a,map_o2a)
    frac_a = map_o2a.mapAV(frac_o)
    
    # clean up atm fraction: must be 1 everywhere
    ones = Numeric.ones( (frac_a.lsize()), Numeric.Float64 )
    frac_a.importRAttr( "afrac", ones )

    # clean up ice faction: must be 0 everywhere:
    zeros = Numeric.zeros( (frac_a.lsize()), Numeric.Float64 )
    frac_a.importRAttr("ifrac",zeros)

    # clean up ocn fraction: must be in [0,1]:
    ofrac,ofrac_size = frac_a.exportRAttr("ofrac")
    debugPrint("ATM ofrac:\n",ofrac[0:400])
    # Numeric.where( condition, x, y )
    # Where condition is true, we set the element to x, if false, y.
    #frac_a.rdata[ko] = Numeric.where( (frac_a.rdata[ko] > 1.0),
    #                                  1.0, frac_a.rdata[ko] )
    #frac_a.rdata[ko] = Numeric.where( (frac_a.rdata[ko] < 0.0),
    #                                  0.0, frac_a.rdata[ko] )
    for i,v in enumerate( ofrac ):
        if( v > 1.0 ):
            ofrac[i] = 1.0
        if( v < 0.0 ):
            ofrac[i] = 0.0
    frac_a.importRAttr("ofrac",ofrac)
    
    # compute land fraction: lnd = 1 - ocn, then clean it up:
    lfrac,lfrac_size = frac_a.exportRAttr("lfrac")
    lfrac = 1.0 - ofrac
    debugPrint("lfrac = 1 - ofrac:\n",lfrac[0:400])
    #frac_a.rdata[kl] = Numeric.where( (frac_a.rdata[kl] > 1.0),
    #                                  1.0, frac_a.rdata[kl] )
    #frac_a.rdata[kl] = Numeric.where( (frac_a.rdata[kl] < 0.0),
    #                                  0.0, frac_a.rdata[kl] )
    for i,v in enumerate( lfrac ):
        if( v > 1.0 ):
            lfrac[i] = 1.0
        if( v < 0.0 ):
            lfrac[i] = 0.0
        """
        NOTE: it's a requirement to elminate land points smaller than .001  ---
        this is to avoid active land cells with tiny land fractions   ---
        NOTE: this may result in some, presumably small, non-conservation   ---
        NOTE: should probably count the number of pnts that get set to zero ---
        """
        #frac_a.rdata[kl] = Numeric.where( (frac_a.rdata[kl] < 0.001),
        #                              0.0, frac_a.rdata[kl] )
        if( v < 0.001 ):
            lfrac[i] = 0.0
            
    frac_a.importRAttr("lfrac",lfrac)
    ## Initialize values on land grid: same as
    ## atm except ice,ocn = 0 
    for field in frac_fields:
        tmp,tmpsize = frac_a.exportRAttr(field)
        frac_l.importRAttr(field,tmp)
        
    frac_l.importRAttr("ifrac",zeros)
    frac_l.importRAttr("ofrac",zeros)

    lfrac,lfrac_size = frac_a.exportRAttr( "lfrac" )
    debugPrint("frac_a lfrac at the end of frac.init()\n",
               lfrac[0:400])
    return

def set(ifrac_i, map_o2a, domain_a, domain_i, domain_l, domain_o):
    """
    ROUTINE: set - set/update the surface fraction bundles

    DESCRIPTION:
    Set/update the fraction bundles based on input (time-variant) ice fraction 
    and time-invariant land fraction. This set/update routine sets the 
    time-variant values.  The companion initialization routine must be called 
    first to set the time-invariant values.
    """
    """
    Note: we assume that indifrac_is are the same for
    ALL fraction vectors/bundles
    """
    global frac_a, frac_i, frac_l, frac_o
    ## Check for erroneous ifrac_i fractions
    #km = domain_i.lgrid.av.indexRA("mask")
    mask,mask_size = domain_i.lgrid.exportRAttr("mask")
    #ki = frac_i.av.indexRA("ifrac")

    # We want to use ifrac_i that is passed in above, not the incorrect values in frac_i.
    #ifrac_i,ifrac_i_size = frac_i.exportRAttr("ifrac")
    
    if( cpl.control.infodbug > 1):
        maxl = 0.0
        minl = 1.0
        # Compute local min and max
        for n in xrange(frac_i.lsize()):
            if( nint( mask[n] ) != 0 ):
                maxl = max( max1, ifrac_i[n] )
                minl = min( minl, ifrac_i[n] )
        # Reduce to find global min and max
        maxg = mpi.reduce( maxl, 1, mpi.MPI_FLOAT,
                           mpi.MPI_MAX, 0, cpl.comm.local_comm )
        ming = mpi.reduce( minl, 1, mpi.MPI_FLOAT,
                           mpi.MPI_MIN, 0, cpl.comm.local_comm )
        maxg,ming = maxg[0],ming[0]
        #
        if( cpl.comm.component_pid == 0 ):
            if( maxg > 1.0 ):
                print "frac: set(): WARNING: global max ifrac =",maxg
            if( ming < 0.0 ):
                print "frac: set(): WARNING: global min ifrac =",ming
            
    ## Set/Update values on IFRAC_I grid, confine values into
    #  [0,1], mask = 0 -> frac = 0
    #print "len( frac_i.rdata[ki] ) = ",len(ifrac_i)
    #print "len( frac_i.rdata[km] ) = ",len(mask)
    #print "size of ifrac_i:",len(ifrac_i)
    
    for n in xrange(frac_i.lsize()):
        tmp = ifrac_i[n]
        ifrac_i[n] = min(1.0,tmp)
        ifrac_i[n] = max(0.0,tmp)
        if(nint( mask[n] )== 0):
            ifrac_i[n] = 0.0
            
    frac_i.importRAttr("ifrac",ifrac_i)
    ## Set/Update values on OCEAN grid
    #  (assumes ice and ocn have same domain)
    #ki = frac_i.av.indexRA("ifrac")
    # ifrac_i, from above
    #ka = frac_o.av.indexRA("afrac")
    ifrac_o,ifo_size = frac_o.exportRAttr("ifrac")
    afrac_o,afo_size = frac_o.exportRAttr("afrac")
    for n in xrange(frac_o.lsize()):
        ifrac_o[n] =       ifrac_i[n]
        afrac_o[n] = 1.0 - ifrac_i[n]
        afrac_o[n] = min( 1.0, afrac_o[n] )
        afrac_o[n] = max( 0.0, afrac_o[n] )
    frac_o.importRAttr("ifrac",ifrac_o)
    frac_o.importRAttr("afrac",afrac_o)
    ## Set/Update vales on ATMOS grid
    # map ifrac onto atm grid (other mapped fracs are not useful )
    #print "*** frac.py: NEED TO CALL CPL_MAP_BUN! ***"
    statusPrint( "*** frac.py: CALLING CPL_MAP_BUN! ***" )
    frac_a = map_o2a.mapAV(frac_o)
    
    #ka = frac_a.av.indexRA("afrac")
    #ki = frac_a.av.indexRA("ifrac")
    #kl = frac_a.av.indexRA("lfrac")
    #ko = frac_a.av.indexRA("ofrac")
    afrac_a,afa_size = frac_a.exportRAttr("afrac")
    ifrac_a,ifa_size = frac_a.exportRAttr("ifrac")
    lfrac_a,lfa_size = frac_a.exportRAttr("lfrac")
    ofrac_a,ofa_size = frac_a.exportRAttr("ofrac")
    lfrac_l,lfl_size = frac_l.exportRAttr("lfrac")
    for n in xrange( frac_a.lsize() ):
        # restore afrac
        afrac_a[n] = 1.0
        # restore lfrac
        lfrac_a[n] = lfrac_l[n]
        # clean up ifrac
        ifrac_a[n] = min(1.0, ifrac_a[n])
        ifrac_a[n] = max(0.0, ifrac_a[n])
        # compute ofrac = 1.0 - ifrac - lfrac
        ofrac_a[n] = ( 1.0 - ifrac_a[n] - lfrac_a[n] )
        # clean up ofrac
        ofrac_a[n] = min( 1.0, ofrac_a[n] )
        ofrac_a[n] = max( 0.0, ofrac_a[n] )
    ## ?? Do we need to do this?
    ## Set accumulation counts to 1
    ## No, we can skip this
    # Clean up:
    frac_a.importRAttr("afrac",afrac_a)
    frac_a.importRAttr("ifrac",ifrac_a)
    frac_a.importRAttr("lfrac",lfrac_a)
    frac_a.importRAttr("ofrac",ofrac_a)
    return
