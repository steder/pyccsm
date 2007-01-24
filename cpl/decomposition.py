"""
decompose.py

Defines some common decompositions used by both the coupler and component models.
"""
import Numeric
from error import CPLException

class DecompositionError(CPLException):   
    def __init__(self, value="No error message provided"):
        Exception.__init__(self)
        self.value = value  
    def __str__(self):
        return repr(self.value)

def decompose( decomp, gi,gj, myid, npes ):
    """
    lsize = recvcontract.decompose(decomp, gi, gj, myid, npes, indx)
    
    This function creats a decomposition given a global grid
    and a decomposition type.  This is used for contract
    initialization on the recv side.
    
    This method takes the following arguments:
    
    decomp   # the decomposition type(integer)
    gi,gj    # global i and j sizes(integer)
    myid     # this processes pe number(integer)
    npes     # total number of pes(integer)
    
    And returns a tuple containing: 
    lsize    # local size of decomposition(integer)
    indx     # global index of decomposition(integer array)
    
    """
    if( decomp == 1 ): # Assume 1D Decomposition in i direction
        lsize,indx = decomp_1d_x( gi, gj, myid, npes )
    elif( decomp == 2 ): # Assume 1D decomposition in the j direction
        raise DecompositionError,"decompose: Invalid or Unavailable Decomposition!"
    elif( decomp == 901 ): # test decomposition, do not use
        raise DecompositionError,"decompose: Invalid or Unavailable Decomposition!"
    else:
        raise DecompositionError,"decompose: Invalid or Unavailable Decomposition!"
    
    # Sort 'indx' for mapping performance?
    # Skipping this for now.
    return lsize,indx

def decomp_1d_x( gi, gj, myid, npes ):
    """
    lsize,indx = decomp_1D_x(gi,gj,myid,npes)
    
    Computes a 1D decomposition in X.
    """
    gsize = gi*gj
    indx = None
    lsize = 0
    n = 0
    nx = gi
    ny = gj / npes
    lsize = ny * nx
    ng = 0
    nglast = 0
    print "decompose: nx = %s\tny = %s\tlsize = %s"%(nx,ny,lsize)
    indx = Numeric.zeros( lsize, Numeric.Float32 )
    for j in xrange(1,ny+1):
        for i in xrange(1,nx+1):
            ig = i
            jg = j + (myid*ny)
            nglast = ng
            ng = ((jg - 1)*gi) + ig
            if False:
                if ( ng > (nglast+1) ):
                    print "nglast = %s - ng = %s" % (nglast, ng)
            indx[n] = ng 
            n += 1
    return lsize,indx
