#coupler.py
"""
A generic coupler object for testing MCT functionality.
"""

# standard python
import sys, os, math

# external dependencies:
import Numeric, mpi
import pycdf

# PyMCT modules:
import MCT
import MCT.World
import MCT.GlobalSegMap
import MCT.AttrVect
import MCT.Router
import MCT.SparseMatrix
import MCT.SparseMatrixPlus

# Constants:
AREA_FIELD = "aream"

class Coupler:
    # Set grid dimensions for atmosphere and ocean grids.
    # MCT could be used for this by defining a GeneralGrid
    # on each component model and sending them to the
    # coupler.
    #   However, for this example we'll just assume
    #   they are known to the coupler:
    nxa = 128
    nya = 64
    nxo = 320
    nyo = 384
    def __init__(self, comm, ncomps,
                 compid, modelid, filename):
        ### Get local rank and size:
        self.rank = mpi.comm_rank( comm )
        self.size = mpi.comm_size( comm )
        root = 0
        if (self.rank==0):
            print "coupler.py: rank:",self.rank
            print "coupler.py: size:",self.size
        ### Initialize MCTworld:
        self.world = MCT.World.World()
        self.world.initd0(ncomps, mpi.MPI_COMM_WORLD,
                          comm, compid)

        if(self.rank==root):
            ### Read matrix weights for interpolation
            # from a file:
            self.sMatGlobal = MCT.SparseMatrix.SparseMatrix()
            self.areasrcGlobal = MCT.AttrVect.AttrVect()
            self.areadstGlobal = MCT.AttrVect.AttrVect()
            self.read( filename )

        ### Initialize a global segment map for the ocean
        # set up a 1D decomposition
        # just one segment per processor
        localsize = (self.nxo * self.nyo) / self.size
        # we'll use the distributed init of GSMap
        # so initialize start and length arrays for
        # this processor.
        start = Numeric.zeros((1),Numeric.Int32)
        length = Numeric.zeros((1),Numeric.Int32)
        start[0] = (self.rank * localsize) + 1
        length[0] = localsize
        # initialize the GSMap
        print "coupler.py: initializing Ocean GSMap..."
        self.OcnGSMap = MCT.GlobalSegMap.GlobalSegMap()
        self.OcnGSMap.initd0( start, length, root,
                           comm, compid )
        ### Initialize the global segment map for the atmosphere
        # set up a 1D decomposition
        # there is just 1 segment per processor
        localsize = (self.nxa * self.nya) / self.size
        # we're using the distributed init of GSMap so
        # initialize start and length array sfor this processor
        start = Numeric.zeros((1),Numeric.Int32)
        length = Numeric.zeros((1),Numeric.Int32)
        start[0] = (self.rank * localsize) + 1
        length[0] = localsize
        # initialize the gsmap:
        print "coupler.py: initializing Atmosphere GSMap..."
        self.AtmGSMap = MCT.GlobalSegMap.GlobalSegMap()
        self.AtmGSMap.initd0( start, length, root,
                              comm, compid )

        ### Get points local to this processor
        self.points = self.AtmGSMap.OrderedPoints(self.rank)

        ### Build a SparseMatrixPlus for doing
        # the interpolation.  Specify matrix
        # decomposition to be by row following the
        # atmosphere's decomposition.
        self.A20MatPlus = MCT.SparseMatrixPlus.SparseMatrixPlus()
        """
        Needs read to work:
        self.A20MatPlus.initFromRoot(self.sMatGlobal,
                              self.AtmGSMap,
                              self.OcnGSMap,
                              #\/Strategy:"Xonly","Yonly","XandY"
                              "Xonly",
                              root,
                              comm,
                              compid)
        """
        ### Initialize an attribute vector for the ATM grid:
        self.aavsize = self.AtmGSMap.lsize(comm)
        if(self.rank==root):
            print "coupler.py: localsize: ATM:",self.aavsize
        self.atmav = MCT.AttrVect.AttrVect()
        self.atmav.init("","field1:field2",self.aavsize)
        ### Initialize an attribute vector for the OCN grid:
        self.oavsize = self.OcnGSMap.lsize(comm)
        if(self.rank==root):
            print "coupler.py: localsize: OCN:",self.oavsize
        self.ocnav = MCT.AttrVect.AttrVect()
        self.ocnav.init("","field1:field2",self.oavsize)
        ### Initialize a router:
        self.router = MCT.Router.Router()
        self.router.init(modelid, self.AtmGSMap, comm)
        return

    def run(self, steps=10):
        for i in xrange(steps):
            self.step(i)
        return

    def step(self, step):
        self.match = True
        # router, tag, sum(boolean)
        self.atmav.recv(self.router, 0, False)
        # The 2nd attribute of atmav
        # has the values of each gridpoint in the
        # index numbering scheme.  Check the received
        # values against the points on this processor.
        #
        # They should match exactly:
        recvpoints,recvsize = self.atmav.exportRAttr("field2")
        for n in xrange(self.aavsize):
            if (recvpoints[n] != self.points[n]):
                print "\t** coupler.py: data element",n,"doesn't match! **"
                self.match=False
        print "\tcoupler.py: received data step",step

        # Interpolate by doing a parallel
        # sparsematrix-attrvect multiply.
        # NOTE: it doesn't make sense to interpolate
        # field2 (grid point indices) but MatVecMul
        # will interpret ALL real attributes.
        #self.atmav.sMatAvMult_sMPlus(self.ocnav,
        #                             self.A20MatPlus)

        # Pass interpolated data on to ocean model
        # and/or do more calculations
        return

    def read( self, filename ):
        """
        Modifies the internal state of this map object to contain 
        all the appropriate data gathered from the netcdf mapping 
        file whose path is provided above as the filename argument.
        """
        print "coupler.py: reading file..."
        return

    def __del__(self):
        print "coupler.py: done."
        return
    
