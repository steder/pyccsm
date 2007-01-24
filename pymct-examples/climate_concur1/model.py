#model.py
"""
A generic model object for testing MCT functionality.
"""

# standard python
import sys, os, math

# external dependencies:
import Numeric, mpi

# PyMCT modules:
import MCT
import MCT.World
import MCT.GlobalSegMap
import MCT.AttrVect
import MCT.Router

class Model:
    nxa = 128
    nya = 64
    def __init__(self, comm, ncomps, compid, cplid):
        self.rank = mpi.comm_rank(comm)
        self.size = mpi.comm_size(comm)
        root = 0

        if(self.rank==0):
            print "model.py: rank:",self.rank
            print "model.py: size:",self.size

        ### initialize MCTWorld:
        self.world = MCT.World.World()
        self.world.initd0(ncomps, mpi.MPI_COMM_WORLD,
                          comm, compid)

        ### Initialize a global segment map:
        #
        # Set up a 1D decomposition
        # 1 segment per processor
        #
        localsize = (self.nxa * self.nya) / self.size
        # we'll use the distributed init of GSMap
        # so we'll need start and length arrays for
        # this processor
        start = Numeric.zeros((1),Numeric.Int32)
        length = Numeric.zeros((1),Numeric.Int32)
        start[0] = (self.rank*localsize) + 1
        length[0] = localsize

        self.GSMap = MCT.GlobalSegMap.GlobalSegMap()
        self.GSMap.initd0(start, length, root, comm, compid)

        # get the points local to this processor in
        # their assumed order:
        self.points = self.GSMap.OrderedPoints( self.rank )

        ### Initialize an Attribute vector:
        avsize = self.GSMap.lsize( comm )
        if( self.rank == 0 ):
            print "model.py: localsize:",avsize
        self.av = MCT.AttrVect.AttrVect()
        print "model.py: initializing attribute vector..."
        self.av.init("","field1:field2",avsize)

        ### Initialize a router to the coupler component:
        print "model.py: initializing router..."
        self.router = MCT.Router.Router()
        self.router.init( cplid, self.GSMap, comm )

        ### Create an array used in run method:
        self.avdata = Numeric.zeros((avsize),Numeric.Float64)
        print "model.py: successfully initialized!"
        return

    def run(self, steps=10):
        # timestep loop:
        print "model.py: running for",steps,"timesteps..."
        for i in xrange(steps):
            print "\tmodel.py: starting step",i,":"
            self.step(i)
        return

    def step(self,step):
        ### load data into av
        # the first field will be a constant real number
        self.avdata[:] = 30.0
        self.av.importRAttr("field1",self.avdata)
        # the second field will be the indices of each
        # grid point in the grid point numbering scheme
        for n in xrange(len(self.avdata)):
            self.avdata[n] = self.points[n]
        self.av.importRAttr("field2",self.avdata)

        ### Send the data
        # this is a synchronization point between
        # the coupler and this model
        if(self.rank==0):
            print "\tmodel.py: sending data on step",step
        # router, tag
        self.av.send(self.router,0)
        return
    
    def __del__(self):
        print "model.py: done!"
        return
