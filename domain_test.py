import sys
import Numeric

import mpi

import cpl
import cpl.domain

import MCT
from MCT import World
from MCT import Router
from MCT import GlobalSegMap

rank,size = mpi.init( len(sys.argv), sys.argv )

print "Shared Variables Initializing..."
mct_world = World.World()
r = Router.Router()
gsm = GlobalSegMap.GlobalSegMap()
start = Numeric.array([0,0],Numeric.Int32)
count = Numeric.array([0,100],Numeric.Int32)
print "Shared Variables Initialized!"

if( rank == 0 ):
    ds = cpl.domain.Domain()
    ds.initGrid( "bob", "q", 100, 10, 10, [],cpl.fields.grid_fields, 100 )
    mycomm = mpi.comm_split( mpi.MPI_COMM_WORLD, 0)
    myrank = mpi.comm_rank( mycomm )
    print "Before MCT World Init..."
    mct_world.initd0( 2, mpi.MPI_COMM_WORLD, mycomm, 1 )
    print "Before GSM Init..."
    gsm.initd0( start, count, myrank, mycomm, 1 )
    print "Before Router Init..."
    r.init( 2, gsm, mycomm )
    print "Before Domain Send..."
    ds.send( r, 7 )
else:
    dr = cpl.domain.Domain()
    dr.initGrid( "bob", "q", 100, 10, 10, [],cpl.fields.grid_fields, 100 )
    mycomm = mpi.comm_split( mpi.MPI_COMM_WORLD, 1)
    myrank = mpi.comm_rank( mycomm )
    print "Before MCT World Init..."
    mct_world.initd0( 2, mpi.MPI_COMM_WORLD, mycomm, 2 )
    print "Before GSM Init..."
    gsm.initd0( start, count, myrank, mycomm, 2 )
    print "Before Router Init..."
    r.init( 1, gsm, mycomm )
    print "Before Domain Send..."
    dr.recv( r, 7 )

print "DONE!"
mpi.finalize()
print "EXITING!"
sys.exit(0)
