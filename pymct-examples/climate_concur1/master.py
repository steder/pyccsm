#master.py
"""
Driver for a simple concurrent coupled model

Provide a simple example of using MCT components
executing concurrently in a single executable.
"""

# standard python modules:
import sys, os

# non-standard external libraries:
import Numeric
import mpi

# PyMCT modules:
import MCT

# Example modules:
import model
import coupler

# CONSTANTS:
NCOMPONENTS=2
ATMID=1
CPLID=2
FILENAME=""

def main():
    """
    This program sets up a single-executable,
    concurrent-execution system.

    Main carves up the MPI_COMM_WORLD and starts
    each component on its own processor set.
    """
    rank,size = mpi.init()

    # set color:
    if( rank < (size/2) ):
        color = 0
    else:
        color = 1

    # create subcommunicators:
    comm = mpi.comm_split( mpi.MPI_COMM_WORLD, color )

    if( color == 0 ):
        mymodel = model.Model(comm, NCOMPONENTS, ATMID,CPLID)
        mymodel.run()
    elif(color == 1):
        mycoupler = coupler.Coupler(comm, NCOMPONENTS, CPLID, ATMID, FILENAME)
        mycoupler.run()

    mpi.finalize()

if __name__=="__main__":
    main()
