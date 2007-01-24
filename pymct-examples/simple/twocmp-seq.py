#! /usr/bin/env python
# twocmp-seq.py
"""
Mike Steder

PyMCT:
two component, sequential example script

We are implementing a single-executable, sequential-execution, multiple
(in this case 2) component system.

In this example communication occurs through main via arguments to the
modelX methods.  Both components share the same processors.
"""

# General Python Imports
import os,sys,math

# Useful External Modules
import Numeric
import mpi

# PyMCT Imports
import MCT
import MCT.World
import MCT.GlobalSegMap
import MCT.AttrVect
import MCT.Rearranger

# Constants:
NGX,NGY = 6,4
NCOMPONENTS = 2
COMPONENT_IDS = [1,2]

def model1(comm, av):
    # Find local rank and size
    size = mpi.comm_size( comm )
    rank = mpi.comm_rank( comm )
    print "model 1 size:",size
    avsize = av.lsize()
    print "model 1 av size:",avsize

    # Fill av with some data:
    data1 = [ i + (20*rank) for i in range(avsize) ]
    data1 = Numeric.array(data1, Numeric.Int32)
    av.importRAttr("field1",data1)
    data2 = [ math.cos( i * (20*rank) ) * 3.14 for i in range(avsize)]
    av.importRAttr("field2",data2)
    
    return av

def model2(comm, av):
    # Find local rank and size
    size = mpi.comm_size( comm )
    rank = mpi.comm_rank( comm )
    print "model 1 size:",size
    avsize = av.lsize()
    print "model 1 av size:",avsize

    # print out AV data:
    print "model 2 data after:"
    for f in ["field1","field2"]:
        d = av.exportRAttr(f)
        print d
    return

def twoseq():
    # --- BEGIN INITIALIZATION ---
    # Init MPI
    rank,size = mpi.init()
    print "creating sub-communicators:"
    comm1 = mpi.comm_dup( mpi.MPI_COMM_WORLD )
    comm2 = mpi.comm_dup( mpi.MPI_COMM_WORLD )
    print "ok!"
    
    # Init World:
    myids = Numeric.array([1,2],Numeric.Int32)
    # print "initializing world:",
    world = MCT.World.World()
    """
     Help on built-in function initd1:

 initd1(...)
     initd1(in int ncomps,
            in int globalcomm,
            in int mycomm,
            in array<int> myids)
     RETURNS
         None
     RAISES
         sidl.RuntimeException

    """
    print "myids:",myids
    world.initd1(NCOMPONENTS, mpi.MPI_COMM_WORLD, comm1, myids)
    print "okay!"
    # Set up a grid and decomposition, first gsmap is the grid
    # decomposed by rows
    # There's 1 segment per processor
    length1 = Numeric.zeros((1),Numeric.Int32)
    length1[0] = NGX * (NGY / size)
    start1 = Numeric.zeros((1),Numeric.Int32)
    start1[0] = rank * (length1[0] + 1)
    print "initializing gsmap: %d,%d,%d"%(rank,length1[0],start1[0])
    GSMap1 = MCT.GlobalSegMap.GlobalSegMap()
    # 0 is the root mpi processor, 1 is the component id
    GSMap1.initd0(start1, length1, 0, comm1, 1)

    # Second GSMap is the grid decomposed by columns:
    length2 = Numeric.zeros((NGY),Numeric.Int32)
    start2 = Numeric.zeros((NGY),Numeric.Int32)
    for i in range(NGY):
        length2[i] = NGX/size
        start2[i] = (i-1)*NGX + 1 + rank*length2[i]
        print "gsmap2: %d,%d,%d,%d"%(rank,i,length2[i],start2[i])
    GSMap2 = MCT.GlobalSegMap.GlobalSegMap()
    # 0 is the root mpi processor, 2 is the component id
    GSMap2.initd0(start2,length2,0,comm2,2)

    # Initialize some AttrVects
    av1 = MCT.AttrVect.AttrVect()
    av1.init("","field1:field2",GSMap1.lsize(comm1))
    av2 = MCT.AttrVect.AttrVect()
    av2.init("","field1:field2",GSMap2.lsize(comm2))

    # Create the rearranger
    Rranger = MCT.Rearranger.Rearranger()
    Rranger.init( GSMap1, GSMap2, mpi.MPI_COMM_WORLD )
    # --- END INITIALIZATION ---

    # Start up model1 which fills av1 with data:
    av1 = model1( comm1, av1 )
    # print out av1 data:
    f1 = av1.exportRAttr("field1")
    f2 = av1.exportRAttr("field2")
    print "model 1 data:"
    print f1
    print f2
    # Rearrange data from model1 so that model2 can use it:
    # source, target, tag, sum
    Rranger.rearrange( av1, av2, 0, False )#??
    
    # pass data to model2 (which will print it out)
    model2( comm2, av2 )
    mpi.finalize()
    return

if __name__=="__main__":
    twoseq()
