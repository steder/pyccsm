#! /usr/bin/env python
# twocmp-con.py
"""
Mike Steder

A simple MCT example using Two Components running Sequentially
"""

# Check all modules are available
#import modulecheck
# All modules are available:
import sys,math
import mpi
import MCT,MCT.World,MCT.AttrVect,MCT.GlobalSegMap,MCT.Router
import MCT.Rearranger
import Numeric

# CONSTANTS:
NGX = 6 # Number of points in the x-direction
NGY = 4 # Number of points in the y-direction

def model1(numcomps, id, otherid, comm):
    """
    These example 'models' simply fill the attribute vectors with
    data.
    """
    print "* model1 starting! *"
    world = MCT.World.World()
    mysize = mpi.comm_size(comm)
    myrank = mpi.comm_rank(comm)
    print "model1 size:",mysize
    if(mysize > 1):
        print "This script assumes only 2 processors, 1 per component model."
    # Init World, first MCT call
    world.initd0(numcomps,mpi.MPI_COMM_WORLD,comm,id)
    # Setup some bookkeeping variables for a Global Segment Map: start,length
    npoints = NGX * NGY
    localsize = npoints / mysize
    start = Numeric.zeros((mysize))
    start[0] = (localsize * myrank) + 1
    length = Numeric.zeros((mysize))
    length[0] = localsize
    # Describe decomposition with an MCT GSmap type:
    GSMap = MCT.GlobalSegMap.GlobalSegMap()
    # 0 = Root of processor of 'comm'
    print "model1: initializing GSMap..."
    #print "model1: start(%s),length(%s),root(0),comm(%d),id(%d)"%(start, length, comm, id)
    GSMap.initd0(start, length, 0, comm, id)
    print "model1(%d) - #global segments(%d) - start indices(%s)"%(myrank, GSMap.ngsegs(), start)
    # Initialize an AV:
    av = MCT.AttrVect.AttrVect()
    av.init("","field1:field2",GSMap.lsize(comm))
    avsize = av.lsize()
    print "model1: av size",avsize
    # Fill AV with some data:
    data = Numeric.array(range(npoints),Numeric.Float64)
    av.importRAttr("field1",data)
    data = [ x/npoints for x in data ]
    av.importRAttr("field2",data)
    # Initialize a router:
    router = MCT.Router.Router()
    router.init(otherid, GSMap, comm)
    # Print Data:
    for f in ["field1","field2"]:
        data = av.exportRAttr(f)
        print data
    # Send Data to Model2:
    av.send(router,7)#7 is just a tag value, it must match the recv
    print "* model1 finished! *"
    return
    
def model2(numcomps,id,otherid,comm):
    """
    This example 'model' simply fills the attribute vector with data.
    One small difference between this and model1 is that this model
    prints it's own data.
    """
    print "* model2 starting! *"
    world = MCT.World.World()
    mysize = mpi.comm_size(comm)
    myrank = mpi.comm_rank(comm)
    print "model2 size:",mysize
    if(mysize > 1):
        print "This script assumes only 2 processors, 1 per component model."
    # Init World, first MCT call
    world.initd0(numcomps,mpi.MPI_COMM_WORLD,comm,id)
    # Setup some bookkeeping variables for a Global Segment Map: start,length
    npoints = NGX * NGY
    localsize = npoints / mysize
    start = Numeric.zeros((mysize))
    start[0] = (localsize * myrank) + 1
    length = Numeric.zeros((mysize))
    length[0] = localsize
    # Describe decomposition with an MCT GSmap type:
    GSMap = MCT.GlobalSegMap.GlobalSegMap()
    # 0 = Root of processor of 'comm'
    print "model2: initializing GSMap..."
    #print "model2: start(%s),length(%s),root(0),comm(%d),id(%d)"%(start, length,comm, id)
    GSMap.initd0(start, length, 0, comm, id)
    print "model2(%d) - #global segments(%d) - start indices(%s)"%(myrank, GSMap.ngsegs(), start)
    # Initialize an Attribute Vector:
    av = MCT.AttrVect.AttrVect()
    av.init("","field1:field2",GSMap.lsize(comm))
    avsize = av.lsize()
    print "model2: av size",avsize
    # ZERO out the AV
    zero = Numeric.zeros((npoints),Numeric.Float64)
    for f in ["field1","field2"]:
        av.importRAttr(f,zero)

    # Initialize a router:
    router = MCT.Router.Router()
    router.init(otherid, GSMap, comm)
    # Export and Print AV Data before recv:
    for f in ["field1","field2"]:
        data = av.exportRAttr(f)
        print data
    av.recv(router,7,False)#(router, tag, bool Sum)
    # Export and Print AV Data after recv:
    for f in ["field1","field2"]:
        data = av.exportRAttr(f)
        print data
    print "* model2 finished! *"
    return
    
def main():
    print "Entering twocmp-seq.py main..."
    myrank,mysize = mpi.init()
    if(mysize != 2):
        print "twocmp-seq.py requires 2 processors to run!"
        sys.exit()
    
    ## Initialize Models:
    # Initialize MCT World
    print "Initializing MCT World..."
    # We use initd1 because initd0 doesn't take a myids array.
    componentID = (myrank%2) + 1
    otherID = (componentID%2) + 1
    # Create communicators for each model
    print "componentID:",componentID
    subcomm = mpi.comm_split(mpi.MPI_COMM_WORLD,componentID)
    if( componentID == 1):
        model1(2,componentID,otherID,subcomm)
    elif(componentID== 2):
        model2(2,componentID,otherID,subcomm)
    else:
        print "Unknown Component ID!"
    
    print "Example has successfully finished!"
    mpi.finalize()
    return

if __name__=="__main__":
    main()

    
