"""
comm.py

Python Version of cpl_comm_mod.F90

The names of these constants have been changed to reflect that
they are defined inside of a python class and namespace information
does not have to be encoded in the variable names.

Hopefully the name changes are clear(er) and more easily understood.

Relies on F2py'd MPH Module: pymph

A call to 'init'

returns a local communicator for all the processors in this model 
component.  The 'init' call invokes a series of MPH routines to 
set all of the module's constants, and finishes by initializing MCT's
World.
"""

# Standard Python Modules
import os
import sys
import string

# Required External Modules
import mpi
#from _pymph import mphmod

# Coupler Library Modules
import fields

#mph_components, mph_global_proc_id, mph_local_proc_id, mph_total_components, mph_comp_id, mph_local_totprocs, mph_global_totprocs, mph_global_id, mph_comp_name, mph_global_world

from MCT import World

argv = ['']

# INIT WITHOUT MPH
def init( name ):
    print "setting up communicators, name =", name
    print name, ": Initializing comm groups via MPI_COMM_SPLIT..."
    ncomponents = 0

    myid, numprocs = mpi.init( len(sys.argv), sys.argv )

    global mct_world, world_pe0
    global world_comm, world_npe, world_pid, local_comm, component_npe
    global component_pid, mph_component_id, mph_component_id_atm 
    global mph_component_id_ice, mph_component_id_lnd, mph_component_id_ocn
    global mph_component_id_cpl, mph_component_id_rtm
    global world_pe0_cpl,world_pe0_atm, world_pe0_ice, world_pe0_lnd
    global world_pe0_rtm, world_pe0_ocn

    mph_component_id_cpl=1
    mph_component_id_atm=2
    mph_component_id_ocn=3
    mph_component_id_ice=4
    mph_component_id_lnd=5
    mph_component_id_rtm=5
    if name == "cpl":
        mph_component_id = 1
        ncomponents = 1
        local_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 0, 0 )
    elif name == "atm":
        mph_component_id = 2
        ncomponents = 1
        local_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 1, 0 )
    elif name == "ocn":
        mph_component_id = 3
        ncomponents = 1
        local_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 2, 0 )
    elif name == "ice":
        mph_component_id = 4
        ncomponents = 1
        local_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 3, 0 )
    elif name == "lnd":
        mph_component_id = 5
        ncomponents = 1
        local_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 4, 0 )
    else:
        print "(%s,%s): Invalid Name!" % ( myid, numprocs )

    print "%s (%s,%s): Comm group initialized!"%(name,myid,numprocs) 
    
    local_rank = mpi.comm_rank( local_comm )
    local_size = mpi.comm_size( local_comm )

    print "%s (%s,%s): creating root communicator!"%(name,myid,numprocs)
    sys.stdout.flush()
    if local_rank == 0:
        non_root_comm = None
        root_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 5, 0 )
        ncomponents = mpi.comm_size( root_comm )
        print "%s (%s,%s): joined root communicator!"%(name,myid,numprocs)
        sys.stdout.flush()
    else:
        root_comm = None
        non_root_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 6, 0 )
        print "%s (%s,%s): Joined non-root communicator!"%(name,myid,numprocs)
        sys.stdout.flush()
        ncomponents = 0
    print "%s (%s,%s): Distributing root communicator!"%(name,myid,numprocs) 
    ncomponents = mpi.bcast( ncomponents, 1, mpi.MPI_INT, 0, local_comm )
    print "%s (%s,%s): Distributed root communicator!"%(name,myid,numprocs) 
    # Get total number of components and distribute to every node:
    #ncomponents = mpi.allreduce( ncomponents, 1, mpi.MPI_INT, mpi.MPI_SUM, root_comm )
    ncomponents = ncomponents[0]
    print "%s(%s,%s): ncomponents = %s"%(name, myid, numprocs, ncomponents )
    # since root_comm has already been broadcast we shouldn't have to
    # broadcast ncomponents.
    #ncomponents = mpi.bcast( ncomponents, 1, mpi.MPI_INT, 0, local_comm )
    #ncomponents = ncomponents[0]
    sys.stdout.flush()
    print "%s(%s,%s): 1 of %s components!"%(name, myid,numprocs, ncomponents)
    # setup component variables:
    component_pid = local_rank
    component_npe = local_size
    sys.stdout.flush()
    print "%s(%s,%s): setting remaining Component values..."%(name,myid,numprocs)
    sys.stdout.flush()
    # Determin CID's and component pe0's for all components.
    world_pe0_cpl = 0
    world_pe0_atm = 0
    world_pe0_ocn = 0
    world_pe0_ice = 0
    world_pe0_lnd = 0
    world_pe0 = 0 # is always 0, doesn't need to be reassigned.
    # We now assign these values and do an allreduce in MPI_COMM_WORLD
    # to distribute them to every node.
    sys.stdout.flush()
    # myid, numprocs are from the mpi.init call
    world_size = numprocs #mpi.comm_size( mpi.MPI_COMM_WORLD )
    world_rank = myid #mpi.comm_rank( mpi.MPI_COMM_WORLD )

    # Setup local/global comm info, and component id's
    world_comm = mpi.MPI_COMM_WORLD
    world_pid = world_rank #mpi.comm_rank( mpi.MPI_COMM_WORLD )
    world_npe = world_size #mpi.comm_size( mpi.MPI_COMM_WORLD )

    if name == "cpl":
        print "(%s,%s) is a cpl node."%(world_rank,world_size)
        if local_rank == 0:
            world_pe0_cpl = world_rank
    if name == "atm":
        print "(%s,%s) is a atm node."%(world_rank,world_size)
        if local_rank == 0:
            world_pe0_atm = world_rank
    if name == "ocn":
        print "(%s,%s) is a ocn node."%(world_rank,world_size)
        if local_rank == 0:
            world_pe0_ocn = world_rank
    if name == "ice":
        print "(%s,%s) is a ice node."%(world_rank,world_size)
        if local_rank == 0:
            world_pe0_ice = world_rank
    if name == "lnd":
        print "(%s,%s) is a lnd node."%(world_rank,world_size)
        if local_rank == 0:
            world_pe0_lnd = world_rank
    world_pe0_cpl = mpi.allreduce( world_pe0_cpl, 1, mpi.MPI_INT,
                                       mpi.MPI_SUM, mpi.MPI_COMM_WORLD )
    world_pe0_atm = mpi.allreduce( world_pe0_atm, 1, mpi.MPI_INT,
                                       mpi.MPI_SUM, mpi.MPI_COMM_WORLD )
    world_pe0_ocn = mpi.allreduce( world_pe0_ocn, 1, mpi.MPI_INT, 
                                       mpi.MPI_SUM, mpi.MPI_COMM_WORLD )
    world_pe0_ice = mpi.allreduce( world_pe0_ice, 1, mpi.MPI_INT, 
                                       mpi.MPI_SUM, mpi.MPI_COMM_WORLD )
    world_pe0_lnd = mpi.allreduce( world_pe0_lnd, 1, mpi.MPI_INT, 
                                       mpi.MPI_SUM, mpi.MPI_COMM_WORLD )
    world_pe0_cpl = world_pe0_cpl[0]
    world_pe0_atm = world_pe0_atm[0]
    world_pe0_ocn = world_pe0_ocn[0]
    world_pe0_ice = world_pe0_ice[0]
    world_pe0_lnd = world_pe0_lnd[0]
    

    if False:
        print "world_pe0_cpl = %s"%(world_pe0_cpl)
        print "world_pe0_atm = %s"%(world_pe0_atm)
        print "world_pe0_ocn = %s"%(world_pe0_ocn)
        print "world_pe0_ice = %s"%(world_pe0_ice)
        print "world_pe0_lnd = %s"%(world_pe0_lnd)
    
    # Initialize MCT:
    print name,": Setting up MCT...",
    mct_world = World.World()
    print "mct_world.initd0(%s,%s,%s,%s)"%(ncomponents, world_comm, local_comm, mph_component_id )
    mct_world.initd0( ncomponents, world_comm, local_comm, mph_component_id )
    print "Comm State:"
    # Document comm groups, pe0's, mph component ids
    print "comm world : comm, npe, pid -", world_comm, world_npe, world_pid
    print "comm component : comm,npe,pid -",local_comm, component_npe, component_pid
    print "comm world pe0 : atm,ice,lnd,ocn,cpl -",world_pe0_atm,world_pe0_ice,
    print world_pe0_lnd, world_pe0_ocn,world_pe0_cpl
    print "component ids : atm,ice,lnd,ocn,cpl,me -",mph_component_id_atm,
    print mph_component_id_ice,mph_component_id_lnd,mph_component_id_ocn,
    print mph_component_id_cpl,mph_component_id
    print name,": comm.init: finished!"
    sys.stdout.flush()
    return local_comm

def finalize():
    mpi.finalize()
    mct_world.clean()
    return
    

# def init( name , ignore_mph=False ):
#     """
#     local_communicator = comm.init( name )
    
#     Where 'name' is usually one of the following
#     ['atm','ice','lnd','ocn','cpl']
    
#     Returns a MPI communicator object containing all the processors
#     with the same processor id(name) as this component.
#     """
    
#     print "pyCPL: sys.argv =",sys.argv, argv
#     mpi.init( len(argv), argv )
#     print "pyCPL: MPI WORLD VALUES:"
#     print "pyCPL: rank = %d"%( mpi.mpi_comm_rank( mpi.MPI_COMM_WORLD ) )
#     print "pyCPL: size = %d"%( mpi.mpi_comm_size( mpi.MPI_COMM_WORLD ) )
    
#     print "pyCPL: setting up communicators, name =", name
    
#     # Make sure assignments go to the global variables above:
#     global mct_world, world_pe0
#     global world_comm, world_npe, world_pid, local_comm, component_npe
#     global component_pid, mph_component_id, mph_component_id_atm 
#     global mph_component_id_ice, mph_component_id_lnd, mph_component_id_ocn
#     global mph_component_id_cpl, mph_component_id_rtm
#     global world_pe0_cpl,world_pe0_atm, world_pe0_ice, world_pe0_lnd
#     global world_pe0_rtm, world_pe0_ocn
#     print "pyCPL: Entering Fortran 'mph_components' call..."
#     comm = mphmod.mph_components( name )
#     print "pyCPL: Leaving 'mph_components' call..."
#     print "pyCPL: Testing comm object returned from mph_components call..."
#     """
#     comm = mphmod.mph_local_world seems to work when the comm returned
#     from mphmod.mph_components seems to be invalid.

#     I'm thinking that this occurs because the conversion done py F2PY
#     is not correct for this mpi type.
#     """
#     temp_comm = mphmod.mph_local_world
#     print "pyCPL: comm = %s"%(comm)
# #rank = mpi.mpi_comm_rank( comm )
# #size = mpi.mpi_comm_size( comm )
# #print "pyCPL: MPH thinks that this component is rank %s, of %s"%(rank,size)
#     # Comm is a fortran object:  
#     # let's use mpi.communicator(an undocumented pyMPI feature)
#     # to generate a python communicator.
#     #comm = mpi.communicator(comm)
#     print "pyCPL: setting global variables..."
#     world_comm = mphmod.mph_global_world
#     #world_comm = mpi.communicator(world_comm)
    
#     world_pid = mphmod.mph_global_proc_id()
#     world_npe = mphmod.mph_global_totprocs()
#     world_pe0 = mphmod.mph_global_id(name, 0)
#     print "pyCPL: world_pid:",world_pid
#     print "pyCPL: world_npe:",world_npe
#     print "pyCPL: world_pe0:",world_pe0
#     local_comm = comm
#     print "pyCPL: local_comm:",local_comm
#     mph_component_id = mphmod.mph_comp_id(name)
#     print "pyCPL: mph_component_id:",mph_component_id
#     component_pid = mphmod.mph_local_proc_id(mph_component_id)
#     print "pyCPL: component_pid:",component_pid
#     component_npe = mphmod.mph_local_totprocs(mph_component_id)
#     print "pyCPL: component_npe:",component_npe
    
#     # Determine MPH Component ID's and comm_world pe0's for all components
#     n = mphmod.mph_total_components()
#     print "pyCPL: total components:",n
#     print "pyCPL: Determining MPH Component ID's for %d components..."%(n)
#     for cid in xrange( 1, n+1 ): # Component ID's start at 1, not 0
#         name = mphmod.mph_comp_name(cid)
#         name = name.strip() # remove leading/trailing whitespace
#         if ( name == fields.atmname ):
#             mph_component_id_atm = cid
#             world_pe0_atm = mphmod.mph_global_id( fields.atmname, 0 )
#             print "pyCPL: world_pe0_atm:",world_pe0_atm
#         elif ( name == fields.icename ):
#             mph_component_id_ice = cid
#             world_pe0_ice = mphmod.mph_global_id( fields.icename, 0 )
#             print "pyCPL: world_pe0_ice:",world_pe0_ice
#         elif ( name == fields.lndname ):
#             mph_component_id_lnd = cid
#             mph_component_id_rtm = cid
#             world_pe0_lnd = mphmod.mph_global_id( fields.lndname, 0 )
#             world_pe0_rtm = mphmod.mph_global_id( fields.lndname, 0 )
#             print "pyCPL: world_pe0_lnd:",world_pe0_lnd
#             print "pyCPL: world_pe0_rtm:",world_pe0_rtm
#         elif ( name == fields.ocnname ):
#             mph_component_id_ocn = cid
#             world_pe0_ocn = mphmod.mph_global_id( fields.ocnname, 0 )
#             print "pyCPL: world_pe0_ocn:",world_pe0_ocn
#         elif ( name == fields.cplname ):
#             mph_component_id_cpl = cid
#             world_pe0_cpl = mphmod.mph_global_id( fields.cplname, 0 )
#             print "pyCPL: world_pe0_cpl:",world_pe0_cpl
#         else:
#             #print "pyCPL: [%s]"%(name)
#             print "pyCPL: mph_component_name error:",name," is not a supported name."
#             # Raise Exception:  MPH_COMPONENT_NAME_ERROR
        
#     """
#     The following is a workaround for invalid values being
#     returned from MPH:
#     """
#     if ignore_mph:
#         print "pyCPL:  IGNORING MPH!"
#         local_comm = mpi.mpi_comm_split( mpi.MPI_COMM_WORLD, mph_component_id, 0)
#         world_pe0_cpl = 1
#         world_pe0_atm = 0
#         print "pyCPL:  SUCCESSFULLY IGNORED MPH!"

#     # Initialize MCT:
#     print "pyCPL: Setting up MCT...",
#     mct_world = World.World()
    
#     try:
#         mpi_world_size = mpi.mpi_comm_size( mpi.MPI_COMM_WORLD )
#     except:
#         print "pyCPL: mpi_world_size=",mpi_world_size

#     mct_world.initd0( 2, world_comm, local_comm, mph_component_id )
#     print " success!"
# #    if( mpi_world_size >= mphmod.mph_total_components() ):
# #        mct_world.initd0( mphmod.mph_total_components(), world_comm, local_comm, mph_component_id )
# #    elif( mpi_world_size == 1 ):
# #        print "pyCPL: ** Initializing MCT world on %s processes" %(mpi_world_size)
# #        print "pyCPL:   This could be quite dangerous ! **"
# #        print "pyCPL: mpi.WORLD.size=",mpi_world_size
# #        print "pyCPL: world_comm",world_comm
# #        print "pyCPL: local_comm",local_comm
# #        print "pyCPL: mph_component_id",mph_component_id
# #        mct_world.initd0( 1, world_comm, local_comm, 1 )
# #    else:
# #        print "pyCPL: ** Warning:  MCT World was NOT initialzed because there are fewer processors then components (no way to initialize without deadlocking) **"
#     print "pyCPL: Comm:  Initialized!"
#     return local_comm


    
