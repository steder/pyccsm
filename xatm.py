"""
pyCPL6 Main Program

---

"""
import os, sys, time, traceback
import math
import Numeric
import mpi

try:
    assert os.getenv("LD_ASSUME_KERNEL")
    assert os.getenv("SIDL_DLL_PATH")
    assert os.getenv("PYTHONPATH")
except:
    print "Check your environment variables!"
#    rank = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
#    print "Printing Environement Variables on node %s" % (rank)
#    for key in os.environ.keys():
#        print "%s - %s: %s" % ( rank, key, os.environ[key] )

import cpl
import cpl.comm
import cpl.control
import cpl.contract
import cpl.infobuffer

class DeadModel:
    def __init__(self, name, nxg, nyg ):
        startTime = time.time()
        print "dead.py: starting @",time.ctime(),"!"
        if name in cpl.fields.component_names:
            local_comm = cpl.comm.init( name )
        else:
            raise ValueError,"Invalid Component Name %s"%(name)
        print "dead.py: OK!"
        self.name = name
        rank = mpi.comm_rank( mpi.MPI_COMM_WORLD )
        size = mpi.comm_size( mpi.MPI_COMM_WORLD )
        print "dead.py: mpi.MPI_COMM_WORLD size=%d,myrank=%d"%(size,rank)
        print "dead.py: Size of Xatmosphere Comm:",mpi.comm_size( local_comm )
        print "dead.py: Rank in Xatmosphere Comm:",mpi.comm_rank( local_comm )
        print "dead.py: world_pid:",cpl.comm.world_pid
        print "dead.py: mph_component_id:",cpl.comm.mph_component_id
        print "dead.py: component_pid:",cpl.comm.component_pid
        print "dead.py: component_npe:",cpl.comm.component_npe
        # How come HOSTNAME is undefined ?
        print "dead.py: Python Job running on %s, MPI Node %s"%(os.getenv("HOSTNAME"),rank)
        
        # fortran: call shr_msg_dirio('cpl')
        # Set flags for vector-friendly mapping    
        print "dead.py: read input namelist file"
        cpl.control.init()
        # cpl.control.readNamelist( "atm.stdin" )
            
        print "Initializing infobuffers..."
        self.infobuf = cpl.infobuffer.Infobuffer()
        """
        atm.stdin:
        atm                  !  myModelName (must match one in cpl_fields_mod)
        96                  !  i-direction global dimension
        48                  !  j-direction global dimension
        1                    !  decomp_type  1=1d-in-i, 2=1d-in-j, 3=2d, 4=segmented
        0                    !  num of pes for i (type 3 only)
        0                    !  length of segments (type 4 only)
        24                   !  ncpl  number of communications w/coupler per day
        0.0                  !  simul time proxy (secs): time between cpl comms
        1                    !  ignore_mph flag
        """
        self.infobuf.ibufSet("ncpl",24)
        self.infobuf.ibufSet("dead",1)
        self.infobuf.ibufSet("userest",1)
        self.infobuf.ibufSet("inimask",1)
        self.infobuf.ibufSet("gisize",nxg)
        self.infobuf.ibufSet("gjsize",nyg)
        self.infobuf.ibufSet("gsize",nxg * nyg)

        ###
        # Define Model Fields:
        ###
        # assuming Dead Atmosphere for now:
        self.fields_m2c_total = cpl.fields.a2c_total
        self.fields_m2c_list = cpl.fields.a2c_fields
        self.fields_c2m_total = cpl.fields.c2a_total
        self.fields_c2m_list = cpl.fields.c2a_fields

        ###
        # Decompose Model - This should become a method.
        # Possibly in the cpl package, since this
        # same decomposition code is used by both models
        # and couplers.
        ###
        # Assume 1D decomposition by latitude:
        # Default Decomposition:
        nx = nxg
        ny = nyg / cpl.comm.component_npe
        lsize = (nx*ny)
        
        self.infobuf.ibufSet("lsize",lsize )
        # Build Global Data Buffer:
        self.data_buffer = Numeric.zeros( (cpl.fields.grid_total, (nx*ny)), Numeric.Float32 )
        n=0
        ng = 0
        nglast = 0
        mype = cpl.comm.component_pid
        try:
            for j in xrange(1,ny+1):#j=1,ny                    ! local  j index
                for i in xrange(1,nx+1):#do i=1,nx                    ! local  i index
                    # 1D decomposition in the y dimension
                    ig = i                # global i index
                    jg = j + (mype*ny)      # global j index
                    nglast = ng
                    ng = (jg-1)*nxg + ig      #! global n index (vector index)
                    if False:
                        if( ng > (nglast+1) ):
                            print "*** ng = %s && nglast = %s ***"%(ng, nglast)
                    self.data_buffer[cpl.fields.grid_indices["lon"]][n] =         (ig-1)*360.0/(nxg)
                    self.data_buffer[cpl.fields.grid_indices["lat"]][n]  = -90.0 + (jg-1)*180.0/(nyg-1)
                    # ng - 1 adjusts from Fortran to C indices
                    self.data_buffer[cpl.fields.grid_indices["index"]][n] = ng 
                    n  = n+1                  # local  n index (vector index)
        except IndexError:
            print "lon: %s"%( cpl.fields.grid_indices["lon"])
            print "lat: %s"%( cpl.fields.grid_indices["lat"])
            print "index: %s"%( cpl.fields.grid_indices["index"])
            print "n: %s"%(n)
            raise
        #print "self.data_buffer['index']:"
        #print self.data_buffer[cpl.fields.grid_indices["index"]]
        # ATM Contract
        # First Parallel Calls
        print "dead.py: Contract Init:  establishes domains & routers (excluding land)"
        print "dead.py: thinks that the world communicator is:",cpl.comm.world_comm
        print "dead.py: thinks that the local communicator is:",cpl.comm.local_comm
        
        self.con_xm2c = cpl.contract.SendingContract( cpl.fields.atmname,
                                                 cpl.fields.cplname,
                                                 cpl.control.decomp_a,
                                                 self.fields_m2c_list,
                                                 self.infobuf,
                                                 self.data_buffer,
                                                 )
        self.con_xc2m = cpl.contract.SendingContract( cpl.fields.atmname,
                                                 cpl.fields.cplname,
                                                 cpl.control.decomp_a,
                                                 self.fields_c2m_list,
                                                 self.infobuf,
                                                 self.data_buffer,
                                                 )
                                                 
        # ICE Contract
        # OCN Contract
        
        # Begin LND Contract 
        # Special domain data Contract -- unique to land
        # End LND Contract
        # End of IF
        print "dead.py: Send domain info to land model? (optional)"
        if( not cpl.control.sendLndDom ):
            print " * lnd DOES NOT request optional domain data exchange"
        else:
            print " * lnd requests optional domain data exchange"
        # End If
        
        print "Receiving initial message from coupler with cdate, etc."
        self.infobuf.recv( )
        dbug = self.infobuf.ibufGet('infobug')

        
        print "Compute IC data:"
        for nf in xrange( 1, self.fields_m2c_total ):
            for j in xrange( 0, ny ):
                for i in xrange( 0, nx ):
                    n = (j * nx) + i #! local 1D index
                    try:
                        lon = self.data_buffer[cpl.fields.grid_indices["lon"]][n]
                        lat = self.data_buffer[cpl.fields.grid_indices["lat"]][n]
                    except IndexError:
                        print "n = %s, i = %s, j = %s, nx = %s"%(n,i,j,nx)
                        raise
                    tmp = ( (nf*100) * math.cos( math.pi * (lat/180.0)) * (math.sin( math.pi * (lon/180.0) ) -
                                                                           (cpl.comm.mph_component_id-1)*(math.pi / 3.0)) +
                            (cpl.comm.mph_component_id * 10.0) )
        print "Sending IC data to the coupler:"
        self.con_xm2c.send()
        
        print "dead.py: Init Done!"
        finishTime = time.time()
        print "dead.py: Elapsed Time:",(finishTime - startTime)
        sys.stdout.flush()

    def run(self):
        # Main Integration Loop: repeat until cpl says stop:
        print self.name,": beginning main integration loop..."
        iterations=0 # outer loop index
        ncpl = self.infobuf.ibufGet('ncpl')
        done = False
        while not done:
            #print "xatm: Top of outer loop:"
            iterations+=1
            # Beginning of NCPL Loop
            for n in xrange(ncpl):
                #print "%s: Starting hour %s..."%(self.name, n)
                # Get initial data from Coupler
                self.con_xc2m.recv()
                stopnow = self.con_xc2m.infobuffer.ibufGet('stopnow')
                if (stopnow!=0):
                    print "CPL SAYS STOP NOW"
                    done = True
                    break

                # Update my infobuffer:
                self.infobuf.ibufSet( 'cdate',
                                      self.con_xc2m.infobuffer.ibufGet('cdate') )
                self.infobuf.ibufSet( 'sec',
                                      self.con_xc2m.infobuffer.ibufGet('sec') )
                # Bunch of computations -- update con_xm2c
                # i'm skipping for now to focus on communications
                # Contract Send:
                self.con_xm2c.send()
                #print "%s: Finished hour %s!"%(self.name, n)
                
            # End of NCPL Loop
            if (done):
                break
        print self.name,": end of main integration loop..."
            
if __name__=="__main__":
    # initialize an atmosphere with the default decomposition on a 96x48 grid.
    model = DeadModel( "atm", 96, 48 )
    model.run()
    cpl.comm.finalize()
    print "xatm done!"
