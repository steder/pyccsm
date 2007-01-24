"""
contract.py

Python version of cpl_contract_mod.F90

Defines a contract class
"""

# Builtin Modules
import sys,time

# External Modules
import mpi
import Numeric
from MCT import Router,GlobalSegMap

# Internal Modules
import attributevector
import bundle
import comm
import control
import decomposition
import domain
import error
import fields
import infobuffer

from debug import newDebugPrint
# Make these local to the module:
ContractError = error.ContractError
DEBUG = control.infodbug

debugPrint = newDebugPrint(False)
statusPrint = newDebugPrint(False)

class Contract:
    def __init__(self, my_name, other_name ):
        """
        Initializes the core of a contract, setting up communicators,
        component and processor ids.
        """
        self.name = my_name
        self.other_name = other_name
        self.__initialized_gsmap = False
        self.__initialized_router = False
        statusPrint("Contract.__init__:Creating a %s to %s contract object @ %s"%(my_name,other_name,time.ctime()))
        self.infobuffer = infobuffer.Infobuffer()
        self.m2c_fields = []
        self.c2m_fields = []
        self.domain = domain.Domain()
        self.router = Router.Router()
        #self.gsMap = GlobalSegMap.GlobalSegMap()
        # Handle Basic Initialization
        self.comm = comm.local_comm
        self.pid = comm.component_pid
        statusPrint("Contract.__init__:%s component id = %s"%(my_name, comm.mph_component_id))
        self.cid = comm.mph_component_id
        if ( my_name == fields.cplname ):
            self.pe0 = comm.world_pe0_cpl
        elif ( my_name == fields.atmname ):
            self.pe0 = comm.world_pe0_atm
        elif ( my_name == fields.ocnname ):
            self.pe0 = comm.world_pe0_ocn
        elif ( my_name == fields.icename ):
            self.pe0 = comm.world_pe0_ice
        elif ( my_name == fields.lndname ):
            self.pe0 = comm.world_pe0_lnd
        else:
            raise AttributeError,"contract.py:  My name is invalid!"
        self.tag = 1001
        statusPrint("Contract.__init__:component ids:")
        components = ""
        try:
            components += "%s:%d\n"%("cpl",comm.mph_component_id_cpl)
            components += "%s:%d\n"%("atm",comm.mph_component_id_atm)
            components += "%s:%d\n"%("lnd",comm.mph_component_id_lnd)
            components += "%s:%d\n"%("ocn",comm.mph_component_id_ocn)
            components += "%s:%d\n"%("ice",comm.mph_component_id_ice)
            #   components += "%s:%d\n"%("rtm",comm.mph_component_id_rtm)
            statusPrint(components)
        except:
            statusPrint(components)
            statusPrint("Contract.__init__:Not all components available")
        if (other_name == fields.atmname):
            cid_other = comm.mph_component_id_atm
            pid_other = comm.world_pe0_atm
            suffix = "a"
        elif( other_name == fields.lndname ):
            cid_other = comm.mph_component_id_lnd
            pid_other = comm.world_pe0_lnd
            suffix = "l"
        elif( other_name == fields.rtmname ):
            cid_other = comm.mph_component_id_lnd
            pid_other = comm.world_pe0_lnd
            suffix = "l"
        elif( other_name == fields.ocnname ):
            cid_other = comm.mph_component_id_ocn
            pid_other = comm.world_pe0_ocn
            suffix = "o"
        elif( other_name == fields.icename ):
            cid_other = comm.mph_component_id_ice
            pid_other = comm.world_pe0_ice
            suffix = "i"
        elif( other_name == fields.cplname ):
            cid_other = comm.mph_component_id_cpl
            pid_other = comm.world_pe0_cpl
            if ( my_name == fields.atmname ):
                suffix = "a"
            elif( my_name == fields.icename ):
                suffix = "i"
            elif( my_name == fields.lndname ):
                suffix = "l"
            elif( my_name == fields.ocnname ):
                suffix = "o"
            else:
                print "Contract.__init__:ERROR: this should never happen!"
                raise ContractError, "Contract.__init__:: Unrecognized my_name = %s"%(my_name)
        else:
            print "Contract.__init__:ERROR:  this should never happen!"
            raise ContractError, "Contract.__init__:: Unrecognized other_name = %s"%(other_name)
        statusPrint("Contract.__init__:%s component id = %s @ %s" %(other_name, cid_other,time.ctime()))
        self.cid_other = cid_other
        self.pid_other = pid_other
        self.suffix = suffix
        return

    def __del__(self):
        if( self.__initialized_gsmap ):
            self.gsMap.clean()
        if( self.__initialized_router ):
            self.router.clean()
        return

    def initAV( self, fields ):
        """
        """
        #print "initAV: Entering initAV"
        lsize = self.domain.lgrid.size()
        #print "initAV: Calling Initialize"
        self.av = attributevector.AttributeVector( [], fields, lsize )
        #print "initAV: returning"
        return

    def execute(self):
        """
        All sending and receiving contracts define execute, but the 
        two types of contracts define it differently.
        """
        return

    def send(self):
        """
        Send data/message to component.
        
        Send message containing integer parameters
        Send message containing time-variant data such as state 
        and forcing fields
        """
        statusPrint(self.name,": sending contract to",self.other_name,"...")
        self.tag = 1003
        #--- send info-buffer data ---
        if (self.pid == 0):
            debugPrint(self.name,": sending infobuffer to",self.other_name,"...")
            self.infobuffer.send( self.pid_other, self.tag, comm.world_comm )
            debugPrint(self.name,": infobuffer sent!")
        #--- send bundle data ---
        debugPrint(self.name,": PyMCT Send...")
        #print "Before PyMCT Send:  Latitudes:"
        #lat,lat_size = self.av.av.exportRAttr("lat")
        #print lat
        self.av.send( self.router, 600 )
        #print self.name,": attributevector.py:send"
        #self.av.send(self.router)
        
        debugPrint(self.name,":",self.other_name,"contract sent!")
        return

    def recv( self ):
        statusPrint(self.name,": receiving contract from",self.other_name,"...")
        self.tag = 1003
        # Initialize bundle to 0:
        # This should probably be default behavior of Bundle() constructor
        # call cpl_bundle_zero(contract%bundle)
        if (self.pid == 0):
            debugPrint("\t* Receiving Infobuffer")
            debugPrint("\t  -self.pid:",self.pid)
            debugPrint("\t  -self.pid_other:",self.pid_other,"\n",
                       "\t  -self.tag:",self.tag,"\n",
                       "\t  -comm.world_comm:",comm.world_comm)
            self.infobuffer.recv( self.pid_other, self.tag, comm.world_comm )
        debugPrint("\t* Broadcasting Infobuffer")
        self.infobuffer.bcast( self.comm, 0 )
        #--- recv bundle data ---
        # print "\t* receiving attribute vector"
        self.av.zero()
        debugPrint("Receiving Contract Bundle:","\n",
                   "*size:",self.av.size(),"\n",
                   "*fields:","\n",self.av)
        self.av.recv( self.router, 600, False )
        debugPrint(self.name,":", self.other_name, "contract received!")
        return 

class ReceivingContract( Contract ):
    def __init__(self, my_name, other_name, decomp, myfields):
        """
        Initializes the contract as a receiving contract
        """
        statusPrint("contract.py: %s->%s Receiving Contract Initializing @ %s..."%(my_name,other_name,time.ctime()))
        Contract.__init__(self, my_name, other_name )

        self.m2c_fields = myfields
        self.c2m_fields = myfields
        
        # Receive infobuffer
#print "contract.py: Contracts think that the world communicator is:",comm.world_comm
#print "contract.py: Contracts think that the local communicator is:",comm.local_comm
        if( self.pid == 0 and ( mpi.comm_size( mpi.MPI_COMM_WORLD ) > 1) ):
            # print "contract.py: %s->%s Receiving infobuffers..."%(my_name,other_name)
            self.infobuffer.recv( self.pid_other, self.tag, comm.world_comm )
            gsize = self.infobuffer.ibufGet("gsize")
            gisize = self.infobuffer.ibufGet("gisize")
            gjsize = self.infobuffer.ibufGet("gjsize")
            # print "Received: gsize=%s\tgisize=%s\tgjsize=%s"%(gsize,gisize,gjsize)
            gsize = self.infobuffer.ibufGet("gsize")
            gisize = self.infobuffer.ibufGet("gisize")
            gjsize = self.infobuffer.ibufGet("gjsize")
            # print "Broadcasted: gsize=%s\tgisize=%s\tgjsize=%s"%(gsize,gisize,gjsize)
        # print "contract.py: %s->%s Broadcasting infobuffers..."%(my_name,other_name)
        self.infobuffer.bcast( comm.local_comm, root=0 )
        """
        print "contract.py: %s->%s Debugging Dump of Infobuffer Values:"%(my_name,other_name)
        ifields,rfields = self.infobuffer.getFields()
        for ifield in ifields:
            print "contract.py: %s: %s"%(ifield,self.infobuffer.ibufGet(ifield))
        for rfield in rfields:
            print "contract.py: %s: %s"%(rfield,self.infobuffer.rbufGet(rfield))
        print "contract.py: %s->%s Completed Infobuffer Dump!"%(my_name,other_name)
        """
        # Get Local index values:
        gsize = self.infobuffer.ibufGet("gsize")
        gisize = self.infobuffer.ibufGet("gisize")
        gjsize = self.infobuffer.ibufGet("gjsize")

        # print "Adjusting gisize, gjsize:(HACK!)"
        # print "before: gsize=%s\tgisize=%s\tgjsize=%s"%(gsize,gisize,gjsize)
        if (gsize != 0):
            if( (gisize==0) and (gjsize!=0) ):
                gisize = gsize / gjsize
                self.infobuffer.ibufSet("gisize",gisize)
            if( (gisize!=0) and (gjsize==0) ):
                gjsize = gsize / gisize
                self.infobuffer.ibufSet("gjsize",gjsize)
        # print "after: gsize=%s\tgisize=%s\tgjsize=%s"%(gsize,gisize,gjsize)

        debugPrint( "contract.py: %s->%s Infobuffers received!"%(my_name,other_name))
        debugPrint( "contract.py: %s->%s Decomposing domain..."%(my_name,other_name))
        if ( decomp != None and type(decomp)==type(0) ):
            lsize,indx = decomposition.decompose( decomp, gisize, gjsize, comm.component_pid, comm.component_npe )
        else:
            lsize,indx = decomposition.decompose( 1, gisize, gjsize, comm.component_pid, comm.component_npe )
        # print "contract.py: %s->%s domain decomposed!"%(my_name, other_name)
        # Initialize gsMap from index data:
#print "contract.py: %s-%s contract Skipping GSMap Init:"%(my_name,other_name)
        # print "contract.py: %s->%s - contract GlobalSegMap Init:"%(my_name,other_name)
        #print "contract.py: %s-%s - initGsmap( %s, %s )"%(my_name,other_name,lsize,indx)

        #self.initGsmap( lsize, indx )
        self.domain.initGsmap( lsize, indx )

        # print "contract.py: %s->%s domain setup:"%(my_name,other_name)
        # Init Contracts, Setup lgrid, receive lgrid
        # print "contract.py: %s->%s: gsize = %d, gisize = %d, gjsize = %d"%(my_name,other_name,gsize,gisize, gjsize)
        self.domain.initGrid(other_name.strip() + " contract domain",
                                self.suffix.strip(), 
                                gsize, 
                                gisize, 
                                gjsize,
                                [],
                                fields.grid_fields,
                                lsize)
        
        # Receive lGrid data( requires a router )
        statusPrint("contract.py: %s->%s router init:"%(my_name,other_name))
        # print "contract.py: cid =",self.cid
        # print "contract.py: cid_other =",self.cid_other

        # print "contract.py: router.init( %s, %s, %s )" %(self.cid_other, self.gsMap, self.comm )
        self.router.init( self.cid_other, self.domain.gsMap, self.comm )       
        self.__initialized_router = True
        # print "contract.py: %s->%s contract receiving lgrid(domain) data:"%(my_name,other_name)
        sys.stdout.flush()
        #self.domain.lgrid.recv( self.router )
        # 600 is a magic number from MCT.
        self.domain.lgrid.recv( self.router,600 )# self.tag ) 
        # Write out some debug/sanity check info:

        self.domain.info()
        # print "Allocate and Initialize Bundle/AttrVect"
        self.initAV( self.c2m_fields )
        statusPrint("contract.py: %s->%s Contract successfully initialized @ %s!"%(my_name,other_name,time.ctime()))
        return
            
    def execute( self ):
        return self.recv( )
        

class SendingContract(Contract):
    def __init__(self, my_name, other_name, decomp, myfields, infobuffer, buf ):
        """
        Initializes this contract as a "sending" contract
        """
        ### Pull in initialization for all contracts
        Contract.__init__(self, my_name, other_name)
        # Handle Arguments:
        self.m2c_fields = myfields

        self.infobuffer = infobuffer
        # print "contract.py: SendingContract - get/set local index values"
        
        ### Real Sending Contract Work:
        gsize = self.infobuffer.ibufGet("gsize")
        gisize = self.infobuffer.ibufGet("gisize")
        gjsize = self.infobuffer.ibufGet("gjsize")
        # Send Infobuffer
        if (self.pid == 0):
            #print "before send: self.pid_other=%s\tself.tag=%s\tcomm.world_comm=%s"%(self.pid_other, self.tag, comm.world_comm )
            debugPrint("contract.py: Sending Infobuffer to %s from %s with tag %s."%(self.pid_other, self.pid,self.tag))
            self.infobuffer.send( self.pid_other, self.tag, comm.world_comm )

        # Get/set local index values
        gsize = self.infobuffer.ibufGet("gsize")
        gisize = self.infobuffer.ibufGet("gisize")
        gjsize = self.infobuffer.ibufGet("gjsize")
        # print "contract.py: SendingContract - gsize,gisize,gjsize = %s,%s,%s"%(gsize,gisize,gjsize)

        lsize = self.infobuffer.ibufGet( "lsize" )
        # if( buf == None ):
        #    print "contract.py: SendingContract - Allocating Buffer..."
        #    buf = Numeric.zeros( (fields.grid_total, gisize*gjsize), Numeric.Float32 )
        indx = buf[ fields.grid_indices["index"] ]

        # print "contract.py: SendingContract - lsize = %s" %(lsize)
        
        # print "contract.py: SendingContract - defining indx"
        self.domain.initGsmap( lsize, indx )
            
        # RECV side does this:
        #    def decompose( self, decomp, gi, gj, myid, npes ):
        #        if ( decomp != None and type(decomp)==type(0) ):
        #            lsize,indx = self.decompose( decomp, gisize, gjsize, comm.component_pid, comm.component_npe )
        #        else:
        #            lsize,indx = self.decompose( 1, gisize, gjsize, comm.component_pid, comm.component_npe )

        
        #print "contract.py: SendingContract - buf=%s\nfields.grid_index=%s\n"%(buf, fields.grid_index)
        
        #print "contract.py: SendingContract - indx = %s" %(indx)
        statusPrint("contract.py: SendingContract %s->%s contract GlobalSegMap Init:"%(my_name,other_name))
        

        # Initialize contracts, setup lGrid
        # On the send side lGrid exists
        # On the recv side, need to recv lGrid
        self.domain.initGrid(other_name.strip() + " contract domain",
                                self.suffix.strip(),
                                self.infobuffer.ibufGet("gsize"),
                                self.infobuffer.ibufGet("gisize"),
                                self.infobuffer.ibufGet("gjsize"),
                                "",
                                fields.grid_fields,
                                lsize)
        # print "contract.py: domain grid initialized!"
        # These fields should all be defined in the fields module
        debugPrint("buf[fields.grid_lon] (",len(buf[fields.grid_lon]),"):",
                   buf[fields.grid_lon][:10],"\n",
                   buf[fields.grid_lon][-10:])
        self.domain.lgrid.importRAttr( "lon", buf[fields.grid_lon] )
        lon = self.domain.lgrid.exportRAttr("lon")
        debugPrint("lon (",len(lon),"):",
                   lon[:10],"\n",
                   lon[-10:])
        # print "contract.py: SendingContract - Imported 'lon' grid!"
        self.domain.lgrid.importRAttr( "lat", buf[fields.grid_lat] )
        # print "contract.py: SendingContract - Imported 'lat' grid!"
        self.domain.lgrid.importRAttr( "mask", buf[fields.grid_mask] )
        # print "contract.py: SendingContract - Imported 'mask' grid!"
        self.domain.lgrid.importRAttr( "area", buf[fields.grid_area] )
        # print "contract.py: SendingContract - Imported 'area' grid!"
        self.domain.lgrid.importRAttr( "index", buf[fields.grid_index] )
        # print "contract.py: SendingContract - Imported 'index' grid!"
        #self.domain.lGrid.importRAttr( "mpad", buf[fields.grid_mpad] )
        #self.domain.lGrid.importRAttr( "mcps", buf[fields.grid_mcps] )
        tmp = Numeric.zeros( lsize, Numeric.Float64 )
        tmp[0] = comm.world_pid
        self.domain.lgrid.importRAttr( "pid", tmp )
        # print "contract.py: domain grid defined!"
            
            
        # Send lGrid (requires a router)
        # Initialize router:
        statusPrint("contract.py: router.init( %s, %s, %s )"%(self.cid_other, self.gsMap, self.comm ))
        self.router.init( self.cid_other, self.gsMap, self.comm )
        self.__initialized_router = True
        # print "contract.py: sending domain..."
        sys.stdout.flush()
        # self.tag replaced with 600 (magic number from MCT )
        self.domain.lgrid.send( self.router, 600 )
        # print "contract.py: domain sent!"
        # Write some debug/sanity-check info to stdout:
        self.domain.info()
        # print "Allocate and Initialize Bundle/AttrVect"
        self.initAV( self.m2c_fields )
        statusPrint("contract.py: SendingContract __init__ complete @ %s!"%(time.ctime()))
        return
        
    def execute( self ):
        return self.send()
        
