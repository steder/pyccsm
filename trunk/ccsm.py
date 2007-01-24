#!/usr/bin/env python2.4

"""
pyCPL6 Main Program

---

"""
# Standard Modules:
import os, sys, time, traceback
import pickle

# Required external modules:
import Numeric
import pycdf
# import numpy as Numeric

# Optional External Modules:
try:
    import matplotlib
    import pylab
    MovieMode = True
except:
    print "*** MatPlotLib is unavailable! ***"
    MovieMode = False
# Our provided modules:
import mpi

try:
    assert os.getenv("SIDL_DLL_PATH")
except:
    print "Check your environment variables!"
#    rank = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
#    print "Printing Environement Variables on node %s" % (rank)
#    for key in os.environ.keys():
#        print "%s - %s: %s" % ( rank, key, os.environ[key] )

# Standard CPL6 Modules
import cpl
import cpl.attributevector
import cpl.comm
import cpl.control
import cpl.contract
import cpl.domain
import cpl.fields
import cpl.infobuffer
import cpl.map

# Coupler Only Modules
import frac
import flux
import areafact

from cpl.debug import *

#debugPrint = cpl.debug.newDebugPrint(True)
#statusPrint = cpl.debug.newDebugPrint(True)
#mapPrint = cpl.debug.newDebugPrint(True)

TIMESTEP = 0

# Data:  contents of data_mod.F90:
# Domains:
# dom_a, dom_i, dom_l, dom_r, dom_o
# self.contracts["c2a"].domain <=> dom_a
# self.contracts["c2a"].av
#NOT: self.contracts["c2a"].bundle.domain
# Contracts:
# these are now defined at initialization
# time in the coupler object.  (And are contained
# in the contracts dictionary.)
#con_Xc2a,con_Xc2l,con_XDc2l,con_Xc2o,con_Xc2i
# Coupler.contracts[""]

# Fundamental Maps:
map_Fa2i,map_Fa2l,map_Sa2i,map_Sa2l = 4*(None,)
map_Xr2o,map_ID = 2*(None,)

# Redundant maps:
"""
map_Fa2i, mapFa2l, map_Sa2i, map_Sa2l

map_Fi2a, map_Fi2l, map_Fi2o, map_Si2a,
map_Si2l, map_Si2o

map_Fl2a, map_Fl2i, map_Fl2o,
map_Sl2a, map_Sl2i, map_Sl2o

map_Fo2i, map_Fo2l, map_So2i, map_So2l
"""


class Coupler:
    def __init__(self,filepath="cpl.nml",computations=False,mapping=False,river=False):
        self.name = "cpl"
        self.doComputations = computations
        self.doMapping = mapping
        self.doRiver = river
        self.contracts={}
        startTime = time.time()
        print "pyCPL Init: starting @",time.ctime(),"!"
        local_comm = cpl.comm.init( cpl.fields.cplname )
        # print "pyCPL: Barrier 1 @",time.ctime(),": ..."
        # mpi.barrier( mpi.MPI_COMM_WORLD )
        print "OK!"
        rank = mpi.comm_rank( mpi.MPI_COMM_WORLD )
        size = mpi.comm_size( mpi.MPI_COMM_WORLD )
        print "mpi.MPI_COMM_WORLD size=%d,myrank=%d"%(size,rank)
        print "Size of Coupler Comm:",mpi.comm_size( local_comm )
        print "Rank in Coupler Comm:",mpi.comm_rank( local_comm )
        print "world_pid:",cpl.comm.world_pid
        print "mph_component_id:",cpl.comm.mph_component_id
        print "component_pid:",cpl.comm.component_pid
        print "component_npe:",cpl.comm.component_npe
        
        debugPrint( "Python Job running on %s, MPI Node %s"%(os.getenv("HOSTNAME"),rank) )
        if( cpl.comm.component_pid == 0 ):
            # fortran: call shr_msg_dirio('cpl')
            print "reading input namelist file:"
            cpl.control.init()
            self.namelist = cpl.control.readNamelist( filepath )
            print self.namelist
        debugPrint( "Broadcasting namelist to all coupler processors:" )
        ## Broadcast namelist to all processors of this
        # component.  (We'll serialize the namelist dict first)
        if ( cpl.comm.component_pid == 0 ):
            s_namelist = pickle.dumps( self.namelist )
            length_namelist = len(s_namelist)
        else:
            s_namelist = ""
            length_namelist = 0
        # broadcast length and data:
        length_namelist = mpi.bcast(length_namelist,1,
                                    mpi.MPI_INT,
                                    0,cpl.comm.local_comm )
        length_namelist = length_namelist[0]
        s_namelist = mpi.bcast( s_namelist, length_namelist,
                                mpi.MPI_CHAR, 0,
                                cpl.comm.local_comm )
        # Every processor in cpl.comm.local_comm
        # has s_namelist defined
        # Convert from array to string:
        s_namelist = "".join(s_namelist)
        self.namelist = pickle.loads( s_namelist )
        debugPrint( "Namelist broadcasted!" )
        # now we can set these values on all coupler processes.
        case_name = self.namelist["case_name"]
        self.date = self.namelist["start_date"]
        # Convert Date from Binary to Decimal Int:
        print "startdate:", self.date
        #self.date = int(self.date,2)
        info_dbug = self.namelist["info_dbug"]
        map_a2of_fn = self.namelist["map_a2of_fn"]
        map_a2os_fn = self.namelist["map_a2os_fn"]
        map_o2af_fn = self.namelist["map_o2af_fn"]
        map_r2o_fn = self.namelist["map_r2o_fn"]

        #cpl.control.ncpl_a = 24
        #cpl.control.ncpl_o = 1
        #cpl.control.ncpl_l = 12
        #cpl.control.ncpl_i = 24
        #cpl.control.ncpl_r = 12

        statusPrint("Setting up Infobuffers:")
        
        self.infobuf = cpl.infobuffer.Infobuffer()


        statusPrint("Contract Init:  establishes domains & routers (excluding land)")
        # Create contracts and
        # Initialize Info-buffers
        # Setup Fields:
        # fields_m2c_total = cpl.fields.a2c_total
        # fields_m2c_list = cpl.fields.a2c_fields
        # fields_c2m_total = cpl.fields.c2a_total
        # fields_c2m_list = cpl.fields.c2a_fields

        # ATM Contract
        # First Parallel Calls
        debugPrint("coupler.py thinks that the world communicator is:",cpl.comm.world_comm)
        debugPrint("coupler.py thinks that the local communicator is:",cpl.comm.local_comm)
        statusPrint("Initializing ATM Contracts...")
        self.contracts["a2c"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.atmname, cpl.control.decomp_a, cpl.fields.a2c_fields)
        self.contracts["c2a"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.atmname, cpl.control.decomp_a, cpl.fields.c2a_fields)
        # OCN contract
        statusPrint("Initializing OCN Contracts...")
        self.contracts["o2c"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.ocnname, cpl.control.decomp_a, cpl.fields.o2c_fields)
        self.contracts["c2o"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.ocnname, cpl.control.decomp_a, cpl.fields.c2o_fields)
        # ICE contract
        statusPrint("Initializing ICE Contracts...")
        self.contracts["i2c"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.icename, cpl.control.decomp_a, cpl.fields.i2c_fields)
        self.contracts["c2i"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.icename, cpl.control.decomp_a, cpl.fields.c2i_fields)

        # Special Domain Data Contract -- unique to lnd --
        statusPrint("Initializing Dc2l contract...")
        #print "fields:",cpl.fields.l2c_fields+cpl.fields.c2lg_fields
        self.contracts["Dc2l"] = cpl.contract.ReceivingContract( cpl.fields.cplname,
                                                                 cpl.fields.lndname,
                                                                 cpl.control.decomp_a,
                                                                 cpl.fields.c2lg_fields)
        self.dom_l = self.contracts["Dc2l"].domain
        if( self.contracts["Dc2l"].infobuffer.ibufGet("inimask") == 0 ):
            cpl.control.sendlnddom = False
        else:
            cpl.control.sendlnddom = True
        statusPrint("Initialized Dc2l contract!")

        ##
        statusPrint("Moved Partial Map Data Init to later in the script...")

        ##
        statusPrint("Send domain info to land model? (optional)")
        #statusPrint("SKIPPING LAND DOMAIN SEND!")
        if( not cpl.control.sendlnddom ):#not cpl.control.sendlnddom ):
            statusPrint(" * lnd DOES NOT request optional domain data exchange")
        else:
            statusPrint(" * lnd requests optional domain data exchange")
            print "--- initialize gData0..."
            gdata0 = cpl.attributevector.AttributeVector()
            ifields,rfields = self.contracts["Dc2l"].av.getFields()
            gdata0.initialize( ifields, rfields,
                               self.contracts["Dc2l"].av.av.lsize() )
            print "--- initialize gData1..."
            gdata1 = cpl.attributevector.AttributeVector()
            ifields,rfields = self.contracts["c2a"].domain.lgrid.getFields()
            gdata1.initialize( ifields, rfields,
                               self.contracts["c2a"].domain.lgrid.av.lsize() )
            print "--- initialize gData2..."
            gdata2 = cpl.attributevector.AttributeVector()
            ifields,rfields = frac.frac_a.getFields()
            gdata2.initialize( ifields, rfields,
                               frac.frac_a.size() )
            print "--- initialized gdata0,gdata1,gdata2! ---"
            
            print "Zeroing \'Dc2l' AttrVect..."
            #call cpl_bundle_zero(con_Dc2l%bundle)
            self.contracts["Dc2l"].av.av.zero()

            #print "Init Gdata0,1,2"
            if( False ):
                statusPrint("Disabling Gather(s) for now:")
            else:
                statusPrint("gData0 gather...")
                gdata0.av.gather(self.contracts["Dc2l"].av.av,
                              self.contracts["Dc2l"].domain.gsMap,
                              0, cpl.comm.component_pid,
                              cpl.comm.local_comm )
                mpi.barrier(cpl.comm.local_comm)
                statusPrint("gData1 gather...")
                gdata1.av.gather(self.contracts["c2a"].domain.lgrid.av,
                              self.contracts["c2a"].domain.gsMap,
                              0, cpl.comm.component_pid,
                              cpl.comm.local_comm )
                statusPrint("gData2 gather...")
                gdata2.av.gather(frac.frac_a.av,
                              self.contracts["c2a"].domain.gsMap,
                              0, cpl.comm.component_pid,
                              cpl.comm.local_comm )
                # --- set gData0 to be scattered as new data in con_Dc2l ---
                if (cpl.comm.component_pid == 0): 
                    # nfld = cpl_mct_aVect_indexRA(gData2,"lfrac",perrWith='gData2 lfrac')
                    # lsize = cpl_mct_aVect_lsize(gData0)
                    lfrac2,lfrac2_size = gdata2.av.exportRAttr("lfrac")
                    lsize = gdata0.av.lsize()
                    # do n=1,lsize
                    print gdata0
                    print gdata1
                    print gdata2
                    x,x_size = gdata1.av.exportRAttr("lon")
                    gdata0.av.importRAttr("lon",x)
                    x,x_size = gdata1.av.exportRAttr("lat")
                    gdata0.av.importRAttr("lat",x)
                    x,x_size = gdata1.av.exportRAttr("area")
                    gdata0.av.importRAttr("area",x)
                    mask,mask_size = gdata1.av.exportRAttr("mask")
                    gdata0.av.importRAttr("maska",x)
                    lmask,lmask_size = gdata0.av.exportRAttr("maskl")
                    lfrac0,lfrac0_size = gdata0.av.exportRAttr("lfrac")

                    for n in range(0,lsize-1):
                        if (lfrac2[n] < 1.0e-06):
                            lmask[n] = 0.0
                            lfrac0[n] = 0.0
                            mask[n] = 0.0
                        else:
                            lmask[n] = 1.0
                            lfrac0[n] = lfrac2[n]
                            mask[n] = 1.0
                        gdata0.av.importRAttr("maskl",lmask)
                        gdata0.av.importRAttr("lfrac",lfrac0)
                        gdata1.av.importRAttr("mask",mask)
                # endif
                        
            # --- reset dom_l based on dom_a ---
            # --- cpl_fields_grid_mask is not from dom_a (see above loop) ---
            # call cpl_mct_aVect_scatter(gData1,dom_l%lGrid, &
            #   dom_l%gsMap,0,cpl_comm_comp,rcode)

            self.dom_l.lgrid.av.scatter( gdata1.av,self.contracts["Dc2l"].domain.gsMap,0, cpl.comm.component_pid, cpl.comm.local_comm )
            
            # --- scatter gData0 to con_Dc2l bundle ---
            # call cpl_mct_aVect_scatter(gData0,con_Dc2l%bundle%data, &
            #   con_Dc2l%bundle%dom%gsMap,0,cpl_comm_comp,rcode)

            self.contracts["Dc2l"].av.av.scatter( gdata0.av,
                                                  self.contracts["Dc2l"].domain.gsMap, 0, cpl.comm.component_pid, cpl.comm.local_comm )

            # --- clean up ---
            # if (cpl_comm_comp_pid == 0) then
            #    call cpl_mct_aVect_clean(gData0)
            #    call cpl_mct_aVect_clean(gData1)
            #    call cpl_mct_aVect_clean(gData2)
            # endif

            # --- send con_Dc2l ---
            #call cpl_interface_contractSend(cpl_fields_lndname,con_Dc2l)
            print "Sending Dc2l Contract..."
            self.contracts["Dc2l"].send()
            print "Dc2l Contract Sent!"
        # END OPTIONAL LAND DOMAIN DATA EXCHANGE
        sys.stdout.flush()
        mpi.barrier( cpl.comm.local_comm )
        ## LND Contract
        statusPrint("Initializing LND Contracts...")
        statusPrint("Initializing LND: l2c")
        self.contracts["l2c"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.lndname, cpl.control.decomp_a, cpl.fields.l2c_fields)
        #self.contracts["l2c"] = cpl.contract.SendingContract( cpl.fields.cplname, cpl.fields.lndname, cpl.control.decomp_a, cpl.fields.l2c_fields)

        statusPrint("Initializing LND: c2l")
        self.contracts["c2l"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.lndname, cpl.control.decomp_a, cpl.fields.c2l_fields)

        ## Special runoff contract:
        statusPrint("Initializing r2c contract...")
        self.contracts["r2c"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.rtmname, cpl.control.decomp_a, cpl.fields.r2c_fields)

        ##
        statusPrint("Partial Map Data Init, Map Init")
        self.mapInit()
        frac.init( self.map_Fo2a,
                   self.contracts["a2c"].domain,
                   self.contracts["i2c"].domain,
                   self.contracts["Dc2l"].domain,
                   self.contracts["o2c"].domain)

        statusPrint("Init Data_Mod Bundles, Tavg History Bundles...")
        self.bundleInit()

        statusPrint("** Ignoring history bundles for now...")
        # history_avbundleInit()

        # --- START KLUDGE ---
        # fix_So2c_a is similar to bun_So2c_a but with a slightly different mapping
        # fix_frac_a is used to re-normalize map_Fo2a mapping for fix_So2c_a only
        # ---
        # call cpl_bundle_initv(fix_So2c_a,"fix_So2c_a",bun_So2c_a,bun_So2c_a%dom)
        # call cpl_bundle_initv(fix_frac_a,"fix_frac_a",bun_frac_a,bun_frac_a%dom)
        self.fix_so2c_a = cpl.attributevector.AttributeVector()
        self.fix_so2c_a.initv( self.av_so2c_a, self.av_so2c_a.size() )
        self.fix_frac_a = cpl.attributevector.AttributeVector()
        debugPrint("frac.frac_a.size():",frac.frac_a.size())
        self.fix_frac_a.initv( frac.frac_a, frac.frac_a.size() )
        # --- END KLUDGE ---
        
       
        # --- Check atm/lnd and ocn/ice model domains for consistency
        statusPrint("Checking ATM/LND & OCN/ICE domains for consistency...")
        statusPrint("Currently consistency checks are disabled!")
        # These appear to hang, I don't think they are
        # being executed on all processors.

        if( False ):
            statusPrint("Comparing a2c - l2c domains: ...")
            self.contracts["a2c"].domain.compare( self.contracts["l2c"].domain,
                                              enforce_grid=True,
                                              enforce_area=True,)
            statusPrint("Comparing o2c - i2c domains: ...")
            self.contracts["o2c"].domain.compare( self.contracts["i2c"].domain, enforce_all=True)


        # --- Send inital message with cday(cdate?), etc. ---
        debugPrint(" setting integer infobuffer fields...")
        self.infobuf.ibufSet('ncpl',24)
        self.infobuf.ibufSet('cdate',0)
        self.infobuf.ibufSet('rcode',0)
        self.infobuf.ibufSet('stopeod',0)
        self.infobuf.ibufSet('stopnow',0)
        self.infobuf.ibufSet('sec',0)
        self.infobuf.ibufSet('infotim',0)
        self.infobuf.ibufSet('infobug',cpl.control.infodbug)
        debugPrint(" setting real infobuffer fields...")
        self.infobuf.rbufSet('spval',0.0)
        self.infobuf.rbufSet('eccen',0.0)
        self.infobuf.rbufSet('obliqr',0.0)
        self.infobuf.rbufSet('lambm0',0.0)
        self.infobuf.rbufSet('mvelpp',0.0)

        statusPrint(" Initial send of infobuffer...")
        # Initial infobuffer send uses tag=1002
        self.infobuf.send( cpl.comm.world_pe0_atm, 1002 )
        self.infobuf.send( cpl.comm.world_pe0_ocn, 1002 )
        self.infobuf.send( cpl.comm.world_pe0_ice, 1002 )
        self.infobuf.send( cpl.comm.world_pe0_lnd, 1002 )
            
        statusPrint("Verifying acceptable coupling intervals:")
        ncpl_a = self.contracts["a2c"].infobuffer.ibufGet("ncpl")
        ncpl_i = self.contracts["i2c"].infobuffer.ibufGet("ncpl")
        ncpl_l = self.contracts["l2c"].infobuffer.ibufGet("ncpl")
        ncpl_r = self.contracts["r2c"].infobuffer.ibufGet("ncpl")
        ncpl_o = self.contracts["o2c"].infobuffer.ibufGet("ncpl")
        if( (ncpl_a < ncpl_o ) or ( (ncpl_a%ncpl_o)!=0 ) ):
            debugPrint("ERROR: Unacceptable ncpl_a, ncpl_o =",ncpl_a, ",",ncpl_o)
        elif((ncpl_a != ncpl_i) or (ncpl_a != ncpl_l)):
            debugPrint("Error Unacceptable ncpl_[ail] =",ncpl_a,",",ncpl_i,",",ncpl_l)
            statusPrint("Fudging ncpl_l value: setting it to ncpl_a(24):")
            ncpl_l = 24
        else:
            debugPrint("ncpl_[ailro] =",ncpl_a,",",ncpl_i,",",ncpl_l,",",ncpl_r,",",ncpl_o)
        statusPrint("Checking for dead models:")
        ndead = 0
        if(self.contracts["a2c"].infobuffer.ibufGet("dead")!=0):
            ndead += 1
            statusPrint("\tdead ATM component")
        if(self.contracts["i2c"].infobuffer.ibufGet("dead")!=0):
            ndead += 1
            statusPrint("\tdead ICE component")
        if(self.contracts["l2c"].infobuffer.ibufGet("dead")!=0):
            ndead += 1
            statusPrint("\tdead LND component")
        if(self.contracts["o2c"].infobuffer.ibufGet("dead")!=0):
            ndead += 1
            statusPrint("\tdead OCN component")
        if(ndead == 0):
            statusPrint("\tno dead components")

        statusPrint("Initializing Area Correction Data:")
        areafact.init(self.contracts['a2c'].domain,
                      self.contracts['i2c'].domain,
                      self.contracts['l2c'].domain,
                      self.contracts['o2c'].domain,
                      self.contracts['r2c'].domain)
        statusPrint("Initialized Area Correction Data!")
        # Initial receive from model(s)
        statusPrint(" Initial receive from Atmosphere Model...")
        self.contracts["a2c"].recv()
        statusPrint(" Initial receive from Ice Model...")
        self.contracts["i2c"].recv()
        statusPrint(" Initial receive from Land Model...")
        self.contracts["l2c"].recv()
        statusPrint(" Initial receive from River Model...")
        self.contracts["r2c"].recv()
        statusPrint(" Initial receive from Ocean Model...")
        self.contracts["o2c"].recv()
        statusPrint(" Initial receives complete!")
        # Read IC data from restart file?
        # Skipping restart stuff for now.

        # PROCESS IC DATA
        statusPrint("Processing Initial Condition Data:")
        statusPrint("Processing ATM IC data...")
        #    call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,bun_Fa2c_a)
        self.av_sa2c_a.copy( self.contracts["a2c"].av )
        self.av_fa2c_a.copy( self.contracts["a2c"].av )
        #    call cpl_bundle_zero (bun_precip_a)
        ifields,rfields = self.av_precip_a.getFields()
        temp_size = self.av_precip_a.size()
        debugPrint("temp_size:",temp_size)
        for f in rfields:
            self.av_precip_a.importRAttr(f,Numeric.zeros((temp_size),
                                                         Numeric.Float64))
        #    call cpl_bundle_add(bun_precip_a,'Faxc_rain',
        #                        bun_Fa2c_a,'Faxa_rainc')
        #    call cpl_bundle_add(bun_precip_a,'Faxc_rain',
        #                        bun_Fa2c_a,'Faxa_rainl')
        #    call cpl_bundle_add(bun_precip_a,'Faxc_snow',
        #                        bun_Fa2c_a,'Faxa_snowc')
        #    call cpl_bundle_add(bun_precip_a,'Faxc_snow',
        #                        bun_Fa2c_a,'Faxa_snowl')
        faxa_rainc = self.av_fa2c_a.exportRAttr( "Faxa_rainc")
        faxa_rainl = self.av_fa2c_a.exportRAttr( "Faxa_rainl")
        faxa_snowc = self.av_fa2c_a.exportRAttr( "Faxa_snowc")
        faxa_snowl = self.av_fa2c_a.exportRAttr( "Faxa_snowl")
        precip_a = faxa_rainc + faxa_rainl 
        self.av_precip_a.importRAttr( "Faxc_rain" , precip_a  )
        precip_a = faxa_snowc + faxa_snowl
        self.av_precip_a.importRAttr( "Faxc_snow" , precip_a )
        
        #cpl_control_fluxAShift = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ashift)
        cpl.control.fluxashift = self.contracts["a2c"].infobuffer.ibufGet("ashift")
        statusPrint("Processing ICE IC data...")
        #    lsize = cpl_mct_aVect_lsize(con_Xi2c%bundle%data)
        lsize = self.contracts["i2c"].av.size()
        #    allocate(ifrac_i(lsize))
        #    call cpl_mct_aVect_getRAttr(con_Xi2c%bundle%data,
        #                                "Si_ifrac",ifrac_i,rcode)
        ifrac_i = self.contracts["i2c"].av.exportRAttr("Si_ifrac")
        #    call set(ifrac_i,map_Fo2a,
        #                  dom_a,dom_i,dom_l,dom_o) ! init all fracs
        frac.set( ifrac_i, self.map_Fo2a,
                       self.contracts["a2c"].domain,
                       self.contracts["i2c"].domain,
                       self.contracts["l2c"].domain,
                       self.contracts["o2c"].domain )
        #    call cpl_bundle_split(con_Xi2c%bundle,bun_Si2c_i,bun_Fi2c_i)
        """
        copy cpl.fields.i2c_fluxes and cpl.fields.i2c_states into seperate bundles
        bundle_Si2c_i(States), bundle_Fi2c_i(Fluxes).

        create 2 av's
        av_Si2c_i = cpl.attributevector.AttributeVector("",cpl.fields.i2c_states,
                                   self.contracts["i2c"].av.lsize() )
        av_Fi2c_i = ...
        av_Si2c_i = self.contracts["i2c"].av.copy( av_Si2c_i )
        av_Fi2c_i =    ..                .av.copy( av_Fi2c_i )
        """
        #    call cpl_map_bun(bun_Si2c_i,bun_Si2c_a,map_So2a, &
        #                     bun_frac_i,'ifrac',bun_frac_a,'ifrac')
        
        #    call cpl_map_bun(bun_Fi2c_i,bun_Fi2c_a,map_Fo2a, &
        #                     bun_frac_i,'ifrac',bun_frac_a,'ifrac')
        statusPrint("Processing LND IC data...")
        #    call cpl_bundle_split(con_Xl2c%bundle,bun_Sl2c_l,bun_Fl2c_l )
        #    !--- map land fields to ocean, not allowed now.
        # !  call cpl_map_bun(bun_Sl2c_l,bun_Sl2c_o,map_Sa2o)
        # !  call cpl_map_bun(bun_Fl2c_l,bun_Fl2c_o,map_Fa2o)
        #    !--- map land fields to ocean, not allowed now.
        statusPrint("Processing OCN IC data...")
        #    !--- process ocn IC data ---
        #    call cpl_bundle_split(con_Xo2c%bundle,bun_So2c_o,bun_Fo2c_o)
        #    call cpl_map_bun(bun_So2c_o,bun_So2c_a,map_So2a, &
        #                     bun_frac_o,'afrac',bun_frac_a,'ofrac')
        #    !*** KLUDGE - start *********
        #    call cpl_map_bun(bun_frac_o,fix_frac_a,map_So2a)
        #    call cpl_map_bun(bun_So2c_o,fix_So2c_a,map_So2a, &
        #                     bun_frac_o,'afrac',fix_frac_a,'afrac')
        #    !*** KLUDGE - end ***********
        #    call cpl_map_bun(bun_Fo2c_o,bun_Fo2c_a,map_Fo2a, &
        #                     bun_frac_o,'afrac',bun_frac_a,'ofrac')
        #    call cpl_map_bun(bun_aoflux_o,bun_aoflux_a,map_Fo2a, &
        #                     bun_frac_o,'afrac',bun_frac_a,'ofrac')
        #    call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
        #                     bun_frac_o,'afrac',bun_frac_a,'ofrac')

        statusPrint("Optional ATM initialization: send albedos")
        if( self.contracts["a2c"].infobuffer.ibufGet("xalbic")==0):
            cpl.control.sendatmalb = False
            statusPrint("* atm component requests NO recalculation of initial solar")
        else:
            cpl.control.sendatmalb = True
            statusPrint("* atm component requests recalculation of initial solar")
            statusPrint("\tMerging ATM inputs:")
            # Only albedo, surface temp, snow, ifrac and ofrac
            # need to be valid.
            # MERGE_ATM? where is this defined?
            # call merge_atm( fix_So2c_a ) # KLUDGE!

            self.merge_atm( self.fix_so2c_a )

            statusPrint("\tSend albedos to ATM, recv new ATM IC's")
            #     call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'cpl2comp',  &
            #               bunlist=cpl_fields_c2a_fluxes)
            #     call cpl_interface_contractSend(cpl_fields_atmname,con_Xc2a)
            self.contracts["c2a"].send()
            #     call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'comp2cpl',  &
            #               bunlist=cpl_fields_c2a_fluxes)
            statusPrint("\twait for new atm IC data...")
            # call cpl_interface_contractRecv(cpl_fields_atmname,con_Xa2c)
            # call cpl_bundle_mult(con_Xa2c%bundle,bun_areafact_a,'comp2cpl',  &
            #               bunlist=cpl_fields_a2c_fluxes)

            statusPrint("\tprocess atm IC data...")
            # call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,bun_Fa2c_a)
            # call cpl_bundle_zero (bun_precip_a)
            # call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainc')
            # call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainl')
            # call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowc')
            # call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowl')
            # cpl_control_fluxAShift = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ashift)
        # End Optional ATM Initialization"

        statusPrint("Create data as necessary for 1st iteration of main event loop")
        # --- map ocn & ice albedos onto atm domain ---
        # call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
        #            bun_frac_o,'afrac',bun_frac_a,'ofrac')
        # call cpl_map_bun(con_Xr2c%bundle,bun_Xr2c_o,map_Xr2o)

        # --- prepare cpl->atm bundle (as necessary for 1st diag calc) ---
        # --- merge atm inputs (necessary for 1st diag calc) ---
        # !  call merge_atm()
        # call merge_atm(fix_So2c_a)  ! KLUDGE

        #statusPrint("Merging ATM (necessary for 1st diagnostic calculation)")
        self.merge_atm(self.fix_so2c_a)
        
        # Clean up
        finishTime = time.time()
        statusPrint("pyCPL Init: Elapsed Time:",(finishTime - startTime))
        sys.stdout.flush()
        return
    
    def mapInit(self):
        statusPrint("Initializing mapping data...")
        ## initialize idenity map
        self.map_id = cpl.map.Map()
        self.map_id.idtype = 1
        ## Initialize maps between atm & ocn:
        map_a2of_fn = self.namelist["map_a2of_fn"]
        map_a2os_fn = self.namelist["map_a2os_fn"]
        map_o2af_fn = self.namelist["map_o2af_fn"]
        map_r2o_fn = self.namelist["map_r2o_fn"]
        
        #mpi.barrier( mpi.MPI_COMM_WORLD )
        debugPrint(map_o2af_fn)
        debugPrint(map_a2os_fn)
        debugPrint(map_a2of_fn)
        debugPrint(map_r2o_fn)

        statusPrint("\tInitializing map_Sa2o:...")
        self.map_Sa2o = cpl.map.Map()
        self.map_Sa2o.initialize( self.contracts["a2c"].domain,
                                  self.contracts["o2c"].domain,
                                  "Sa2o",
                                  os.path.join("data",map_a2of_fn),
                                  "src")
        statusPrint("\tInitializing map_Sa2o:OK")
        
        statusPrint("\tInitializing map_Fa2o:...")
        self.map_Fa2o = cpl.map.Map()
        self.map_Fa2o.initialize( self.contracts["a2c"].domain,
                                  self.contracts["o2c"].domain,
                                  "Fa2o",
                                  os.path.join("data",map_a2of_fn),
                                  "src")
        statusPrint("\tInitializing map_Fa2o:OK")

        statusPrint("\tInitializing map_So2a:...")
        self.map_So2a = cpl.map.Map()
        self.map_So2a.initialize( self.contracts["o2c"].domain,
                                  self.contracts["a2c"].domain,
                                  "So2a",
                                  os.path.join("data",map_o2af_fn),
                                  "dst")
        statusPrint("\tInitializing map_So2a:OK")

        statusPrint("\tInitializing map_Fo2a:...")
        self.map_Fo2a = cpl.map.Map()
        self.map_Fo2a.initialize( self.contracts["o2c"].domain,
                                  self.contracts["a2c"].domain,
                                  "Fo2a",
                                  os.path.join("data",map_o2af_fn),
                                  "dst" )
        statusPrint("\tInitializing map_Fo2a:OK")

        if (self.doRiver):
            statusPrint("\tInitializing map_Xr2o:...")
            self.map_Xr2o = cpl.map.Map()
            self.map_Xr2o.initialize( self.contracts["l2c"].domain,
                                      self.contracts["o2c"].domain,
                                      "Xr2o",
                                      os.path.join("data",map_r2o_fn),
                                      "dst")
            statusPrint("\tInitializing map_Xr2o:OK")
        else:
            statusPrint("\t** Skipping MAP_XR2O Init to speed debugging **")
        return
    
    def bundleInit(self):
        """
        -------------------------------------------------------------------------------
        METHOD:
        1) assign the domain (an input) to the bundle
        2) set size of bundle aVect to be the same size as that of the domain
        3) use input str to define real aVect fields (must be a valid aVect list str)
        4) hard-coded st the aVect has *no* integer data
        NOTE:
            o memory is allocated for bun%data, but data values are left undefined
        -------------------------------------------------------------------------------
        """
        # Bundles:
        # Initialize these:
        size_o = self.contracts["o2c"].domain.lsize()
        size_a = self.contracts["a2c"].domain.lsize()
        size_l = self.contracts["l2c"].domain.lsize()
        size_r = self.contracts["r2c"].domain.lsize()
        size_i = self.contracts["i2c"].domain.lsize()

        self.av_xc2osnap_o = cpl.attributevector.AttributeVector()
        self.av_xc2opsum_o = cpl.attributevector.AttributeVector()

        self.av_aoflux_o = cpl.attributevector.AttributeVector()
        self.av_aoflux_a = cpl.attributevector.AttributeVector()

        self.av_oalbedo_o = cpl.attributevector.AttributeVector()
        self.av_oalbedo_a = cpl.attributevector.AttributeVector()

        self.av_precip_o = cpl.attributevector.AttributeVector()
        self.av_precip_a = cpl.attributevector.AttributeVector()

        self.av_sa2c_a = cpl.attributevector.AttributeVector()
        self.av_fa2c_a = cpl.attributevector.AttributeVector()
        self.av_sa2c_o = cpl.attributevector.AttributeVector()
        self.av_fa2c_o = cpl.attributevector.AttributeVector()
        self.av_sl2c_l = cpl.attributevector.AttributeVector()
        self.av_fl2c_l = cpl.attributevector.AttributeVector()
        self.av_xr2c_o = cpl.attributevector.AttributeVector()
        self.av_so2c_o = cpl.attributevector.AttributeVector()
        self.av_fo2c_o = cpl.attributevector.AttributeVector()
        self.av_so2c_a = cpl.attributevector.AttributeVector()
        self.av_fo2c_a = cpl.attributevector.AttributeVector()
        self.av_si2c_i = cpl.attributevector.AttributeVector()
        self.av_fi2c_i = cpl.attributevector.AttributeVector()
        self.av_si2c_a = cpl.attributevector.AttributeVector()
        self.av_fi2c_a = cpl.attributevector.AttributeVector()
        # Initialize:
        self.av_xc2osnap_o.init( [], cpl.fields.c2o_fields, size_o )
        self.av_xc2opsum_o.init( [], cpl.fields.c2o_fields, size_o )

        aoflux_fields = ["Faoc_sen",
                         "Faoc_lat",
                         "Faoc_lwup",
                         "Faoc_evap",
                         "Faoc_taux",
                         "Faoc_tauy",
                         "Faoc_tref",
                         "Faoc_qref",
                         "Faoc_duu10n",
                         "Faoc_swnet",
                         # Added these for fluxEpbal call
                         "aream",
                         "mask"]
        self.av_aoflux_o.init([], aoflux_fields, size_o )
        self.av_aoflux_a.init([], aoflux_fields, size_a )
        albedo_fields = ["So_avsdr","So_anidr","So_avsdf","So_anidf"]
        self.av_oalbedo_o.init([], albedo_fields, size_o )
        self.av_oalbedo_a.init([], albedo_fields, size_a )
        precip_fields = ["Faxc_rain","Faxc_snow"]
        self.av_precip_o.init([],precip_fields,size_o)
        self.av_precip_a.init([],precip_fields,size_a)
        #
        self.av_sa2c_a.init([],cpl.fields.a2c_states,size_a)
        self.av_fa2c_a.init([],cpl.fields.a2c_fluxes,size_a)
        self.av_sa2c_o.init([],cpl.fields.a2c_states,size_o)
        self.av_fa2c_o.init([],cpl.fields.a2c_fluxes,size_o)
        #
        self.av_sl2c_l.init([],cpl.fields.l2c_states,size_l)
        self.av_fl2c_l.init([],cpl.fields.l2c_fluxes,size_l)
        #self.av_sl2c_o.init([],cpl.fields.l2c_states,size_o)
        #self.av_fl2c_o.init([],cpl.fields.l2c_fluxes,size_o)
        #
        self.av_xr2c_o.init([],cpl.fields.r2c_fields,size_o)
        #
        self.av_so2c_o.init([],cpl.fields.o2c_states,size_o)
        self.av_fo2c_o.init([],cpl.fields.o2c_fluxes,size_o)
        self.av_so2c_a.init([],cpl.fields.o2c_states,size_a)
        self.av_fo2c_a.init([],cpl.fields.o2c_fluxes,size_a)
        #
        self.av_si2c_i.init([],cpl.fields.i2c_states,size_i)
        self.av_fi2c_i.init([],cpl.fields.i2c_fluxes,size_i)
        self.av_si2c_a.init([],cpl.fields.i2c_states,size_a)
        self.av_fi2c_a.init([],cpl.fields.i2c_fluxes,size_a)
        return

        
    def run( self, days, LogFile=None ):
        """
        Timestep the model for 'days' days.

        (Maybe we should seperate out all the 'one time'
        code in here and put it into the init method.

        Then we can consider making just a simple
        'step' routine that contains all the code that
        needs to occur daily.  That way this method ('run')
        could be simplified down to simply a

        for i in days:
            self.step()
        )
        """
        # Main Integration Loop: repeat until cpl says stop:
        print self.name,": beginning main integration loop..."
        print self.name,"Running for",days,"days."
        debugPrint("NCPL:\tVALUE")
        debugPrint("atm:\t",cpl.control.ncpl_a)
        debugPrint("ocn:\t",cpl.control.ncpl_o)
        debugPrint("ice:\t",cpl.control.ncpl_i)
        debugPrint("lnd:\t",cpl.control.ncpl_l)
        debugPrint("rtm:\t",cpl.control.ncpl_r)
        iterations=0 # outer loop index
        done = False
        stopnow = 0
        ncpl = self.infobuf.ibufGet('ncpl')
        while not cpl.control.stopeod:
            #print "cpl: Top of outer loop:"
            iterations+=1
            if (iterations >= days):
                print "Setting Stop Now!"
                cpl.control.stopeod=1
                stopeod = 1
                for c in ["c2a","c2i","c2o","c2l","r2c"]:
                    self.contracts[c].infobuffer.ibufSet('stopeod',stopeod)
            # Beginning of NCPL Loop
            # 1 to cpl.control.ncpl+1(excludes end, so 1-24)
            for n in xrange(1,cpl.control.ncpl_a+1):
                print "cpl: Starting hour %s, %s..."%(n,cpl.control.ncpl_a)
                a = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_a)
                o = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_o)
                i = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_i)
                l = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_l)
                r = (n)%(cpl.control.ncpl_a/cpl.control.ncpl_r)
                print "a,o,i,l,r = ",a,o,i,l,r
                
                # Send message to ocn
                """
                if(mod(n-1,ncpl_a/ncpl_o) == 0 ) then
                if (.not. cpl_control_lagOcn) then
                call cpl_bundle_mult(con_Xc2o%bundle,bun_areafact_o,'cpl2comp',  &
                bunlist=cpl_fields_c2o_fluxes)
                call cpl_interface_contractSend(cpl_fields_ocnname,con_Xc2o)
                call cpl_bundle_mult(con_Xc2o%bundle,bun_areafact_o,'comp2cpl',  &
                bunlist=cpl_fields_c2o_fluxes)
                endif
                call cpl_bundle_zero(bun_Xc2oPSUM_o) ! zero-out partial sum
                endif
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_o) == 0 ):
                    if (not cpl.control.lagocn):
                        ifields,rfields=self.contracts["c2o"].av.getFields()
                        cpl2comp=areafact.av_areafact_o.exportRAttr("cpl2comp")
                        comp2cpl=areafact.av_areafact_o.exportRAttr("comp2cpl")
                        # Mult by cpl2comp
                        for r in rfields:
                            temp = self.contracts["c2o"].av.exportRAttr(r)
                            temp *= cpl2comp
                            self.contracts["c2o"].av.importRAttr(r,temp)
                        # Send OCN
                        self.contracts["c2o"].send()
                        # Mult by comp2cpl
                        for r in rfields:
                            temp = self.contracts["c2o"].av.exportRAttr(r)
                            temp *= comp2cpl
                            self.contracts["c2o"].av.importRAttr(r,temp)
                    ifields,rfields = self.av_xc2opsum_o.getFields()
                    size = self.av_xc2opsum_o.size()
                    # Zero-out partial sum:
                    for r in rfields:
                        self.av_xc2opsum_o.importRAttr(r,
                                                       Numeric.zeros((size),
                                                            Numeric.Float64))
                        
                # Send message to land
                """
                call cpl_bundle_gather(con_Xc2l%bundle, bun_Sa2c_a, bun_Fa2c_a,
                &                      bun_Sl2c_l, bun_Fl2c_l,
                &                      bun_So2c_a, bun_Fo2c_a,
                &                      bun_Si2c_a, bun_Fi2c_a  ) 
                call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,'cpl2comp',
                                 bunlist=cpl_fields_c2l_fluxes)
                call cpl_interface_contractSend(cpl_fields_lndname,con_Xc2l)
                call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,'comp2cpl',
                                     bunlist=cpl_fields_c2l_fluxes)
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_l) == 0 ):
                    # Gather:
                    self.contracts["c2l"].av.copy( self.av_sa2c_a )
                    self.contracts["c2l"].av.copy( self.av_fa2c_a )
                    self.contracts["c2l"].av.copy( self.av_sl2c_l )
                    self.contracts["c2l"].av.copy( self.av_fl2c_l )
                    self.contracts["c2l"].av.copy( self.av_so2c_a )
                    self.contracts["c2l"].av.copy( self.av_fo2c_a )
                    self.contracts["c2l"].av.copy( self.av_si2c_a )
                    self.contracts["c2l"].av.copy( self.av_fi2c_a )
                    # Mult, send, mult
                    ifields,rfields=self.contracts["c2l"].av.getFields()
                    cpl2comp=areafact.av_areafact_l.exportRAttr("cpl2comp")
                    comp2cpl=areafact.av_areafact_l.exportRAttr("comp2cpl")
                    # Mult by cpl2comp
                    for r in cpl.fields.c2l_fluxes:
                        temp = self.contracts["c2l"].av.exportRAttr(r)
                        temp *= cpl2comp
                        self.contracts["c2l"].av.importRAttr(r,temp)
                    # Send
                    print "** calling self.contracts[\"c2l\"].send()... **"

                    #self.contracts["Dc2l"].send()
                    self.contracts["c2l"].send()


                    print "** returned from self.contracts[\"c2l\"].send()! **"
                    # Mult by comp2cpl
                    for r in cpl.fields.c2l_fluxes:
                        temp = self.contracts["c2l"].av.exportRAttr(r)
                        temp *= comp2cpl
                        self.contracts["c2l"].av.importRAttr(r,temp)
                        
                print "** Compute, map, merge ice inputs"
                if (True):
                    # This is the first time bun_precip_[a|o] is really used:
                    # call cpl_map_bun(bun_precip_a,bun_precip_o,
                    #                  map_Fa2o,mvector=a2ovector)
                    self.av_precip_o = self.map_Fa2o.mapAV( self.av_precip_a )

                    # call cpl_map_bun(bun_Sa2c_a,
                    #                  bun_Sa2c_o,
                    #                  map_Sa2o,mvector=a2ovector)
                    self.av_sa2c_o = self.map_Sa2o.mapAV( self.av_sa2c_a )
                    # call cpl_map_bun(bun_Fa2c_a,
                    #                  bun_Fa2c_o,map_Fa2o,mvector=a2ovector)
                    self.av_fa2c_o = self.map_Fa2o.mapAV( self.av_fa2c_a )

                    
                # Force zero net water flux into ocn+ice ?
                if (self.doComputations):
                    if (cpl.control.fluxepbal != 'off'): 
                        # cpl_control_fluxepfac = con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_precAdj)
                        cpl.control.fluxepfac = self.contracts["o2c"].infobuffer.ibufGet("precadj")
                        # cpl_control_fluxEPfac = cpl_control_fluxEPfac * 1.0e-6_r8
                        cpl.control.fluxepfac = cpl.control.fluxepfac * 1.0e-6
                        # call flux_epbal(date,bun_aoflux_o,con_Xi2c%bundle, 
                        #                bun_precip_o,bun_Xr2c_o,bun_frac_o)
                        print "*** flux.fluxEpbal: av_aoflux_o is unitialized!"
                        print "*** flux.fluxEpbal is being skipped!"
                        #results = flux.fluxEpbal(self.date,
                        #               self.av_aoflux_o,
                        #               self.contracts["i2c"].av,
                        #               self.av_precip_o,
                        #               self.contracts["r2c"].av,
                        #               frac.frac_o)
                        #av_precip_o, self.contracts["r2c"].av = results
                    else:
                        print "*** Skipping call to flux.fluxEpbal (...) !!"
                # Merge total snow and precip for ice input
                if (self.doComputations):
                    # call cpl_bundle_zero (con_Xc2i%bundle,'Faxc_rain')
                    # call cpl_bundle_zero (con_Xc2i%bundle,'Faxc_snow')
                    zero=Numeric.zeros((self.contracts["c2i"].av.size()),Numeric.Float64)
                    self.contracts["c2i"].av.importRAttr("Faxc_rain",zero)
                    self.contracts["c2i"].av.importRAttr("Faxc_snow",zero)

                    # call cpl_bundle_copy(bun_precip_o,bunrList='Faxc_rain',&
                    #                     bunTrList='Faxc_rain',outbun=con_Xc2i%bundle)
                    # call cpl_bundle_copy(bun_precip_o,bunrList='Faxc_snow',&
                    #                     bunTrList='Faxc_snow',outbun=con_Xc2i%bundle)
                    print "av_precip_o size:",self.av_precip_o.size()
                    print "ice av size:",self.contracts["c2i"].av.size()
                    self.contracts["c2i"].av.copy(self.av_precip_o)
                    
                # Correct a->o vector mapping near North Pole(NP)
                if (self.doMapping):
                    print "*** Skipping NPFix for now!"
                    # call cpl_map_npfix(bun_Sa2c_a,bun_Sa2c_o,'Sa_u','Sa_v')
                
                # Send message to ice
                """
                call cpl_bundle_gather(con_Xc2i%bundle, bun_Sa2c_o, bun_Fa2c_o, &
!--- not allowed now ------------------------------ bun_Sl2c_o, bun_Fl2c_o, &
            &                                       bun_So2c_o, bun_Fo2c_o, &
            &                                       bun_Si2c_i, bun_Fi2c_i, &
            &                                       fcopy = .true.  )
            call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,'cpl2comp', &
                                 bunlist=cpl_fields_c2i_fluxes)
            call cpl_interface_contractSend(cpl_fields_icename,con_Xc2i)
            call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,'comp2cpl', &
                                 bunlist=cpl_fields_c2i_fluxes)
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_i) == 0 ):
                    # Gather:
                    self.contracts["c2i"].av.copy( self.av_sa2c_o )
                    self.contracts["c2i"].av.copy( self.av_fa2c_o )
                    self.contracts["c2i"].av.copy( self.av_so2c_o )
                    self.contracts["c2i"].av.copy( self.av_fo2c_o )
                    self.contracts["c2i"].av.copy( self.av_si2c_i )
                    self.contracts["c2i"].av.copy( self.av_fi2c_i )
                    # Mult, send, mult
                    ifields,rfields=self.contracts["c2i"].av.getFields()
                    cpl2comp=areafact.av_areafact_i.exportRAttr("cpl2comp")
                    comp2cpl=areafact.av_areafact_i.exportRAttr("comp2cpl")
                    # Mult by cpl2comp
                    for r in cpl.fields.c2i_fluxes:
                        temp = self.contracts["c2i"].av.exportRAttr(r)
                        temp *= cpl2comp
                        self.contracts["c2i"].av.importRAttr(r,temp)
                    # Send:
                    self.contracts["c2i"].send()
                    # Mult by comp2cpl
                    for r in cpl.fields.c2i_fluxes:
                        temp = self.contracts["c2i"].av.exportRAttr(r)
                        temp *= comp2cpl
                        self.contracts["c2i"].av.importRAttr(r,temp)
                
                # Compute net solar flux into ocn
                # CALL flux_solar( bun_Fa2c_o, bun_oalbedo_o, bun_aoflux_o )
                # Merge ocn inputs
                self.merge_ocn()
                
                # Form partial sum of tavg ocn inputs (virtual "send" to ocn)
                # CALL cpl_bundle_accum( bun_Xc2oSNAP_o, outbun=bun_Xc2oPSUM_o)
                
                # Write bit-for-bit check info
                # do diagnositcs
                # Skipping Diagnostics and B4B check!
                
                # Compute ocn albedos (virtual "recv" from ocn)
                # CALL flux_albo(date, bun_oalbedo_o)
                
                # Compute atm/ocn fluxes
                # CALL flux_atmOcn(con_Xo2c%bundle, bun_Sa2c_o, cpl_control_dead_ao, bun_aoflux_o )
                
                # Recv msg from ice
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_i) == 0 ):
                    self.contracts["i2c"].recv()
                    # Mult:
                    ifields,rfields=self.contracts["c2i"].av.getFields()
                    comp2cpl=areafact.av_areafact_i.exportRAttr("comp2cpl")
                    for r in cpl.fields.c2i_fluxes:
                        temp = self.contracts["c2i"].av.exportRAttr(r)
                        temp *= comp2cpl
                        self.contracts["c2i"].av.importRAttr(r,temp)
                    # update surface fracs wrt new ice frac
                    ifrac_i = self.contracts["i2c"].av.exportRAttr("Si_ifrac")
                    frac.set(ifrac_i, self.map_Fo2a,
                                  self.contracts["a2c"].domain,
                                  self.contracts["i2c"].domain,
                                  self.contracts["l2c"].domain,
                                  self.contracts["o2c"].domain)

                # Export the ice temperature data to a logfile:
                global TIMESTEP
                if(LogFile):
                    if (cpl.comm.world_pid == cpl.comm.world_pe0):
                        ice_temperature = self.contracts["i2c"].av.exportRAttr("Si_t")
                        # print ice_temperature
                        x = self.contracts["i2c"].domain.ni
                        y = self.contracts["i2c"].domain.nj
                        if(len(ice_temperature) == (x*y)):
                            #LogFile.var("ice_temp").put(ice_temperature)#[,start=(time,0)])
                            icet = Numeric.resize(ice_temperature,(x,y))
                            LogFile.var("ice_temp").put(icet,start=(TIMESTEP,0,0))
                            LogFile.sync()
                            #LogFile.close()
                            print "NetCDF file written!"
                            
                if(MovieMode):
                    if (cpl.comm.world_pid == cpl.comm.world_pe0):
                        ice_temperature = self.contracts["i2c"].av.exportRAttr("Si_t")
                        # print ice_temperature
                        x = self.contracts["i2c"].domain.ni
                        y = self.contracts["i2c"].domain.nj
                        if(len(ice_temperature) == (x*y)):
                            #LogFile.var("ice_temp").put(ice_temperature)#[,start=(time,0)])
                            icet = Numeric.resize(ice_temperature,(x,y))
                            image = pylab.figimage(icet)
                            path = os.path.join("movies","im"+'%03d'%TIMESTEP+".png")
                            image.write_png(path)
                            pylab.clf()
                            print path
                TIMESTEP += 1
                
                # Add diurnal cycle to ice albedos
                # CALL flux_albi( date, con_Xi2c%bundle )
                # CALL cpl_bundle_split( con_Xi2c%bundle, bun_Si2c_i,
                #                        bun_Fi2c_i )
                
                # Map ocn states/fluxes to atm -- BEGIN
                if( self.doMapping or True ):
                    print "\t<<< Mapping: av_si2c_i -> av_si2c_a >>>"
                    mapPrint( "\t<<< Checking Inputs: >>>" )
                    mapPrint( "\t<<< frac.frac_i size = %s >>>"%(frac.frac_i.size()) )
                    mapPrint( "\t<<< domain_i info: >>>" )
                    #self.contracts["c2i"].domain.info()
                    mapPrint( "\t<<< frac.frac_a size = %s >>>"%(frac.frac_a.size()) )
                    mapPrint( "\t<<< domain_a info: >>>" )
                    #self.contracts["c2a"].domain.info()
                    mapPrint( "\t<<< input AV: >>>" )
                    mapPrint( self.av_si2c_i )
                    mapPrint( "\t<<< calling mapAV >>>" )
                    av_si2c_a = self.map_So2a.mapAV( self.av_si2c_i,
                                                     frac.frac_i,"ifrac",
                                                     frac.frac_a,"ifrac")
                    print "\t<<< Mapping: av_so2c_o -> av_so2c_a >>>"
                    av_so2c_a = self.map_So2a.mapAV( self.av_so2c_o,
                                                     frac.frac_o,"afrac",
                                                     frac.frac_a,"ofrac")
                    # KLUDGE - START
                    print "\t<<< Mapping: frac.frac_o -> fix_frac_a >>>"
                    fix_frac_a = self.map_So2a.mapAV( frac.frac_o )
                    print "\t<<< Mapping: av_so2c_o -> fix_so2c_a >>>"
                    fix_so2c_a = self.map_So2a.mapAV( self.av_so2c_o,
                                                      frac.frac_o,"afrac",
                                                      fix_frac_a,"afrac" )
                    # KLUDGE - END
                    """
                    call cpl_map_bun(bun_Fi2c_i,bun_Fi2c_a,map_Fo2a, &
                    bun_frac_i,'ifrac',bun_frac_a,'ifrac',oi2avector)
                    call cpl_map_bun(bun_Fo2c_o,bun_Fo2c_a,map_Fo2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac',oi2avector)
                    call cpl_map_bun(bun_aoflux_o,bun_aoflux_a,map_Fo2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac',oi2avector)
                    call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac',oi2avector)
                    """
                    print "\t<<< Mapping: av_fi2c_i -> av_fi2c_a >>>"
                    self.av_fi2c_a = self.map_Fo2a.mapAV(self.av_fi2c_i,
                                                         frac.frac_i,"ifrac",
                                                         frac.frac_a,"ifrac")
                    print "\t<<< Mapping: av_fo2c_o -> av_fo2c_a >>>"
                    self.av_fo2c_a = self.map_Fo2a.mapAV(self.av_fo2c_o,
                                                         frac.frac_o,"afrac",
                                                         frac.frac_a,"ofrac")
                    print "\t<<< Mapping: av_aoflux_o -> av_aoflux_a >>>"
                    self.av_aoflux_a = self.map_Fo2a.mapAV(self.av_aoflux_o,
                                                           frac.frac_o,"afrac",
                                                           frac.frac_a,"ofrac")
                    print "\t<<< Mapping: av_oalbedo_a -> av_oalbedo_o >>>"
                    self.av_oalbedo_a = self.map_Fo2a.mapAV(self.av_oalbedo_o,
                                                            frac.frac_o,"afrac",
                                                            frac.frac_a,"ofrac")
                # Map ocn states/fluxes to atm -- END
                    
                # Recv message from lnd
                """
                call cpl_interface_contractRecv(cpl_fields_lndname,con_Xl2c)
                call cpl_bundle_mult(con_Xl2c%bundle,bun_areafact_l,'comp2cpl',
                                     bunlist=cpl_fields_l2c_fluxes)
                call cpl_bundle_split(con_Xl2c%bundle,bun_Sl2c_l,bun_Fl2c_l )
   
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_l) == 0 ):
                    self.contracts["l2c"].recv()
                    # mult by comp2cpl
                    ifields,rfields=self.contracts["l2c"].av.getFields()
                    comp2cpl=areafact.av_areafact_l.exportRAttr("comp2cpl")
                    for r in cpl.fields.c2l_fluxes:
                        temp = self.contracts["c2l"].av.exportRAttr(r)
                        temp *= comp2cpl
                        self.contracts["c2l"].av.importRAttr(r,temp)
                    # split:
                    self.av_sl2c_l.copy( self.contracts["c2l"].av )
                    self.av_fl2c_l.copy( self.contracts["c2l"].av )
                    
                # recv from runoff next:
                # possible TYPO:  n is used here instead of n-1
                """
                call cpl_interface_contractRecv(cpl_fields_lndname,con_Xr2c)
                call cpl_bundle_mult(con_Xr2c%bundle,bun_areafact_r,'comp2cpl',
                                 bunlist=cpl_fields_r2c_fluxes)
                """
                if ( (n)%(cpl.control.ncpl_a/cpl.control.ncpl_r) == 0 ):
                    self.contracts["r2c"].recv()
                    # mult by comp2cpl
                    ifields,rfields=self.contracts["r2c"].av.getFields()
                    comp2cpl=areafact.av_areafact_r.exportRAttr("comp2cpl")
                    for r in cpl.fields.r2c_fluxes:
                        temp = self.contracts["r2c"].av.exportRAttr(r)
                        temp *= comp2cpl
                        self.contracts["r2c"].av.importRAttr(r,temp)
                
                # diagnostics: verify net solar calcs are coordinated
                # Merge atm states and fluxes

                print "Merging ATM states and fluxes:"
                # CALL MERGE_ATM(fix_so2c_a) #+KLUDGE (fix_so2c_a is a kludge)
                #print "Calling Merge OCN inplace of Merge ATM!"
                #self.merge_ocn()
                self.merge_atm( self.fix_so2c_a )
                
                # Send message to atm
                """
                call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'cpl2comp',
                                 bunlist=cpl_fields_c2a_fluxes)
                call cpl_interface_contractSend(cpl_fields_atmname,con_Xc2a)
                call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'comp2cpl',
                                 bunlist=cpl_fields_c2a_fluxes)
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_a) == 0 ):
                    # Mult, send, mult
                    ifields,rfields=self.contracts["c2a"].av.getFields()
                    cpl2comp=areafact.av_areafact_a.exportRAttr("cpl2comp")
                    comp2cpl=areafact.av_areafact_a.exportRAttr("comp2cpl")
                    # Mult by cpl2comp
                    for r in cpl.fields.c2a_fluxes:
                        temp = self.contracts["c2a"].av.exportRAttr(r)
                        temp *= cpl2comp
                        self.contracts["c2a"].av.importRAttr(r,temp)
                    # Send
                    self.contracts["c2a"].send()
                    # Mult by comp2cpl
                    for r in cpl.fields.c2a_fluxes:
                        temp = self.contracts["c2a"].av.exportRAttr(r)
                        temp *= comp2cpl
                        self.contracts["c2a"].av.importRAttr(r,temp)
                # Map land to ocean, NOT allowed now
                
                
                if (self.doMapping):
                    """call cpl_map_bun(con_Xr2c%bundle,bun_Xr2c_o,
                    #                 map_Xr2o,mvector=r2ovector)
                    """
                    pass
                
                # Create history files
                # recv msg from OCN
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_o) == 0 ):
                    # --- form tavg of ocn inputs:
                    #call cpl_bundle_avg(bun_Xc2oPSUM_o)
                    # Average a bundle:  basically this divides all fields in
                    # a bndle by the value of the bundle counter (bun%cnt)
                    # This computes the average HOURLY ocean forcing
                    # from the DAILY sum in this AV.
                    # (So division by 24 (ncpl) should
                    # do the same thing as division by bun%cnt).
                    ifields, rfields = self.av_xc2opsum_o.getFields()
                    for r in rfields:
                        temp = self.av_xc2opsum_o.exportRAttr(r)
                        temp /= ncpl
                        self.av_xc2opsum_o.importRAttr(r,temp)
                    #call cpl_bundle_copy(bun_Xc2oPSUM_o,outbun=con_Xc2o%bundle)
                    self.contracts["c2o"].av.copy( self.av_xc2opsum_o )

                    if (not cpl.control.lagocn):
                        self.contracts["o2c"].recv()
                        # call cpl_bundle_mult( con_Xo2c%bundle,
                        #                       bun_areafact_o, "comp2cpl",
                        #                       bunlist=cpl_fields_o2c_fluxes)
                        factor = areafact.av_areafact_o.exportRAttr("comp2cpl")
                        for field in cpl.fields.o2c_fluxes:
                            temp = self.contracts["o2c"].av.exportRAttr(field)
                            temp *= factor
                            self.contracts["o2c"].av.importRAttr(field,temp)
                            
                        self.av_so2c_o.copy(self.contracts["o2c"].av)
                        self.av_fo2c_o.copy(self.contracts["o2c"].av)
                    # --- Start normal interaction with ocn
                    if( cpl.control.lagocn ):
                        print "\nStart of time coordinated integration\n"
                        cpl.control.lagocn = False
                    
                # Recv message from atm
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_a) == 0 ):
                    self.contracts["a2c"].recv()
                    """
                    call cpl_bundle_mult(con_Xa2c%bundle,bun_areafact_a,
                                         'comp2cpl',
                                         bunlist=cpl_fields_a2c_fluxes)
                    call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,
                                          bun_Fa2c_a)
                    cpl_control_fluxAShift = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ashift)
                    """
                    factor = areafact.av_areafact_a.exportRAttr("comp2cpl")
                    ifields,rfields=self.contracts["a2c"].av.getFields()
                    for r in rfields:
                        temp = self.contracts["a2c"].av.exportRAttr(r)
                        temp *= factor
                        self.contracts["a2c"].av.importRAttr(r,temp)
                    self.av_sa2c_a.copy( self.contracts["a2c"].av )
                    self.av_fa2c_a.copy( self.contracts["a2c"].av )
                    cpl.control.fluxashift = self.contracts["a2c"].infobuffer.ibufGet( "ashift" )
                    #---form total rain and snow
                    """
                    call cpl_bundle_zero (bun_precip_a)
                    call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainc')
                    call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainl')
                    call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowc')
                    call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowl')
                    """
                    faxa_rainc = self.av_fa2c_a.exportRAttr("Faxa_rainc")
                    faxa_rainl = self.av_fa2c_a.exportRAttr("Faxa_rainl")
                    faxc_rain = faxa_rainc + faxa_rainl
                    self.av_precip_a.importRAttr("Faxc_rain",faxc_rain)
                    faxa_snowc = self.av_fa2c_a.exportRAttr("Faxa_snowc")
                    faxa_snowl = self.av_fa2c_a.exportRAttr("Faxa_snowl")
                    faxc_snow = faxa_snowc + faxa_snowl
                    self.av_precip_a.importRAttr("Faxc_snow",faxc_snow)

                
                # --- advance date & update control flags ---
                """
                call shr_date_adv1step(date)
                call shr_date_getCDate(date,cDate,sec)
                call shr_date_getYMD(date,year,month,day,sec)
                if (cpl_control_infodbug >= 2 .or. sec == 0) then
                call tStamp_write("cpl",year,month,day,sec)
                endif
                call cpl_control_update(date)
                """
                # --- set infobuf flags --- KLUDGE: need to move this elsewhere
                """
                con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restEOD) = 0
                con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restNow) = 0
                if (cpl_control_restEOD) con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restEOD) = 1
                if (cpl_control_stopEOD) con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_stopEOD) = 1
                if (cpl_control_restNow) con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restNow) = 1
                con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_cdate)   = cDate
                con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_sec)     = sec
                con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_infobug) = cpl_control_infodbug
                con_Xc2l%infobuf = con_Xc2a%infobuf
                con_Xc2i%infobuf = con_Xc2a%infobuf
                con_Xc2o%infobuf = con_Xc2a%infobuf
                call shr_sys_flush(6)
                
                call shr_timer_stop (t26) ; call shr_timer_start(t27)
                """
                # --- create a restart file ---
                """
                call restart_write(date)
                call shr_timer_stop (t27)
                """
                print "cpl: Finished hour %s!"%(n)
                #if (stopeod!=0):
                #    print "CPL STOPPING..."
                #    done = True
                #    break
            # End of NCPL Loop
            if (done):
                break
        print self.name,": end of main integration loop..."
        # Send last message with a stopnow signal:
        for c in ['c2a','c2i','c2l','c2o']:
            self.contracts[c].infobuffer.ibufSet("stopnow",1)
            self.contracts[c].send()
        return

    def merge_atm( self, fix_so2c_a):
        """
        IROUTINE: merge_atm -- merge bundles to form atm input bundle
        DESCRIPTION:
            merge bundles to form atm input bundle
        """
        
        # --- merge atm states & fluxes ---
        # call cpl_bundle_zero(con_Xc2a%bundle) ! zero before adding
        
        # call cpl_bundle_gather(con_Xc2a%bundle, bun_Sa2c_a, bun_Fa2c_a, &
        # &                                       bun_Sl2c_l, bun_Fl2c_l, &
        # &                                       bun_So2c_a, bun_Fo2c_a, &
        # &                                       bun_Si2c_a, bun_Fi2c_a )
        copyables = [self.av_sa2c_a, self.av_fa2c_a,
                     self.av_sl2c_l, self.av_fl2c_l,
                     self.av_so2c_a, self.av_fo2c_a,
                     self.av_si2c_a, self.av_fi2c_a ]
        for c in copyables:
            #print c
            self.contracts["c2a"].av.copy( c )

        #print self.contracts["c2a"].av

        # Factors in frac.frac_a
        ofrac,    n_ofrac     = frac.frac_a.av.exportRAttr("ofrac")
        lfrac,    n_lfrac     = frac.frac_a.av.exportRAttr("lfrac")
        ifrac,    n_ifrac     = frac.frac_a.av.exportRAttr("ifrac")

        Target = "Faxx"
        FieldList = ["taux","tauy","lat","sen","lwup","evap"]
        Sources = ["Faoc","Fall","Faii"]

        for field in FieldList:
            # Target:
            faxx,size_faxx = self.contracts["c2a"].av.av.exportRAttr(Target+"_"+field)
            # Sources in av_aoflux_a
            faoc,size_faoc = self.av_aoflux_a.av.exportRAttr("Faoc"+"_"+field)
            fall,size_faoc = self.av_fl2c_l.av.exportRAttr("Fall"+"_"+field)
            faii,size_faoc = self.av_fi2c_a.av.exportRAttr("Faii"+"_"+field)
            # Sum:
            # Bundle add works like this:
            # First argument is target (A in the equation A = A + B)
            # 2nd Argument is added to the first
            # 3rd->Nth arguments are multiplicative factors applied to the 2nd argument.
            # 1st = 1st + (2nd * 3rd * ... Nth)
            # So this one line is really 3 "cpl_bundle_add" methods together.
            faxx = (faoc * ofrac) + (fall * lfrac) + (faii * ifrac)
            self.contracts["c2a"].av.av.importRAttr(Target+"_"+field, faxx)

        # Tref:
        sx,size_sx = self.contracts["c2a"].av.av.exportRAttr("Sx_tref")
        faoc,size_faoc = self.av_aoflux_a.av.exportRAttr("Faoc_tref")
        sl,size_sl = self.av_sl2c_l.av.exportRAttr("Sl_tref")
        si,size_si = self.av_si2c_a.av.exportRAttr("Si_tref")
        sx = (faoc * ofrac) + (sl * lfrac) + (si * ifrac)
        self.contracts["c2a"].av.av.importRAttr("Sx_tref",sx)
        # Qref:
        sx,size_sx = self.contracts["c2a"].av.av.exportRAttr("Sx_qref")
        faoc,size_faoc = self.av_aoflux_a.av.exportRAttr("Faoc_qref")
        sl,size_sl = self.av_sl2c_l.av.exportRAttr("Sl_qref")
        si,size_si = self.av_si2c_a.av.exportRAttr("Si_qref")
        sx = (faoc * ofrac) + (sl * lfrac) + (si * ifrac)
        self.contracts["c2a"].av.av.importRAttr("Sx_qref",sx)


        # Just Sx fields:
        Target = "Sx"
        FieldList = ["avsdr","anidr","avsdf","anidf"]
        Sources = ["So","Sl","Si"]

        for field in FieldList:
            #target:
            sx,size_sx = self.contracts["c2a"].av.av.exportRAttr(Target+"_"+field)
            # sources:
            so,size_so = self.av_oalbedo_a.av.exportRAttr("So"+"_"+field)
            sl,size_sl = self.av_sl2c_l.av.exportRAttr("Sl"+"_"+field)
            si,size_si = self.av_si2c_a.av.exportRAttr("Si"+"_"+field)

            #Sum:
            sx = (so * ofrac) + (sl * lfrac) + (si * ifrac)
            self.contracts["c2a"].av.av.importRAttr(Target+"_"+field, sx)
            
        # Temperature:        
        sx,size_sx = self.contracts["c2a"].av.av.exportRAttr("Sx_t")
        so,size_so = self.av_so2c_a.av.exportRAttr("So_t")
        sl,size_sl = self.av_sl2c_l.av.exportRAttr("Sl_t")
        si,size_si = self.av_si2c_a.av.exportRAttr("Si_t")
        sx = (so * ofrac) + (sl * lfrac) + (si * ifrac)
        self.contracts["c2a"].av.av.importRAttr("Sx_t",sx)

        # Snow Height:
        sx,size_sx = self.contracts["c2a"].av.av.exportRAttr("Sx_snowh")
        #print self.av_sl2c_l
        # Here "Sl_snowh" is supposed to be used, but since only "Sl_snow"
        # was in the bundle I'm using it instead (MS)
        sl,size_sl = self.av_sl2c_l.av.exportRAttr("Sl_snow")
        sx = (sl * lfrac)
        self.contracts["c2a"].av.av.importRAttr("Sx_t",sx)

        # Kludge:  So_t values have errors due to mapping,
        # can be less then freezing!
        so_t,size_so_t = fix_so2c_a.av.exportRAttr("So_t")
        self.contracts["c2a"].av.av.importRAttr("So_t",so_t)

        # Done with merge_atm:
        return

    def merge_ocn(self):
        """
        !===============================================================================
        !BOP ===========================================================================
        !
        ! !IROUTINE: merge_ocn -- merge bundles to form ocn input bundle
        !
        ! !DESCRIPTION:
        !    merge bundles to form ocn input bundle
        !
        ! !REMARKS:
        !
        ! !REVISION HISTORY:
        !     2002-Jun-06 - B. Kauffman - initial version.
        !
        ! !INTERFACE:  -----------------------------------------------------------------
        """
        
        # call cpl_bundle_zero(bun_Xc2oSNAP_o) ! zero before adding
        # This should default to zero when it's initialized, right?
        # Copy all the fields from the listed bundles into bun_Xc2oSNAP:
        # Not allowed now:
        #self.av_xc2osnap_o.av.Copy( self.av_sl2c_o.av )
        #self.av_xc2osnap_o.av.Copy( self.av_fl2c_o.av )
        # Still allowed:
        self.av_xc2osnap_o.av.Copy( self.av_so2c_o.av )
        self.av_xc2osnap_o.av.Copy( self.av_fo2c_o.av )
        self.av_xc2osnap_o.av.Copy( self.av_si2c_i.av )
        self.av_xc2osnap_o.av.Copy( self.av_fi2c_i.av )
        self.av_xc2osnap_o.av.Copy( self.av_xr2c_o.av )
        self.av_xc2osnap_o.av.Copy( self.av_aoflux_o.av )

        npts = self.av_xc2osnap_o.size()


        foxx_taux = self.av_xc2osnap_o.exportRAttr("Foxx_taux")
        foxx_tauy = self.av_xc2osnap_o.exportRAttr("Foxx_tauy")
        foxx_swnet = self.av_xc2osnap_o.exportRAttr("Foxx_swnet")
        foxx_lat = self.av_xc2osnap_o.exportRAttr("Foxx_lat")
        foxx_sen = self.av_xc2osnap_o.exportRAttr("Foxx_sen")
        foxx_lwup = self.av_xc2osnap_o.exportRAttr("Foxx_lwup")
        foxx_evap = self.av_xc2osnap_o.exportRAttr("Foxx_evap")
        foxx_lwdn = self.av_xc2osnap_o.exportRAttr("Foxx_lwdn")
        foxx_rain = self.av_xc2osnap_o.exportRAttr("Foxx_rain")
        foxx_snow = self.av_xc2osnap_o.exportRAttr("Foxx_snow")
        foxx_prec = self.av_xc2osnap_o.exportRAttr("Foxx_prec")
        foxx_melth = self.av_xc2osnap_o.exportRAttr("Foxx_melth")
        foxx_meltw = self.av_xc2osnap_o.exportRAttr("Foxx_meltw")
        foxx_salt = self.av_xc2osnap_o.exportRAttr("Foxx_salt")
        si_ifrac = self.av_xc2osnap_o.exportRAttr("Si_ifrac")

        faoc_taux = self.av_aoflux_o.exportRAttr("Faoc_taux")
        faoc_tauy = self.av_aoflux_o.exportRAttr("Faoc_tauy")
        faoc_swnet = self.av_aoflux_o.exportRAttr("Faoc_swnet")
        faoc_lat = self.av_aoflux_o.exportRAttr("Faoc_lat")
        faoc_sen = self.av_aoflux_o.exportRAttr("Faoc_sen")
        faoc_lwup = self.av_aoflux_o.exportRAttr("Faoc_lwup")
        faoc_evap = self.av_aoflux_o.exportRAttr("Faoc_evap")

        faxc_rain = self.av_precip_o.exportRAttr("Faxc_rain")
        faxc_snow = self.av_precip_o.exportRAttr("Faxc_snow")
        
        fioi_taux = self.av_fi2c_i.exportRAttr("Fioi_taux")
        fioi_tauy = self.av_fi2c_i.exportRAttr("Fioi_tauy")
        fioi_swpen = self.av_fi2c_i.exportRAttr("Fioi_swpen")
        faxa_lwdn = self.av_fa2c_o.exportRAttr("Faxa_lwdn")
        fioi_melth = self.av_fi2c_i.exportRAttr("Fioi_melth")
        fioi_meltw = self.av_fi2c_i.exportRAttr("Fioi_meltw")
        fioi_salt = self.av_fi2c_i.exportRAttr("Fioi_salt")

        afrac = frac.frac_o.exportRAttr("afrac")
        ifrac = frac.frac_o.exportRAttr("ifrac")
        

        """
        We can do this without the for loop by using the builtin numeric methods:
        """
        # do n = 1,npts
        # for n in range(npts):
        #   afrac = bun_frac_o%data%rAttr(k_afrac,n)
        #   ifrac = bun_frac_o%data%rAttr(k_ifrac,n)
        #   bun_Xc2oSNAP_o%data%rAttr(k_Foxx_taux,n) = &
        #   bun_aoflux_o%data%rAttr(k_Faoc_taux,n)*afrac + &
        #   bun_Fi2c_i%data%rAttr(k_Fioi_taux,n)*ifrac

        try:
            # Suspicion:  this is slower
            foxx_taux = faoc_taux * afrac + fioi_taux * ifrac
            # then:
            #faoc_taux *= afrac
            #fioi_taux *= ifrac
            #foxx_taux += faoc_taux
            #foxx_taux += fioi_taux

            foxx_tauy = faoc_tauy * afrac + fioi_tauy * ifrac

            foxx_swnet = faoc_swnet * afrac + fioi_swpen * ifrac

            foxx_lat = faoc_lat * afrac

            foxx_sen = faoc_sen * afrac 

            foxx_lwup = faoc_lwup * afrac

            foxx_evap = faoc_evap * afrac

            foxx_lwdn = faxa_lwdn * afrac

            foxx_rain = faxc_rain * afrac

            foxx_snow = faxc_snow * afrac

            foxx_prec += foxx_rain
            foxx_prec += foxx_snow

            foxx_melth = fioi_melth * ifrac

            foxx_meltw = fioi_meltw * ifrac

            foxx_salt = fioi_salt * ifrac

            si_ifrac = ifrac
        except ValueError:
            print "len(): foxx_rain -- faxc_rain -- afrac -- npts"
            print "len():",len( foxx_rain ), len(faxc_rain), len(afrac),npts
            raise 
        #    bun_Xc2oSNAP_o%cnt = 1

        # --- Import everything back in ---
        self.av_xc2osnap_o.importRAttr("Foxx_taux",foxx_taux)
        self.av_xc2osnap_o.importRAttr("Foxx_tauy",foxx_tauy)
        self.av_xc2osnap_o.importRAttr("Foxx_swnet",foxx_swnet)
        self.av_xc2osnap_o.importRAttr("Foxx_lat",foxx_lat)
        self.av_xc2osnap_o.importRAttr("Foxx_sen",foxx_sen)
        self.av_xc2osnap_o.importRAttr("Foxx_lwup",foxx_lwup)
        self.av_xc2osnap_o.importRAttr("Foxx_evap",foxx_evap)
        self.av_xc2osnap_o.importRAttr("Foxx_lwdn",foxx_lwdn)
        self.av_xc2osnap_o.importRAttr("Foxx_rain",foxx_rain)
        self.av_xc2osnap_o.importRAttr("Foxx_snow",foxx_snow)
        self.av_xc2osnap_o.importRAttr("Foxx_prec",foxx_prec)
        self.av_xc2osnap_o.importRAttr("Foxx_melth",foxx_melth)
        self.av_xc2osnap_o.importRAttr("Foxx_meltw",foxx_meltw)
        self.av_xc2osnap_o.importRAttr("Foxx_salt",foxx_salt)
        self.av_xc2osnap_o.importRAttr("Si_ifrac",si_ifrac)

        return #end subroutine merge_ocn

def timed():
    print "pyCPL: Start!"
    global TIMESTEP
    # name, doComputations, doMapping, initializeRiverMapping )
    c = Coupler("cpl.nml",True,False,True)
    #c = Coupler(os.path.join("cpl","cpl.nml"),True,False,False)
    starttime = time.time()
    if (cpl.comm.world_pid == cpl.comm.world_pe0):
        ice_temperature = c.contracts["i2c"].av.exportRAttr("Si_t")
        # print ice_temperature
        x = c.contracts["i2c"].domain.ni
        y = c.contracts["i2c"].domain.nj
        if(len(ice_temperature) == (x*y)):
            print "ICE Temperature:"
            if(os.path.exists("ice-t.nc")):
               os.remove("ice-t.nc")
            icetfile = pycdf.CDF("ice-t.nc",(pycdf.NC.WRITE|pycdf.NC.CREATE))
            print "Created NetCDF file..."
            icetfile.definemode()
            icetfile.def_dim("time",pycdf.NC.UNLIMITED)
            icetfile.def_dim("length",len(ice_temperature))
            icetfile.def_dim("x",x)
            icetfile.def_dim("y",y)
            icetfile.def_var("ice_temp",pycdf.NC.DOUBLE,("time","x","y"))
            print "Defined Dims and Vars..."
            icetfile.datamode()
            icetfile.sync()
            print "Synced..."
    c.run(2)
    if (cpl.comm.world_pid == cpl.comm.world_pe0):
        ice_temperature = c.contracts["i2c"].av.exportRAttr("Si_t")
        # print ice_temperature
        x = c.contracts["i2c"].domain.ni
        y = c.contracts["i2c"].domain.nj
        if(len(ice_temperature) == (x*y)):
            #icetfile.var("ice_temp").put(ice_temperature)#[,start=(time,0)])
            icet = Numeric.resize(ice_temperature,(x,y))
            print "TIMESTEP:",TIMESTEP
            icetfile.var("ice_temp").put( icet, start =(TIMESTEP,0,0) )
            icetfile.sync()
            icetfile.close()
            print "NetCDF file written!"
    endtime = time.time()
    cpl.comm.finalize()
    #print "StartDate:",cpl.date
    print "pyCPL: Done!"
    print "Elapsed Time:",endtime-starttime
    
if __name__=="__main__":
    timed()

    

