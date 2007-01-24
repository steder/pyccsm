#!/usr/bin/env python2.4

"""
pyCPL6 Coupler Class

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
import cpl.timer

# Coupler Only Modules
import frac
import flux
import areafact

from cpl.debug import *

debugPrint = cpl.debug.newDebugPrint(True)
mergeocnPrint = cpl.debug.newDebugPrint(False)
avsizePrint = cpl.debug.newDebugPrint(False)
fracPrint = cpl.debug.newDebugPrint(True)
statusPrint = cpl.debug.newDebugPrint(True)
mapPrint = cpl.debug.newDebugPrint(True)
runPrint = cpl.debug.newDebugPrint(True)

TIMESTEP = 0

TIMERS = ["total","run","init","run_init","step"] 
TIMER = cpl.timer.Timer(TIMERS)


abspath = os.path.abspath(__file__)
basedir = os.path.dirname(abspath)
ROOT=basedir


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
    def __init__(self,filepath="cpl.nml",computations=False,mapping=False,river=False, logging=False, logpath="~/exe/pyccsm/"):
        self.name = "cpl"
        self.contracts={}
        self.doComputations = computations
        self.doMapping=mapping
        self.doRiver=river
        self.doLogging=logging
        self.logPath = logpath
        startTime = time.time()
        statusPrint("pyCPL Init: starting @",time.ctime(),"!")
        local_comm = cpl.comm.init( cpl.fields.cplname )
        # print "pyCPL: Barrier 1 @",time.ctime(),": ..."
        # mpi.barrier( mpi.MPI_COMM_WORLD )
        statusPrint("OK!")
        rank = mpi.comm_rank( mpi.MPI_COMM_WORLD )
        size = mpi.comm_size( mpi.MPI_COMM_WORLD )
        statusPrint( "mpi.MPI_COMM_WORLD size=%d,myrank=%d"%(size,rank))
        statusPrint( "Size of Coupler Comm:",mpi.comm_size( local_comm ))
        statusPrint( "Rank in Coupler Comm:",mpi.comm_rank( local_comm ))
        statusPrint( "world_pid:",cpl.comm.world_pid)
        statusPrint( "mph_component_id:",cpl.comm.mph_component_id)
        statusPrint( "component_pid:",cpl.comm.component_pid)
        statusPrint( "component_npe:",cpl.comm.component_npe)
        
        debugPrint( "Python Job running on %s, MPI Node %s"%(os.getenv("HOSTNAME"),rank) )

        self.readNameList( filepath )
        # now we can set these values on all coupler processes.
        self.stepsPerDay = 24
        case_name = self.namelist["case_name"]
        #self.namelist["start_date"]
        self.date = cpl.shr.date.Date()
        self.date.initCDate( int(self.namelist["start_date"]), self.stepsPerDay  )
        self.orb_year  = self.namelist["orb_year"]
        self.flx_epbal = self.namelist["flx_epbal"]
        self.flx_albav = self.namelist["flx_albav"]
        info_dbug = self.namelist["info_dbug"]
        map_a2of_fn = self.namelist["map_a2of_fn"]
        map_a2os_fn = self.namelist["map_a2os_fn"]
        map_o2af_fn = self.namelist["map_o2af_fn"]
        map_r2o_fn = self.namelist["map_r2o_fn"]

        statusPrint("Updating control module with namelist values:")
        cpl.control.casename = case_name
        cpl.control.cdate_a = self.date
        cpl.control.cdate_i = self.date
        cpl.control.cdate_l = self.date
        cpl.control.cdate_o = self.date
        cpl.control.fluxalbav = self.flx_albav
        cpl.control.fluxepbal = self.flx_epbal
        cpl.control.mapfn_a2of = map_a2of_fn
        cpl.control.mapfn_a2os = map_a2os_fn
        cpl.control.mapfn_o2af = map_o2af_fn
        cpl.control.mapfn_r2o = map_r2o_fn

        print "\"%s\""%(self.namelist["start_type"])
        if( self.namelist["start_type"] == "initial" ):
            cpl.control.lagocn = True
            print "*** Setting LAGOCN to TRUE! ***"
        

        # Defaults:
        #cpl.control.ncpl_a = 0
        #cpl.control.ncpl_o = 0
        #cpl.control.ncpl_l = 0
        #cpl.control.ncpl_i = 0
        #cpl.control.ncpl_r = 0

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
        statusPrint("Partial Map Data Init, Map Init")
        self.partialMapInit()
        frac.init( self.map_Fo2a,
                   self.contracts["a2c"].domain,
                   self.contracts["i2c"].domain,
                   self.contracts["Dc2l"].domain,
                   self.contracts["o2c"].domain)

        ##
        if( self.doLogging ):
            statusPrint("Checking LFRAC and AFRAC immediately after frac.init!")
            statusPrint("Outputting LFRAC fields to NetCDF...")
            ni = self.contracts["Dc2l"].domain.ni
            nj = self.contracts["Dc2l"].domain.nj
            n = self.contracts["Dc2l"].domain.n
            if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                frac.frac_l.writeNC(os.path.join(self.logPath,"lfrac-init.nc"),ni,nj,n)
            ni = self.contracts["a2c"].domain.ni
            nj = self.contracts["a2c"].domain.nj
            n = self.contracts["a2c"].domain.n
            if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                frac.frac_a.writeNC(os.path.join(self.logPath,"afrac-init.nc"),ni,nj,n)
        ##
        statusPrint("Send domain info to land model? (optional)")
        #statusPrint("SKIPPING LAND DOMAIN SEND!")
        if( not cpl.control.sendlnddom ):#not cpl.control.sendlnddom ):
            statusPrint(" * lnd DOES NOT request optional domain data exchange")
        else:
            statusPrint(" * lnd requests optional domain data exchange")
            statusPrint( "\t--- initialize gData0...")
            gdata0 = cpl.attributevector.AttributeVector()
            ifields,rfields = self.contracts["Dc2l"].av.getFields()
            gdata0.initialize( ifields, rfields,
                               self.contracts["Dc2l"].av.lsize() )
            statusPrint( "\t--- initialize gData1...")
            gdata1 = cpl.attributevector.AttributeVector()
            ifields,rfields = self.contracts["c2a"].domain.lgrid.getFields()
            gdata1.initialize( ifields, rfields,
                               self.contracts["c2a"].domain.lgrid.lsize() )
            statusPrint( "\t--- initialize gData2...")
            gdata2 = cpl.attributevector.AttributeVector()
            ifields,rfields = [],frac.frac_fields
            gdata2.initialize( ifields, rfields,
                               frac.frac_a.lsize() )
            statusPrint( "\t--- initialized gdata0,gdata1,gdata2! ---")
            
            statusPrint( "\tZeroing \"Dc2l\" AttrVect...")
            #call cpl_bundle_zero(con_Dc2l%bundle)
            self.contracts["Dc2l"].av.zero()

            #statusPrint( "Init Gdata0,1,2"
            
            statusPrint("\t* gData0 gather...")
            gdata0.gather(self.contracts["Dc2l"].av,
                          self.contracts["Dc2l"].domain.gsMap,
                          0, cpl.comm.component_pid,
                          cpl.comm.local_comm )
            mpi.barrier(cpl.comm.local_comm)
            statusPrint("\t* gData1 gather...")
            gdata1.gather(self.contracts["c2a"].domain.lgrid,
                          self.contracts["c2a"].domain.gsMap,
                          0, cpl.comm.component_pid,
                          cpl.comm.local_comm )
            statusPrint("\t* gData2 gather...")
            
            gdata2.gather(frac.frac_a,
                          self.contracts["c2a"].domain.gsMap,
                          0, cpl.comm.component_pid,
                          cpl.comm.local_comm )
            # --- set gData0 to be scattered as new data in con_Dc2l ---
            if (cpl.comm.component_pid == 0): 
                # nfld = cpl_mct_aVect_indexRA(gData2,"lfrac",perrWith="gData2 lfrac")
                # lsize = cpl_mct_aVect_lsize(gData0)
                lfrac2,lfrac2_size = gdata2.exportRAttr("lfrac")
                lsize = gdata0.lsize()
                # do n=1,lsize
                lon,lon_size = gdata1.exportRAttr("lon")
                gdata0.importRAttr("lon",lon)
                lat,lat_size = gdata1.exportRAttr("lat")
                gdata0.importRAttr("lat",lat)
                area,area_size = gdata1.exportRAttr("area")
                gdata0.importRAttr("area",area)
                mask,mask_size = gdata1.exportRAttr("mask")
                gdata0.importRAttr("maska",mask)
                lmask,lmask_size = gdata0.exportRAttr("maskl")
                lfrac0,lfrac0_size = gdata0.exportRAttr("lfrac")

                for n in range(0,lsize-1):
                    if (lfrac2[n] < 1.0e-06):
                        lmask[n] = 0.0
                        lfrac0[n] = 0.0
                        mask[n] = 0.0
                    else:
                        lmask[n] = 1.0
                        lfrac0[n] = lfrac2[n]
                        mask[n] = 1.0
                    gdata0.importRAttr("maskl",lmask)
                    gdata0.importRAttr("lfrac",lfrac0)
                    gdata1.importRAttr("mask",mask)
                                            
            # --- reset dom_l based on dom_a ---
            # --- cpl_fields_grid_mask is not from dom_a (see above loop) ---
            # call cpl_mct_aVect_scatter(gData1,dom_l%lGrid, &
            #   dom_l%gsMap,0,cpl_comm_comp,rcode)

            self.dom_l.lgrid.scatter( gdata1,self.contracts["Dc2l"].domain.gsMap,0, cpl.comm.component_pid, cpl.comm.local_comm )

            # --- scatter gData0 to con_Dc2l bundle ---
            # call cpl_mct_aVect_scatter(gData0,con_Dc2l%bundle%data, &
            #   con_Dc2l%bundle%dom%gsMap,0,cpl_comm_comp,rcode)

            self.contracts["Dc2l"].av.scatter( gdata0,
                                               self.contracts["Dc2l"].domain.gsMap,
                                               0, cpl.comm.component_pid, cpl.comm.local_comm )

            # --- clean up ---
            # if (cpl_comm_comp_pid == 0) then
            #    call cpl_mct_aVect_clean(gData0)
            #    call cpl_mct_aVect_clean(gData1)
            #    call cpl_mct_aVect_clean(gData2)
            # endif

            # --- send con_Dc2l ---
            #call cpl_interface_contractSend(cpl_fields_lndname,con_Dc2l)
            statusPrint("Sending Dc2l Contract...")
            #self.contracts["Dc2l"].infobuffer.send()
            #self.contracts["Dc2l"].av.av.send(self.contracts["Dc2l"].router, 600)
            self.contracts["Dc2l"].send()
            statusPrint("Dc2l Contract Sent!")
        # END OPTIONAL LAND DOMAIN DATA EXCHANGE
        sys.stdout.flush()
        ## LND Contract
        statusPrint("Initializing LND Contracts...")
        statusPrint("Initializing LND: l2c")
        self.contracts["l2c"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.lndname, cpl.control.decomp_a, cpl.fields.l2c_fields)

        sys.stdout.flush()
        statusPrint("Initializing LND: c2l")
        self.contracts["c2l"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.lndname, cpl.control.decomp_a, cpl.fields.c2l_fields)
        statusPrint("Updating land contract domains...")
        self.contracts["l2c"].domain = self.dom_l
        self.contracts["c2l"].domain = self.dom_l
        

        ## Special runoff contract:
        sys.stdout.flush()
        statusPrint("Initializing r2c contract...")
        self.contracts["r2c"] = cpl.contract.ReceivingContract( cpl.fields.cplname, cpl.fields.rtmname, cpl.control.decomp_a, cpl.fields.r2c_fields)

        sys.stdout.flush()
        statusPrint("Init Data_Mod Bundles, Tavg History Bundles...")
        self.bundleInit()
        # Full initialization of mappings.
        self.mapInit()
        
        statusPrint("** Ignoring history bundles for now...")
        # history_avbundleInit()

        # --- START KLUDGE ---
        # fix_So2c_a is similar to bun_So2c_a but with a slightly different mapping
        # fix_frac_a is used to re-normalize map_Fo2a mapping for fix_So2c_a only
        # ---
        # call cpl_bundle_initv(fix_So2c_a,"fix_So2c_a",bun_So2c_a,bun_So2c_a%dom)
        # call cpl_bundle_initv(fix_frac_a,"fix_frac_a",bun_frac_a,bun_frac_a%dom)
        ifields,rfields = self.av_so2c_a.getFields()
        self.fix_so2c_a = cpl.attributevector.AttributeVector(ifields,
                                                              rfields,
                                                              self.av_so2c_a.lsize())
        self.fix_so2c_a.initv( self.av_so2c_a, self.av_so2c_a.lsize() )

        ifields,rfields = frac.frac_a.getFields()
        self.fix_frac_a = cpl.attributevector.AttributeVector(ifields,
                                                              rfields,
                                                              frac.frac_a.lsize())
        debugPrint("frac.frac_a.size():",frac.frac_a.lsize())
        self.fix_frac_a.initv( frac.frac_a, frac.frac_a.lsize() )
        # --- END KLUDGE ---
        
       
        # --- Check atm/lnd and ocn/ice model domains for consistency
        statusPrint("Checking ATM/LND & OCN/ICE domains for consistency...")
        
        # These appear to hang, I don"t think they are
        # being executed on all processors.
        if( False ):
            statusPrint("Comparing a2c - l2c domains: ...")
            self.contracts["a2c"].domain.compare( self.contracts["l2c"].domain,
                                                  enforce_grid=True,
                                                  enforce_area=True,)
            statusPrint("Comparing o2c - i2c domains: ...")
            self.contracts["o2c"].domain.compare( self.contracts["i2c"].domain, enforce_all=True)
        else:
            statusPrint("Currently consistency checks are disabled!")

        # --- Send inital message with cday(cdate?), etc. ---
        debugPrint(" setting integer infobuffer fields...")
        self.infobuf.ibufSet("ncpl",24)
        # getCDate returns (CDate, Seconds), so just the date is [0]
        CDate,sec = self.date.getCDate()
        self.infobuf.ibufSet("cdate",CDate)
        self.infobuf.ibufSet("rcode",0)
        self.infobuf.ibufSet("stopeod",0)
        self.infobuf.ibufSet("stopnow",0)
        self.infobuf.ibufSet("sec",0)
        self.infobuf.ibufSet("infotim",0)
        self.infobuf.ibufSet("infobug",cpl.control.infodbug)
        debugPrint(" setting real infobuffer fields...")
        self.infobuf.rbufSet("spval",cpl.const.spval)
        self.infobuf.rbufSet("eccen",1.670772e-2)
        self.infobuf.rbufSet("obliqr",4.091238e-01)
        self.infobuf.rbufSet("lambm0",-3.250364e-02)
        self.infobuf.rbufSet("mvelpp",4.935568e0)

        statusPrint(" Initial send of infobuffer...")
        # Initial infobuffer send uses tag=1002
        self.infobuf.send( cpl.comm.world_pe0_atm, 1002 )
        self.infobuf.send( cpl.comm.world_pe0_ocn, 1002 )
        self.infobuf.send( cpl.comm.world_pe0_ice, 1002 )
        self.infobuf.send( cpl.comm.world_pe0_lnd, 1002 )

        statusPrint("Verifying acceptable coupling intervals:")
        cpl.control.ncpl_a = self.contracts["a2c"].infobuffer.ibufGet("ncpl")
        cpl.control.ncpl_i = self.contracts["i2c"].infobuffer.ibufGet("ncpl")
        cpl.control.ncpl_l = self.contracts["l2c"].infobuffer.ibufGet("ncpl")
        cpl.control.ncpl_r = self.contracts["r2c"].infobuffer.ibufGet("ncpl")
        cpl.control.ncpl_o = self.contracts["o2c"].infobuffer.ibufGet("ncpl")
        if( (cpl.control.ncpl_a < cpl.control.ncpl_o ) or ( (cpl.control.ncpl_a%cpl.control.ncpl_o)!=0 ) ):
            debugPrint("ERROR: Unacceptable cpl.control.ncpl_a, cpl.control.ncpl_o -",cpl.control.ncpl_a, ",",cpl.control.ncpl_o)
        elif((cpl.control.ncpl_a != cpl.control.ncpl_i) or (cpl.control.ncpl_a != cpl.control.ncpl_l)):
            debugPrint("Error Unacceptable cpl.control.ncpl_[ail] -",cpl.control.ncpl_a,",",cpl.control.ncpl_i,",",cpl.control.ncpl_l)
        else:
            debugPrint("cpl.control.ncpl_[ailro] =",cpl.control.ncpl_a,",",cpl.control.ncpl_i,",",cpl.control.ncpl_l,",",cpl.control.ncpl_r,",",cpl.control.ncpl_o)
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

        statusPrint("Initialize model start date and control flags...")
        cpl.control.init( self.date )
        cpl.control.update(self.date)
        statusPrint("cpl.control.init/update called but are currently null functions!")

        if( True ):
            statusPrint("Resetting bundle domain areas equal to mapping domain areas...")
            print self.map_Fo2a.areasrc
            print self.map_Fo2a.areadst
            src_aream,src_aream_s = self.map_Fo2a.areasrc.exportRAttr("aream")
            dst_aream,dst_aream_s = self.map_Fo2a.areadst.exportRAttr("aream")
            self.contracts["c2a"].domain.lgrid.importRAttr("aream",dst_aream)
            self.contracts["c2l"].domain.lgrid.importRAttr("aream",dst_aream)
            self.contracts["Dc2l"].domain.lgrid.importRAttr("aream",dst_aream)
            self.contracts["c2o"].domain.lgrid.importRAttr("aream",src_aream)
            self.contracts["c2i"].domain.lgrid.importRAttr("aream",src_aream)
            self.contracts["a2c"].domain.lgrid.importRAttr("aream",dst_aream)
            self.contracts["l2c"].domain.lgrid.importRAttr("aream",dst_aream)
            self.contracts["o2c"].domain.lgrid.importRAttr("aream",src_aream)
            self.contracts["i2c"].domain.lgrid.importRAttr("aream",src_aream)
            if( self.doRiver ):
                r_src_aream,r_src_aream_s = self.map_Xr2o.areasrc.exportRAttr("aream")
                self.contracts["r2c"].domain.lgrid.importRAttr("aream",r_src_aream)
        else:
            statusPrint("\t *** Skipping bundle domain reset! ***")
        statusPrint("Initializing Area Correction Data:")
        areafact.init(self.contracts["a2c"].domain,
                      self.contracts["i2c"].domain,
                      self.contracts["l2c"].domain,
                      self.contracts["o2c"].domain,
                      self.contracts["r2c"].domain)
        if( self.doLogging ):
            if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                ni_a = self.contracts["a2c"].domain.ni
                nj_a = self.contracts["a2c"].domain.nj
                ni_l = self.contracts["l2c"].domain.ni
                nj_l = self.contracts["l2c"].domain.nj
                ni_i = self.contracts["i2c"].domain.ni
                nj_i = self.contracts["i2c"].domain.nj
                ni_o = self.contracts["o2c"].domain.ni
                nj_o = self.contracts["o2c"].domain.nj
                ni_r = self.contracts["r2c"].domain.ni
                nj_r = self.contracts["r2c"].domain.nj
                areafact.av_areafact_a.writeNC(os.path.join(self.logPath,"areafact_a.nc"),ni_a,nj_a,ni_a*nj_a)
                areafact.av_areafact_l.writeNC(os.path.join(self.logPath,"areafact_l.nc"),ni_l,nj_l,ni_l*nj_l)
                areafact.av_areafact_i.writeNC(os.path.join(self.logPath,"areafact_i.nc"),ni_i,nj_i,ni_i*nj_i)
                areafact.av_areafact_o.writeNC(os.path.join(self.logPath,"areafact_o.nc"),ni_o,nj_o,ni_o*nj_o)
                if(self.doRiver):
                    areafact.av_areafact_r.writeNC(os.path.join(self.logPath,"areafact_r.nc"),ni_r,nj_r,ni_r*nj_r)
        # Initial receive from model(s)
        statusPrint(" Initial receive from models...")
        self.contracts["a2c"].recv()
        self.contracts["a2c"].av.mult( areafact.av_areafact_a, "comp2cpl",
                                       cpl.fields.a2c_fluxes )
        cpl.control.icdata_a = False
        if self.contracts["a2c"].infobuffer.ibufGet("userest") != 0:
            cpl.control.icdata_a = True
            
        self.contracts["i2c"].recv()
        self.contracts["i2c"].av.mult( areafact.av_areafact_i, "comp2cpl",
                                       cpl.fields.i2c_fluxes )
        self.contracts["i2c"].av = flux.fluxAlbi( self.date, self.contracts["i2c"].av, self.contracts["i2c"].domain )
        cpl.control.icdata_i = False
        if self.contracts["i2c"].infobuffer.ibufGet("userest") != 0:
            cpl.control.icdata_i = True
        self.contracts["l2c"].recv()
        self.contracts["l2c"].av.mult( areafact.av_areafact_l, "comp2cpl",
                                       cpl.fields.l2c_fluxes )
        cpl.control.icdata_l = False
        if self.contracts["l2c"].infobuffer.ibufGet("userest") != 0:
            cpl.control.icdata_l = True
        self.contracts["r2c"].recv()
        self.contracts["r2c"].av.mult( areafact.av_areafact_r, "comp2cpl",
                                       cpl.fields.r2c_fluxes )
        cpl.control.icdata_r = False
        if self.contracts["r2c"].infobuffer.ibufGet("userest") != 0:
            cpl.control.icdata_r = True
        self.contracts["o2c"].recv()
        self.contracts["o2c"].av.mult( areafact.av_areafact_o, "comp2cpl",
                                       cpl.fields.o2c_fluxes )
        self.av_oalbedo_o = flux.fluxAlbo( self.date, self.av_oalbedo_o, self.contracts["o2c"].domain )
        cpl.control.icdata_o = False
        if self.contracts["o2c"].infobuffer.ibufGet("userest") != 0:
            cpl.control.icdata_o = True
        
        # Read IC data from restart file?
        statusPrint("Skipping restart IC read!")

        # PROCESS Initial Condition DATA
        statusPrint("**  Processing Initial Condition Data:")
        statusPrint("* Processing ATM IC data...")
        #    call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,bun_Fa2c_a)
        self.av_sa2c_a.copy( self.contracts["a2c"].av )
        self.av_fa2c_a.copy( self.contracts["a2c"].av )
        #    call cpl_bundle_zero (bun_precip_a)
        self.av_precip_a.zero()
        #    call cpl_bundle_add(bun_precip_a,"Faxc_rain",
        #                        bun_Fa2c_a,"Faxa_rainc")
        #    call cpl_bundle_add(bun_precip_a,"Faxc_rain",
        #                        bun_Fa2c_a,"Faxa_rainl")
        #    call cpl_bundle_add(bun_precip_a,"Faxc_snow",
        #                        bun_Fa2c_a,"Faxa_snowc")
        #    call cpl_bundle_add(bun_precip_a,"Faxc_snow",
        #                        bun_Fa2c_a,"Faxa_snowl")
        faxa_rainc,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_rainc")
        faxa_rainl,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_rainl")
        faxa_snowc,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_snowc")
        faxa_snowl,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_snowl")
        precip_a = faxa_rainc + faxa_rainl 
        self.av_precip_a.importRAttr( "Faxc_rain" , precip_a  )
        precip_a = faxa_snowc + faxa_snowl
        self.av_precip_a.importRAttr( "Faxc_snow" , precip_a )

        cpl.control.fluxashift = self.contracts["a2c"].infobuffer.ibufGet("ashift")

        statusPrint("* Processing ICE IC data...")
        lsize = self.contracts["i2c"].av.lsize()
        ifrac_i, size = self.contracts["i2c"].av.exportRAttr("Si_ifrac")
        frac.set( ifrac_i, self.map_Fo2a,
                       self.contracts["a2c"].domain,
                       self.contracts["i2c"].domain,
                       self.contracts["l2c"].domain,
                       self.contracts["o2c"].domain )
        
        #    call cpl_bundle_split(con_Xi2c%bundle,bun_Si2c_i,bun_Fi2c_i)
        mapPrint("+ ICE IC prepping avs:\n")
        self.av_si2c_i.copy(self.contracts["i2c"].av)
        self.av_fi2c_i.copy(self.contracts["i2c"].av)
        #    call cpl_map_bun(bun_Si2c_i,bun_Si2c_a,map_So2a, &
        #                     bun_frac_i,"ifrac",bun_frac_a,"ifrac")
        mapPrint("+ ICE IC mapping So2a:\n")
        if( self.doLogging ):
            if ( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                ni = self.contracts["i2c"].domain.ni
                nj = self.contracts["i2c"].domain.nj
                self.av_si2c_i.writeNC(os.path.join(self.logPath,"ice-before-ic-mapping.nc"),ni,nj, (ni*nj) )
                
        self.av_si2c_a = self.map_So2a.mapAV( self.av_si2c_i,
                                              frac.frac_i,"ifrac",
                                              frac.frac_a,"ifrac")
        
        # BEGIN: ICE IC Mapping Fo2a
        mapPrint("+ICE IC mapping Fo2a:\n")
        self.av_fi2c_a = self.map_Fo2a.mapAV( self.av_fi2c_i,
                                              frac.frac_i,"ifrac",
                                              frac.frac_a,"ifrac")
        # END: ICE IC Mapping Fo2a
        
        statusPrint("* Processing LND IC data...")
        #    call cpl_bundle_split(con_Xl2c%bundle,bun_Sl2c_l,bun_Fl2c_l )
        self.av_sl2c_l.copy( self.contracts["l2c"].av )
        self.av_fl2c_l.copy( self.contracts["l2c"].av )
        #  map land fields to ocean, not allowed now.
        # !  call cpl_map_bun(bun_Sl2c_l,bun_Sl2c_o,map_Sa2o)
        statusPrint("- Mapping land fields to ocean, but this is not allowed?")
        # self.av_sl2c_o = self.map_Sa2o.mapAV( self.av_sl2c_l )
        # !  call cpl_map_bun(bun_Fl2c_l,bun_Fl2c_o,map_Fa2o)
        # self.av_fl2c_o = self.map_Fa2o.mapAV( self.av_fl2c_l )
        #    !--- map land fields to ocean, not allowed now.
        statusPrint("- Mapping land fields to ocean, skipped!")
        statusPrint("Outputting LND fields to NetCDF...")
        """
        Let's check the atmosphere contracts fields: 
        """
        if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
            ni = self.contracts["l2c"].domain.ni
            nj = self.contracts["l2c"].domain.nj
            n = self.contracts["l2c"].domain.n
            self.contracts["l2c"].av.writeNC(os.path.join(self.logPath,"land-ic.nc"),ni,nj,n)
            frac.frac_l.writeNC(os.path.join(self.logPath,"lfrac.nc"),ni,nj,n)
            frac.frac_a.writeNC(os.path.join(self.logPath,"afrac.nc"),ni,nj,n)
        
        statusPrint("Processing OCN IC data...")
        #    call cpl_bundle_split(con_Xo2c%bundle,bun_So2c_o,bun_Fo2c_o)
        self.av_so2c_o.copy( self.contracts["o2c"].av )
        self.av_fo2c_o.copy( self.contracts["o2c"].av )
        #    call cpl_map_bun(bun_So2c_o,bun_So2c_a,map_So2a, &
        #                     bun_frac_o,"afrac",bun_frac_a,"ofrac")
        self.av_so2c_a = self.map_So2a.mapAV( self.av_so2c_o,
                                              frac.frac_o, "afrac",
                                              frac.frac_a, "ofrac")
        #    !*** KLUDGE - start *********
        #    call cpl_map_bun(bun_frac_o,fix_frac_a,map_So2a)
        self.fix_frac_a = self.map_So2a.mapAV( frac.frac_o ) 
        #    call cpl_map_bun(bun_So2c_o,fix_So2c_a,map_So2a, &
        #                     bun_frac_o,"afrac",fix_frac_a,"afrac")
        self.fix_so2c_a = self.map_So2a.mapAV( self.av_so2c_o,
                                              frac.frac_o, "afrac",
                                              self.fix_frac_a, "afrac")
        #    !*** KLUDGE - end ***********
        #    call cpl_map_bun(bun_Fo2c_o,bun_Fo2c_a,map_Fo2a, &
        #                     bun_frac_o,"afrac",bun_frac_a,"ofrac")
        self.av_fo2c_a = self.map_Fo2a.mapAV( self.av_fo2c_o,
                                              frac.frac_o,"afrac",
                                              frac.frac_a,"ofrac")
        #    call cpl_map_bun(bun_aoflux_o,bun_aoflux_a,map_Fo2a, &
        #                     bun_frac_o,"afrac",bun_frac_a,"ofrac")
        self.av_aoflux_a = self.map_Fo2a.mapAV( self.av_aoflux_o,
                                              frac.frac_o,"afrac",
                                              frac.frac_a,"ofrac")
        #    call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
        #                     bun_frac_o,"afrac",bun_frac_a,"ofrac")
        self.av_oalbedo_a = self.map_So2a.mapAV( self.av_oalbedo_o,
                                              frac.frac_o,"afrac",
                                              frac.frac_a,"ofrac")

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
            self.merge_atm( self.fix_so2c_a )

            statusPrint("Outputting ATM fields to NetCDF...")
            """
            Let"s check the atmosphere contracts fields: 
            """
            
            if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                ni = self.contracts["c2a"].domain.ni
                nj = self.contracts["c2a"].domain.nj
                n = self.contracts["c2a"].domain.n
                self.contracts["c2a"].av.writeNC(os.path.join(self.logPath,"atm-ic.nc"),ni,nj,n)
                

            statusPrint("\tSend albedos to ATM, recv new ATM IC\"s")
            #     call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,
            #                       "cpl2comp",bunlist=cpl_fields_c2a_fluxes)
            self.contracts["c2a"].av.extremes("c2a")
            self.contracts["c2a"].av.mult( areafact.av_areafact_a,
                                           "cpl2comp",
                                           cpl.fields.c2a_fluxes )
                
            #     call cpl_interface_contractSend(cpl_fields_atmname,con_Xc2a)
            self.contracts["c2a"].av.nanCheck("c2a bundle")
            
            self.contracts["c2a"].send()
            #     call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,"comp2cpl",  &
            #               bunlist=cpl_fields_c2a_fluxes)
            self.contracts["c2a"].av.mult( areafact.av_areafact_a,
                                           "comp2cpl",
                                           cpl.fields.c2a_fluxes )
            statusPrint("\twait for new atm IC data...")
            # call cpl_interface_contractRecv(cpl_fields_atmname,con_Xa2c)
            cpl.contract.debugPrint = cpl.debug.newDebugPrint(True)
            self.contracts["a2c"].recv()
            self.contracts["a2c"].av.nanCheck("a2c bundle")
            if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                ni = self.contracts["a2c"].domain.ni
                nj = self.contracts["a2c"].domain.nj
                n = ni * nj
                self.contracts["a2c"].av.writeNC(os.path.join(self.logPath,"atm-ic-recvd.nc"),ni,nj,n)
            cpl.contract.debugPrint = cpl.debug.newDebugPrint(False)
            statusPrint(self.contracts["a2c"].av.getFields())
            statusPrint(cpl.fields.a2c_fluxes)
            # call cpl_bundle_mult(con_Xa2c%bundle,bun_areafact_a,"comp2cpl",  &
            #               bunlist=cpl_fields_a2c_fluxes)
            self.contracts["a2c"].av.mult( areafact.av_areafact_a,
                                           "comp2cpl",
                                           cpl.fields.a2c_fluxes )
            if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                ni = self.contracts["a2c"].domain.ni
                nj = self.contracts["a2c"].domain.nj
                n = ni * nj
                self.contracts["a2c"].av.writeNC(os.path.join(self.logPath,"atm-ic-recvd-mult.nc"),ni,nj,n)
            statusPrint("\tprocess atm IC data...")
            # call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,bun_Fa2c_a)
            self.av_sa2c_a.copy( self.contracts["a2c"].av )
            self.av_fa2c_a.copy( self.contracts["a2c"].av )
            # call cpl_bundle_zero (bun_precip_a)
            self.av_precip_a.zero()
            # call cpl_bundle_add(bun_precip_a,"Faxc_rain",bun_Fa2c_a,"Faxa_rainc")
            # call cpl_bundle_add(bun_precip_a,"Faxc_rain",bun_Fa2c_a,"Faxa_rainl")
            # call cpl_bundle_add(bun_precip_a,"Faxc_snow",bun_Fa2c_a,"Faxa_snowc")
            # call cpl_bundle_add(bun_precip_a,"Faxc_snow",bun_Fa2c_a,"Faxa_snowl")
            faxa_rainc,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_rainc")
            faxa_rainl,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_rainl")
            faxa_snowc,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_snowc")
            faxa_snowl,tmpsize = self.av_fa2c_a.exportRAttr( "Faxa_snowl")
            precip_a = faxa_rainc + faxa_rainl 
            self.av_precip_a.importRAttr( "Faxc_rain" , precip_a  )
            precip_a = faxa_snowc + faxa_snowl
            self.av_precip_a.importRAttr( "Faxc_snow" , precip_a )
            
            ### 
            # cpl_control_fluxAShift = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ashift)
            cpl.control.fluxashift = self.contracts["a2c"].infobuffer.ibufGet("ashift")
        # End Optional ATM Initialization"

        statusPrint("Create data as necessary for 1st iteration of main event loop")
        # --- map ocn & ice albedos onto atm domain ---
        # call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
        #            bun_frac_o,"afrac",bun_frac_a,"ofrac")
        self.av_oalbedo_a = self.map_So2a.mapAV( self.av_oalbedo_o,
                                                 frac.frac_o,"afrac",
                                                 frac.frac_a,"ofrac")
        # call cpl_map_bun(con_Xr2c%bundle,bun_Xr2c_o,map_Xr2o)
        if self.doRiver:
            self.av_xr2c_o = self.map_Xr2o.mapAV( self.contracts["r2c"].av )

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

    def readNameList(self, filepath ):
        if( cpl.comm.component_pid == 0 ):
            # fortran: call shr_msg_dirio("cpl")
            statusPrint( "reading input namelist file:")
            #cpl.control.init()
            self.namelist = cpl.control.readNamelist( filepath )
            statusPrint( self.namelist)
        debugPrint( "Broadcasting namelist to all coupler processors:" )
        ## Broadcast namelist to all processors of this
        # component.  (We"ll serialize the namelist dict first)
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

    def partialMapInit(self):
        """
        Simply initializes the Fo2a mapping:
        """
        map_o2af_fn = self.namelist["map_o2af_fn"]
        debugPrint(map_o2af_fn)
        statusPrint("Partial mapping init...")
        statusPrint("\tInitializing map_Fo2a:...")
        self.map_Fo2a = cpl.map.Map()
        self.map_Fo2a.initialize( self.contracts["o2c"].domain,
                                  self.contracts["a2c"].domain,
                                  "Fo2a",
                                  os.path.join(ROOT,"data",map_o2af_fn),
                                  "dst" )
        statusPrint("\tInitializing map_Fo2a:OK")
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
                                  os.path.join(ROOT,"data",map_a2of_fn),
                                  "src")
        statusPrint("\tInitializing map_Sa2o:OK")
        
        statusPrint("\tInitializing map_Fa2o:...")
        self.map_Fa2o = cpl.map.Map()
        self.map_Fa2o.initialize( self.contracts["a2c"].domain,
                                  self.contracts["o2c"].domain,
                                  "Fa2o",
                                  os.path.join(ROOT,"data",map_a2of_fn),
                                  "src")
        statusPrint("\tInitializing map_Fa2o:OK")

        statusPrint("\tInitializing map_So2a:...")
        cpl.map.initPrint = cpl.debug.newDebugPrint(True)
        cpl.map.testPrint = cpl.debug.newDebugPrint(True)
        
        self.map_So2a = cpl.map.Map()
        mapPrint("map_So2a details:")
        mapPrint("* self.contracts[\"o2c\"].domain.lsize():",self.contracts["o2c"].domain.lsize())
        mapPrint("* self.contracts[\"a2c\"].domain.lsize():",self.contracts["a2c"].domain.lsize())
        mapPrint("* datafile:",os.path.join(ROOT,"data",map_o2af_fn))
        self.map_So2a.initialize( self.contracts["o2c"].domain,
                                  self.contracts["a2c"].domain,
                                  "So2a",
                                  os.path.join(ROOT,"data",map_o2af_fn),
                                  "dst")
        cpl.map.initPrint = cpl.debug.newDebugPrint(False)
        cpl.map.testPrint = cpl.debug.newDebugPrint(False)
        statusPrint("\tInitializing map_So2a:OK")

        statusPrint("\tInitializing map_Fo2a:...")
        self.map_Fo2a = cpl.map.Map()
        self.map_Fo2a.initialize( self.contracts["o2c"].domain,
                                  self.contracts["a2c"].domain,
                                  "Fo2a",
                                  os.path.join(ROOT,"data",map_o2af_fn),
                                  "dst" )
        statusPrint("\tInitializing map_Fo2a:OK")

        if (self.doRiver):
            statusPrint("\tInitializing map_Xr2o:...")
            self.map_Xr2o = cpl.map.Map()
            self.map_Xr2o.initialize( self.contracts["r2c"].domain,
                                      self.contracts["o2c"].domain,
                                      "Xr2o",
                                      os.path.join(ROOT,"data",map_r2o_fn),
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
        avsizePrint("size o avs:",size_o)
        avsizePrint("size a avs:",size_a)
        avsizePrint("size l avs:",size_l)
        avsizePrint("size r avs:",size_r)
        avsizePrint("size i avs:",size_i)

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
        self.av_xc2osnap_o.initialize( [], cpl.fields.c2o_fields, size_o )
        self.av_xc2opsum_o.initialize( [], cpl.fields.c2o_fields, size_o )

        aoflux_fields = flux.AOFLUX_FIELDS
        self.av_aoflux_o.initialize([], aoflux_fields, size_o )
        self.av_aoflux_a.initialize([], aoflux_fields, size_a )
        albedo_fields = ["So_avsdr","So_anidr","So_avsdf","So_anidf"]
        self.av_oalbedo_o.initialize([], albedo_fields, size_o )
        self.av_oalbedo_a.initialize([], albedo_fields, size_a )
        precip_fields = ["Faxc_rain","Faxc_snow"]
        self.av_precip_o.initialize([],precip_fields,size_o)
        self.av_precip_a.initialize([],precip_fields,size_a)
        #
        self.av_sa2c_a.initialize([],cpl.fields.a2c_states,size_a)
        self.av_fa2c_a.initialize([],cpl.fields.a2c_fluxes,size_a)
        self.av_sa2c_o.initialize([],cpl.fields.a2c_states,size_o)
        self.av_fa2c_o.initialize([],cpl.fields.a2c_fluxes,size_o)
        #
        self.av_sl2c_l.initialize([],cpl.fields.l2c_states,size_l)
        self.av_fl2c_l.initialize([],cpl.fields.l2c_fluxes,size_l)
        #self.av_sl2c_o.initialize([],cpl.fields.l2c_states,size_o)
        #self.av_fl2c_o.initialize([],cpl.fields.l2c_fluxes,size_o)
        #
        self.av_xr2c_o.initialize([],cpl.fields.r2c_fields,size_o)
        #
        self.av_so2c_o.initialize([],cpl.fields.o2c_states,size_o)
        self.av_fo2c_o.initialize([],cpl.fields.o2c_fluxes,size_o)
        self.av_so2c_a.initialize([],cpl.fields.o2c_states,size_a)
        self.av_fo2c_a.initialize([],cpl.fields.o2c_fluxes,size_a)
        #
        self.av_si2c_i.initialize([],cpl.fields.i2c_states,size_i)
        self.av_fi2c_i.initialize([],cpl.fields.i2c_fluxes,size_i)
        self.av_si2c_a.initialize([],cpl.fields.i2c_states,size_a)
        self.av_fi2c_a.initialize([],cpl.fields.i2c_fluxes,size_a)
        return

        
    def run( self, days ):
        """
        Timestep the model for "days" days.

        (Maybe we should seperate out all the "one time"
        code in here and put it into the init method.

        Then we can consider making just a simple
        "step" routine that contains all the code that
        needs to occur daily.  That way this method ("run")
        could be simplified down to simply a

        for i in days:
            self.step()
        )
        """
        TIMER.start("run")
        # Main Integration Loop: repeat until cpl says stop:
        statusPrint(self.name,": beginning main integration loop...")
        statusPrint(self.name," Running for ",days," days.")
        debugPrint("NCPL:\tVALUE")
        debugPrint("atm:\t",cpl.control.ncpl_a)
        debugPrint("ocn:\t",cpl.control.ncpl_o)
        debugPrint("ice:\t",cpl.control.ncpl_i)
        debugPrint("lnd:\t",cpl.control.ncpl_l)
        debugPrint("rtm:\t",cpl.control.ncpl_r)
        iterations=0 # outer loop index
        done = False
        stopnow = 0
        ncpl = self.infobuf.ibufGet("ncpl")
        while not cpl.control.stopeod:
            debugPrint("cpl: Top of outer loop:")
            iterations+=1
            if (iterations >= days):
                debugPrint("Setting Stop Now!")
                cpl.control.stopeod=1
                stopeod = 1
                for c in ["c2a","c2i","c2o","c2l","r2c"]:
                    self.contracts[c].infobuffer.ibufSet("stopeod",stopeod)
            # Beginning of NCPL Loop
            # 1 to cpl.control.ncpl+1(excludes end, so 1-24)
            
            for n in xrange(1,cpl.control.ncpl_a+1):
                TIMER.start("step")
                debugPrint( "cpl: Starting hour %s, %s..."%(n+(cpl.control.ncpl_a*(iterations-1)),cpl.control.ncpl_a*days) )
                a = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_a)
                o = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_o)
                i = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_i)
                l = (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_l)
                r = (n)%(cpl.control.ncpl_a/cpl.control.ncpl_r)
                debugPrint( "a,o,i,l,r = ",a,o,i,l,r )
                
                # Send message to ocn
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_o) == 0 ):
                    if (not cpl.control.lagocn):
                        runPrint("Sending OCN (c2o)...")
                        """ Use av.mult and av.zero to clean this up! """
                        self.contracts["c2o"].av.mult( areafact.av_areafact_o,
                                                       "cpl2comp", cpl.fields.c2o_fluxes)
                        self.contracts["c2o"].send()
                        self.contracts["c2o"].av.mult( areafact.av_areafact_o,
                                                       "comp2cpl", cpl.fields.c2o_fluxes)
                    # Zero-out partial sum:    
                    self.av_xc2opsum_o.zero()
                                            
                # Send message to land
                """
                call cpl_bundle_gather(con_Xc2l%bundle, bun_Sa2c_a, bun_Fa2c_a,
                &                      bun_Sl2c_l, bun_Fl2c_l,
                &                      bun_So2c_a, bun_Fo2c_a,
                &                      bun_Si2c_a, bun_Fi2c_a  ) 
                call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,"cpl2comp",
                                 bunlist=cpl_fields_c2l_fluxes)
                call cpl_interface_contractSend(cpl_fields_lndname,con_Xc2l)
                call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,"comp2cpl",
                                     bunlist=cpl_fields_c2l_fluxes)
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_l) == 0 ):
                    runPrint("Sending LND (c2l) contract...")
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
                    self.contracts["c2l"].av.mult( areafact.av_areafact_l,
                                                   "cpl2comp", cpl.fields.c2l_fluxes )
                    #self.contracts["c2l"].av.extremes("c2l in run method")
                    #self.contracts["c2l"].av.nanCheck("c2l in run method")
                    # Output data and check for errors:
                    if(self.doLogging):
                        if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                            ni = self.contracts["c2l"].domain.ni
                            nj = self.contracts["c2l"].domain.nj
                            self.contracts["c2l"].av.writeNC(os.path.join(self.logPath,"c2l-run.nc"),ni,nj,(ni*nj))
                    # Send
                    debugPrint("** calling self.contracts[\"c2l\"].send()... **")
                    #self.contracts["Dc2l"].send()
                    self.contracts["c2l"].send()
                    debugPrint( "** returned from self.contracts[\"c2l\"].send()! **" )
                    # Mult by comp2cpl
                    self.contracts["c2l"].av.mult( areafact.av_areafact_l,
                                                   "comp2cpl", cpl.fields.c2l_fluxes )
                        
                statusPrint( "** Compute, map, merge ice inputs" )
                
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
                    if (cpl.control.fluxepbal != "off"): 
                        # cpl_control_fluxepfac = con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_precAdj)
                        cpl.control.fluxepfac = self.contracts["o2c"].infobuffer.ibufGet("precadj")
                        # cpl_control_fluxEPfac = cpl_control_fluxEPfac * 1.0e-6_r8
                        cpl.control.fluxepfac = cpl.control.fluxepfac * 1.0e-6
                        # call flux_epbal(date,bun_aoflux_o,con_Xi2c%bundle, 
                        #                bun_precip_o,bun_Xr2c_o,bun_frac_o)
                        result_epbal = flux.fluxEpbal(self.date,
                                       self.av_aoflux_o,
                                       self.contracts["i2c"].av,
                                       self.av_precip_o,
                                       self.contracts["r2c"].av,
                                       frac.frac_o)
                        if( result_epbal != None ):
                            self.av_precip_o.copy( result_epbal[0] )
                            self.contracts["r2c"].av.copy( result_epbal[1] )
                    else:
                        statusPrint("*** FluxEPBal calculation is turned off!")
                # Merge total snow and precip for ice input
                if (self.doComputations):
                    statusPrint("computing ice input by merging total snow and precip:")
                    # call cpl_bundle_zero (con_Xc2i%bundle,"Faxc_rain")
                    # call cpl_bundle_zero (con_Xc2i%bundle,"Faxc_snow")
                    """ Zeroing is unnecessary, we'll just overwrite the old data """

                    # call cpl_bundle_copy(bun_precip_o,bunrList="Faxc_rain",&
                    #                     bunTrList="Faxc_rain",outbun=con_Xc2i%bundle)
                    # call cpl_bundle_copy(bun_precip_o,bunrList="Faxc_snow",&
                    #                     bunTrList="Faxc_snow",outbun=con_Xc2i%bundle)
                    statusPrint("av_precip_o size:",self.av_precip_o.size())
                    statusPrint("ice av size:",self.contracts["c2i"].av.size())
                    faxc_rain,faxc_rain_size = self.av_precip_o.exportRAttr("Faxc_rain")
                    faxc_snow,faxc_snow_size = self.av_precip_o.exportRAttr("Faxc_snow")
                    self.contracts["c2i"].av.importRAttr("Faxc_rain",faxc_rain)
                    self.contracts["c2i"].av.importRAttr("Faxc_snow",faxc_snow)
                    # self.contracts["c2i"].av.copy(self.av_precip_o)
                    """
                    Rather then the av.copy( ) function we"ll just import
                    the interesting fields manually.  I"m concerned that
                    the copy might clobber fields other then Faxc_rain/snow.
                    """
                    statusPrint("av_precip_o copied into self.contracts[\"c2i\"].av!")
                    
                # Correct a->o vector mapping near North Pole(NP)

                if( True ):
                    statusPrint( "*** Skipping NPFix for now!" )
                    # call cpl_map_npfix(bun_Sa2c_a,bun_Sa2c_o,"Sa_u","Sa_v")
                else:
                    statusPrint( "*** Calling NPFix ! ***" )
                    # Note this doesn't overwrite the real av_sa2c_o yet
                    av_sa2c_o = cpl.map.npFix( self.av_sa2c_a, self.av_sa2c_o, 
                                               self.contracts["c2a"].domain,
                                               self.contracts["c2o"].domain,
                                               "Sa_u","Sa_v" )
                # Send message to ice
                """
                call cpl_bundle_gather(con_Xc2i%bundle, bun_Sa2c_o, bun_Fa2c_o, &
!--- not allowed now ------------------------------ bun_Sl2c_o, bun_Fl2c_o, &
            &                                       bun_So2c_o, bun_Fo2c_o, &
            &                                       bun_Si2c_i, bun_Fi2c_i, &
            &                                       fcopy = .true.  )
            call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,"cpl2comp", &
                                 bunlist=cpl_fields_c2i_fluxes)
            call cpl_interface_contractSend(cpl_fields_icename,con_Xc2i)
            call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,"comp2cpl", &
                                 bunlist=cpl_fields_c2i_fluxes)
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_i) == 0 ):
                    runPrint("Sending ICE(c2i) contract...")
                    # Gather:
                    self.contracts["c2i"].av.copy( self.av_sa2c_o )
                    self.contracts["c2i"].av.copy( self.av_fa2c_o )
                    self.contracts["c2i"].av.copy( self.av_so2c_o )
                    self.contracts["c2i"].av.copy( self.av_fo2c_o )
                    self.contracts["c2i"].av.copy( self.av_si2c_i )
                    self.contracts["c2i"].av.copy( self.av_fi2c_i )
                    # Mult, send, mult
                    self.contracts["c2i"].av.mult( areafact.av_areafact_i,"cpl2comp",
                                                   cpl.fields.c2i_fluxes )
                    # Send:
                    self.contracts["c2i"].send()
                    # Mult by comp2cpl
                    self.contracts["c2i"].av.mult( areafact.av_areafact_i,"comp2cpl",
                                                   cpl.fields.c2i_fluxes )
                # Compute net solar flux into ocn
                # CALL flux_solar( bun_Fa2c_o, bun_oalbedo_o, bun_aoflux_o )
                statusPrint("Calling flux_solar to compute net solar flux into OCN ...")
                self.av_aoflux_o = flux.fluxSolar( self.av_fa2c_o, self.av_oalbedo_o, self.av_aoflux_o )
                # Merge ocn inputs
                self.merge_ocn()
                
                # Form partial sum of tavg ocn inputs (virtual "send" to ocn)
                # CALL cpl_bundle_accum( bun_Xc2oSNAP_o, outbun=bun_Xc2oPSUM_o)
                
                # Write bit-for-bit check info
                # do diagnositcs
                # Skipping Diagnostics and B4B check!
                
                # Compute ocn albedos (virtual "recv" from ocn)
                # CALL flux_albo(date, bun_oalbedo_o)
                """
                Added the ocean domain to this call because av"s don"t contain domain info like
                bundles do.  This is one of the few instances where the domain reference in the bundle
                is actually useful.
                """
                self.av_oalbedo_o = flux.fluxAlbo( self.date, self.av_oalbedo_o, self.contracts["o2c"].domain )
                # Compute atm/ocn fluxes
                # CALL flux_atmOcn(con_Xo2c%bundle, bun_Sa2c_o, cpl_control_dead_ao, bun_aoflux_o )
                self.av_aoflux_o = flux.fluxAtmOcn( self.contracts["o2c"].av, self.contracts["o2c"].domain,
                                                    self.av_sa2c_o, self.contracts["o2c"].domain )
                
                # Recv msg from ice
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_i) == 0 ):
                    runPrint("Receiving data from ICE(c2i) contract...")
                    self.contracts["i2c"].recv()
                    # Mult:
                    ifields,rfields=self.contracts["c2i"].av.getFields()
                    self.contracts["c2i"].av.mult( areafact.av_areafact_i, "comp2cpl",
                                                   cpl.fields.c2i_fluxes)
                    # update surface fracs wrt new ice frac
                    ifrac_i,size = self.contracts["i2c"].av.exportRAttr("Si_ifrac")
                    frac.set(ifrac_i, self.map_Fo2a,
                                  self.contracts["a2c"].domain,
                                  self.contracts["i2c"].domain,
                                  self.contracts["l2c"].domain,
                                  self.contracts["o2c"].domain)

                
                # Add diurnal cycle to ice albedos
                # CALL flux_albi( date, con_Xi2c%bundle )
                self.contracts["i2c"].av = flux.fluxAlbi( self.date, self.contracts["i2c"].av, self.contracts["i2c"].domain )
                # CALL cpl_bundle_split( con_Xi2c%bundle, bun_Si2c_i,
                #                        bun_Fi2c_i )
                self.av_si2c_i.copy( self.contracts["i2c"].av )
                self.av_fi2c_i.copy( self.contracts["i2c"].av )
                
                # Map ocn states/fluxes to atm -- BEGIN
                self.av_si2c_a = self.map_So2a.mapAV( self.av_si2c_i,
                                                      frac.frac_i,"ifrac",
                                                      frac.frac_a,"ifrac")
                cpl.map.mapPrint = cpl.debug.newDebugPrint(False)
                av_so2c_a = self.map_So2a.mapAV( self.av_so2c_o,
                                                 frac.frac_o,"afrac",
                                                 frac.frac_a,"ofrac")
                cpl.map.mapPrint = cpl.debug.newDebugPrint(False)
                # KLUDGE - START
                fix_frac_a = self.map_So2a.mapAV( frac.frac_o )
                fix_so2c_a = self.map_So2a.mapAV( self.av_so2c_o,
                                                  frac.frac_o,"afrac",
                                                  fix_frac_a,"afrac" )
                # KLUDGE - END
                self.av_fi2c_a = self.map_Fo2a.mapAV(self.av_fi2c_i,
                                                     frac.frac_i,"ifrac",
                                                     frac.frac_a,"ifrac")
                self.av_fo2c_a = self.map_Fo2a.mapAV(self.av_fo2c_o,
                                                     frac.frac_o,"afrac",
                                                     frac.frac_a,"ofrac")
                self.av_aoflux_a = self.map_Fo2a.mapAV(self.av_aoflux_o,
                                                       frac.frac_o,"afrac",
                                                       frac.frac_a,"ofrac")
                self.av_oalbedo_a = self.map_Fo2a.mapAV(self.av_oalbedo_o,
                                                        frac.frac_o,"afrac",
                                                        frac.frac_a,"ofrac")
                # Map ocn states/fluxes to atm -- END
                    
                # Recv message from lnd
                """
                call cpl_interface_contractRecv(cpl_fields_lndname,con_Xl2c)
                call cpl_bundle_mult(con_Xl2c%bundle,bun_areafact_l,"comp2cpl",
                                     bunlist=cpl_fields_l2c_fluxes)
                call cpl_bundle_split(con_Xl2c%bundle,bun_Sl2c_l,bun_Fl2c_l )
   
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_l) == 0 ):
                    runPrint("Receiving LND(l2c) contract...")
                    print "l2c size:",self.contracts["l2c"].av.lsize()
                    print "c2l size:",self.contracts["c2l"].av.lsize()
                    self.contracts["l2c"].recv()
                    # mult by comp2cpl
                    self.contracts["l2c"].av.mult( areafact.av_areafact_l, "comp2cpl",
                                                   cpl.fields.l2c_fluxes )
                    # split:
                    self.av_sl2c_l.copy( self.contracts["l2c"].av )
                    BADLAND = False
                    if( not BADLAND ):
                        self.av_fl2c_l.copy( self.contracts["l2c"].av )
                    else:
                        print "*** NOT COPYING POTENTIALLY BAD FLUXES FROM THE LAND MODEL! ***"
                    if(self.doLogging):
                        if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                            ni = self.contracts["l2c"].domain.ni
                            nj = self.contracts["l2c"].domain.nj
                            self.av_fl2c_l.writeNC(os.path.join(self.logPath,"fl2c-early-integration.nc"),ni,nj,ni*nj)
                    
                    
                    
                # recv from runoff next:
                # possible TYPO:  n is used here instead of n-1
                """
                call cpl_interface_contractRecv(cpl_fields_lndname,con_Xr2c)
                call cpl_bundle_mult(con_Xr2c%bundle,bun_areafact_r,"comp2cpl",
                                 bunlist=cpl_fields_r2c_fluxes)
                """
                if ( (n)%(cpl.control.ncpl_a/cpl.control.ncpl_r) == 0 ):
                    self.contracts["r2c"].recv()
                    # mult by comp2cpl
                    self.contracts["r2c"].av.mult( areafact.av_areafact_r, "comp2cpl",
                                                   cpl.fields.r2c_fluxes )
                
                # diagnostics: verify net solar calcs are coordinated
                # Merge atm states and fluxes

                statusPrint( "Merging ATM states and fluxes:")
                # CALL MERGE_ATM(fix_so2c_a) #+KLUDGE (fix_so2c_a is a kludge)
                #print "Calling Merge OCN inplace of Merge ATM!"
                #self.merge_ocn()
                self.merge_atm( self.fix_so2c_a )
                
                # Send message to atm
                """
                call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,"cpl2comp",
                                 bunlist=cpl_fields_c2a_fluxes)
                call cpl_interface_contractSend(cpl_fields_atmname,con_Xc2a)
                call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,"comp2cpl",
                                 bunlist=cpl_fields_c2a_fluxes)
                """
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_a) == 0 ):
                    runPrint("Sending ATM(c2a) contract...")
                    # Mult, send, mult
                    self.contracts["c2a"].av.mult( areafact.av_areafact_a, "cpl2comp",
                                                   cpl.fields.c2a_fluxes )
                    # Send
                    if(self.doLogging):
                        if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
                            ni = self.contracts["c2a"].domain.ni
                            nj = self.contracts["c2a"].domain.nj
                            self.av_fl2c_l.writeNC(os.path.join(self.logPath,"fl2c-integration.nc"),ni,nj,ni*nj)
                            self.contracts["c2a"].av.writeNC(os.path.join(self.logPath,"atm-integration.nc"),ni,nj,ni*nj)
                    self.contracts["c2a"].send()
                    # Mult by comp2cpl
                    self.contracts["c2a"].av.mult( areafact.av_areafact_a, "comp2cpl",
                                                   cpl.fields.c2a_fluxes )
                        
                # Map land to ocean, NOT allowed now
                """
                !        call shr_timer_start(tm1)
                !        call cpl_map_bun(bun_Sl2c_l,bun_Sl2c_o,map_Sa2o)
                !        call shr_timer_stop(tm1) ; call shr_timer_start(tm2)
                !        call cpl_map_bun(bun_Fl2c_l,bun_Fl2c_o,map_Fa2o)
                !        call shr_timer_stop(tm2)
                """

                """call cpl_map_bun(con_Xr2c%bundle,bun_Xr2c_o,map_Xr2o,mvector=r2ovector)"""
                if( self.doRiver ):
                    self.av_xr2c_o =  self.map_Xr2o.mapAV( self.contracts["r2c"].av )
                # Create history files
                # recv msg from OCN
                if ( (n)%(cpl.control.ncpl_a/cpl.control.ncpl_o) == 0 ):
                    runPrint("Receiving OCN(o2c) contract...")
                    runPrint("lagocn:",cpl.control.lagocn)
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
                        temp,tmpsize = self.av_xc2opsum_o.exportRAttr(r)
                        temp /= ncpl
                        self.av_xc2opsum_o.importRAttr(r,temp)
                    
                    #call cpl_bundle_copy(bun_Xc2oPSUM_o,outbun=con_Xc2o%bundle)
                    avsizePrint("self.contracts[\"c2o\"].av:",self.contracts["c2o"].av.size())
                    avsizePrint("self.av_xc2opsum_o:",self.av_xc2opsum_o.size())
                    self.contracts["c2o"].av.copy( self.av_xc2opsum_o )

                    if (not cpl.control.lagocn):
                        cpl.contract.debugPrint = cpl.debug.newDebugPrint(True)
                        cpl.infobuffer.debugPrint = cpl.debug.newDebugPrint(True)
                        self.contracts["o2c"].recv()
                        cpl.contract.debugPrint = cpl.debug.newDebugPrint(False)
                        cpl.infobuffer.debugPrint = cpl.debug.newDebugPrint(False)
                        # call cpl_bundle_mult( con_Xo2c%bundle,
                        #                       bun_areafact_o, "comp2cpl",
                        #                       bunlist=cpl_fields_o2c_fluxes)
                        #factor,tmpsize = areafact.av_areafact_o.exportRAttr("comp2cpl")
                        #for field in cpl.fields.o2c_fluxes:
                        #    temp,tmpsize = self.contracts["o2c"].av.exportRAttr(field)
                        #    temp *= factor
                        #    self.contracts["o2c"].av.importRAttr(field,temp)
                        self.contracts["o2c"].av.mult( areafact.av_areafact_o,
                                                       "comp2cpl",
                                                       cpl.fields.o2c_fluxes)
                            
                        self.av_so2c_o.copy(self.contracts["o2c"].av)
                        self.av_fo2c_o.copy(self.contracts["o2c"].av)
                    # --- Start normal interaction with ocn
                    if( cpl.control.lagocn ):
                        statusPrint( "\nStart of time coordinated integration\n")
                        cpl.control.lagocn = False
                    runPrint("Received OCN(o2c) contract!")

                sys.stdout.flush()
                sys.stderr.flush()

                # Recv message from atm
                if ( (n-1)%(cpl.control.ncpl_a/cpl.control.ncpl_a) == 0 ):
                    runPrint("Receiving ATM(a2c) contract ...")
                    sys.stdout.flush()
                    sys.stderr.flush()
                    self.contracts["a2c"].recv()
                    self.contracts["a2c"].av.mult( areafact.av_areafact_a, "comp2cpl",
                                                   cpl.fields.a2c_fluxes )
                    self.av_sa2c_a.copy( self.contracts["a2c"].av )
                    self.av_fa2c_a.copy( self.contracts["a2c"].av )
                    cpl.control.fluxashift = self.contracts["a2c"].infobuffer.ibufGet( "ashift" )
                    #---form total rain and snow
                    """
                    call cpl_bundle_zero (bun_precip_a)
                    call cpl_bundle_add(bun_precip_a,"Faxc_rain",bun_Fa2c_a,"Faxa_rainc")
                    call cpl_bundle_add(bun_precip_a,"Faxc_rain",bun_Fa2c_a,"Faxa_rainl")
                    call cpl_bundle_add(bun_precip_a,"Faxc_snow",bun_Fa2c_a,"Faxa_snowc")
                    call cpl_bundle_add(bun_precip_a,"Faxc_snow",bun_Fa2c_a,"Faxa_snowl")
                    """
                    faxa_rainc,size = self.av_fa2c_a.exportRAttr("Faxa_rainc")
                    faxa_rainl,size = self.av_fa2c_a.exportRAttr("Faxa_rainl")
                    faxc_rain = faxa_rainc + faxa_rainl
                    self.av_precip_a.importRAttr("Faxc_rain",faxc_rain)
                    faxa_snowc,size = self.av_fa2c_a.exportRAttr("Faxa_snowc")
                    faxa_snowl,size = self.av_fa2c_a.exportRAttr("Faxa_snowl")
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
                statusPrint("*** Warning:  date not updated! ***")
                
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
                statusPrint( "cpl: Finished hour %s of %s!"%(n+(cpl.control.ncpl_a*(iterations-1)), (cpl.control.ncpl_a*days)) )
                #if (stopeod!=0):
                #    statusPrint( "CPL STOPPING..."
                #    done = True
                #    break
                TIMER.stop("step")
            # End of NCPL Loop
            if (done):
                break
        statusPrint( self.name,": end of main integration loop...")
        # Send last message with a stopnow signal:
        for c in ["c2a","c2i","c2l","c2o"]:
            self.contracts[c].infobuffer.ibufSet("stopnow",1)
            self.contracts[c].send()
        TIMER.stop("run")
        return
    # End of run method

    def merge_atm( self, fix_so2c_a):
        """
        IROUTINE: merge_atm -- merge bundles to form atm input bundle
        DESCRIPTION:
            merge bundles to form atm input bundle
        """
        ni = self.contracts["c2a"].domain.ni
        nj = self.contracts["c2a"].domain.nj
        n = ni * nj
        # --- merge atm states & fluxes ---
        # call cpl_bundle_zero(con_Xc2a%bundle) ! zero before adding
        self.contracts["c2a"].av.zero()
        
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
        ofrac,    n_ofrac     = frac.frac_a.exportRAttr("ofrac")
        lfrac,    n_lfrac     = frac.frac_a.exportRAttr("lfrac")
        ifrac,    n_ifrac     = frac.frac_a.exportRAttr("ifrac")

        Target = "Faxx"
        FieldList = ["taux","tauy","lat","sen","lwup","evap"]
        Sources = ["Faoc","Fall","Faii"]
        #if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
        #    self.av_aoflux_a.writeNC(os.path.join(self.logPath,"merge-aoflux-a.nc"),ni,nj,n)
        #    self.av_fl2c_l.writeNC(os.path.join(self.logPath,"merge-fl2c-l.nc"),ni,nj,n)
        #    self.av_fi2c_a.writeNC(os.path.join(self.logPath,"merge-fi2c-a.nc"),ni,nj,n)
        #self.av_aoflux_a.extremes("av_aoflux_a")
        #self.av_fl2c_l.extremes("av_fl2c_l")
        #self.av_fi2c_a.extremes("av_fi2c_a")
        
        for field in FieldList:
            # Target:
            faxx,size_faxx = self.contracts["c2a"].av.exportRAttr(Target+"_"+field)
            # Sources in av_aoflux_a
            faoc,size_faoc = self.av_aoflux_a.exportRAttr("Faoc"+"_"+field)
            fall,size_fall = self.av_fl2c_l.exportRAttr("Fall"+"_"+field)
            faii,size_faii = self.av_fi2c_a.exportRAttr("Faii"+"_"+field)
            # Sum:
            # Bundle add works like this:
            # First argument is target (A in the equation A = A + B)
            # 2nd Argument is added to the first
            # 3rd->Nth arguments are multiplicative factors applied to the 2nd argument.
            # 1st = 1st + (2nd * 3rd * ... Nth)
            # So this one line is really 3 "cpl_bundle_add" methods together.
            faxx = (faoc * ofrac) + (fall * lfrac) + (faii * ifrac)
            self.contracts["c2a"].av.importRAttr(Target+"_"+field, faxx)

        # Tref:
        sx,size_sx = self.contracts["c2a"].av.exportRAttr("Sx_tref")
        faoc,size_faoc = self.av_aoflux_a.exportRAttr("Faoc_tref")
        sl,size_sl = self.av_sl2c_l.exportRAttr("Sl_tref")
        si,size_si = self.av_si2c_a.exportRAttr("Si_tref")
        sx = (faoc * ofrac) + (sl * lfrac) + (si * ifrac)
        self.contracts["c2a"].av.importRAttr("Sx_tref",sx)
        # Qref:
        sx,size_sx = self.contracts["c2a"].av.exportRAttr("Sx_qref")
        faoc,size_faoc = self.av_aoflux_a.exportRAttr("Faoc_qref")
        sl,size_sl = self.av_sl2c_l.exportRAttr("Sl_qref")
        si,size_si = self.av_si2c_a.exportRAttr("Si_qref")
        sx = (faoc * ofrac) + (sl * lfrac) + (si * ifrac)
        self.contracts["c2a"].av.importRAttr("Sx_qref",sx)


        # Just Sx fields:
        Target = "Sx"
        FieldList = ["avsdr","anidr","avsdf","anidf"]
        Sources = ["So","Sl","Si"]

        for field in FieldList:
            #target:
            sx,size_sx = self.contracts["c2a"].av.exportRAttr(Target+"_"+field)
            # sources:
            so,size_so = self.av_oalbedo_a.exportRAttr("So"+"_"+field)
            sl,size_sl = self.av_sl2c_l.exportRAttr("Sl"+"_"+field)
            si,size_si = self.av_si2c_a.exportRAttr("Si"+"_"+field)

            #Sum:
            sx = (so * ofrac) + (sl * lfrac) + (si * ifrac)
            self.contracts["c2a"].av.importRAttr(Target+"_"+field, sx)
            
        # Temperature:
        """Debugging temp:  dumping inputs to netcdf"""
        #if( cpl.comm.component_pid == cpl.comm.world_pe0_cpl ):
        #    frac.frac_a.writeNC(os.path.join(self.logPath,"merge-frac.nc"),ni,nj,n)
        #    self.av_so2c_a.writeNC(os.path.join(self.logPath,"merge-t-ocn.nc"),ni,nj,n)
        #    self.av_sl2c_l.writeNC(os.path.join(self.logPath,"merge-t-lnd.nc"),ni,nj,n)
        #    self.av_si2c_a.writeNC(os.path.join(self.logPath,"merge-t-ice.nc"),ni,nj,n)
        
        sx_t,size_sx = self.contracts["c2a"].av.exportRAttr("Sx_t")
        #self.av_so2c_a.extremes("av_so2c_a")
        so,size_so = self.av_so2c_a.exportRAttr("So_t")
        #self.av_sl2c_l.extremes("av_sl2c_l")
        sl,size_sl = self.av_sl2c_l.exportRAttr("Sl_t")
        #self.av_si2c_a.extremes("av_si2c_a")
        si,size_si = self.av_si2c_a.exportRAttr("Si_t")

        so_part = (so * ofrac)
        sl_part = (sl * lfrac)
        si_part = (si * ifrac)
                    
        sx_t = so_part + sl_part + si_part
        self.contracts["c2a"].av.importRAttr("Sx_t",sx_t)

        # Snow Height:
        snowh,size_snowh = self.contracts["c2a"].av.exportRAttr("Sx_snowh")
        #print self.av_sl2c_l
        # Here "Sl_snowh" is supposed to be used, but since only "Sl_snow"
        # was in the bundle I"m using it instead (MS)
        sl,size_sl = self.av_sl2c_l.exportRAttr("Sl_snow")
        snowh = (sl * lfrac)
        self.contracts["c2a"].av.importRAttr("Sx_snowh",snowh)

        # Kludge:  So_t values have errors due to mapping,
        # can be less then freezing!
        # Using an alternate calculation of So_t from a bundle
        # passed in as an arg:
        debugPrint( fix_so2c_a)
        so_t,size_so_t = fix_so2c_a.exportRAttr("So_t")
        self.contracts["c2a"].av.importRAttr("So_t",so_t)

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
        mergeocnPrint("merge_ocn: copying fields: ...")
        self.av_xc2osnap_o.copy( self.av_so2c_o )
        self.av_xc2osnap_o.copy( self.av_fo2c_o )
        self.av_xc2osnap_o.copy( self.av_si2c_i )
        self.av_xc2osnap_o.copy( self.av_fi2c_i )
        self.av_xc2osnap_o.copy( self.av_xr2c_o )
        self.av_xc2osnap_o.copy( self.av_aoflux_o )
        mergeocnPrint("merge_ocn: copying fields: finished!")

        npts = self.av_xc2osnap_o.lsize()

        foxx_taux,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_taux")
        foxx_tauy,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_tauy")
        foxx_swnet,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_swnet")
        foxx_lat,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_lat")
        foxx_sen,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_sen")
        foxx_lwup,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_lwup")
        foxx_evap,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_evap")
        foxx_lwdn,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_lwdn")
        foxx_rain,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_rain")
        foxx_snow,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_snow")
        foxx_prec,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_prec")
        foxx_melth,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_melth")
        foxx_meltw,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_meltw")
        foxx_salt,tmpsize = self.av_xc2osnap_o.exportRAttr("Foxx_salt")
        si_ifrac,tmpsize = self.av_xc2osnap_o.exportRAttr("Si_ifrac")

        faoc_taux,tmpsize = self.av_aoflux_o.exportRAttr("Faoc_taux")
        faoc_tauy,tmpsize = self.av_aoflux_o.exportRAttr("Faoc_tauy")
        faoc_swnet,tmpsize = self.av_aoflux_o.exportRAttr("Faoc_swnet")
        faoc_lat,tmpsize = self.av_aoflux_o.exportRAttr("Faoc_lat")
        faoc_sen,tmpsize = self.av_aoflux_o.exportRAttr("Faoc_sen")
        faoc_lwup,tmpsize = self.av_aoflux_o.exportRAttr("Faoc_lwup")
        faoc_evap,tmpsize = self.av_aoflux_o.exportRAttr("Faoc_evap")

        faxc_rain,tmpsize = self.av_precip_o.exportRAttr("Faxc_rain")
        faxc_snow,tmpsize = self.av_precip_o.exportRAttr("Faxc_snow")
        
        fioi_taux,tmpsize = self.av_fi2c_i.exportRAttr("Fioi_taux")
        fioi_tauy,tmpsize = self.av_fi2c_i.exportRAttr("Fioi_tauy")
        fioi_swpen,tmpsize = self.av_fi2c_i.exportRAttr("Fioi_swpen")
        faxa_lwdn,tmpsize = self.av_fa2c_o.exportRAttr("Faxa_lwdn")
        fioi_melth,tmpsize = self.av_fi2c_i.exportRAttr("Fioi_melth")
        fioi_meltw,tmpsize = self.av_fi2c_i.exportRAttr("Fioi_meltw")
        fioi_salt,tmpsize = self.av_fi2c_i.exportRAttr("Fioi_salt")

        afrac,frac_o_size = frac.frac_o.exportRAttr("afrac")
        ifrac,frac_o_size = frac.frac_o.exportRAttr("ifrac")
        

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
            statusPrint( "len(): foxx_rain -- faxc_rain -- afrac -- npts")
            statusPrint( "len():",len( foxx_rain ), len(faxc_rain), len(afrac),npts)
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

def __del__(self):
    cpl.comm.finalize()

def main():
    statusPrint( "pyCPL: Start!")
    global TIMESTEP
    # Coupler Path, doComputations, doMapping, initializeRiverMapping,Logging)
    TIMER.start("total")
    TIMER.start("init")
    statusPrint("pyCPL: Coupler root is ",ROOT)
    c = Coupler("cpl.nml",True,True,True,True,"~/exe/pyccsm/all/")
    TIMER.stop("init")
    c.run(5)
    TIMER.stop("total")
    #print "StartDate:",cpl.date
    statusPrint( "pyCPL: Done!" )
    print TIMER
    
if __name__=="__main__":
    main()

    

