# Builtin Python Modules
import math

# External Modules
import mpi
import Numeric as pyarray
import MCT
from MCT import SparseMatrix
# REQUIRES NetCDF
import pycdf #provides interface to netcdf

# Internal Modules
import bundle
import attributevector
import comm
import control
import const
import domain
import fields
import rearranger

import debug

debugPrint = debug.newDebugPrint(False)
initPrint = debug.newDebugPrint(True)
mapPrint = debug.newDebugPrint(False)
testPrint = debug.newDebugPrint(True)

DEBUG      = control.infodbug
BFBFLAG    = control.bfbflag 
AREA_FIELD = "aream"

NPFIX_FIRSTCALL=True
global_gathered_grid = attributevector.AttributeVector()

def nint( num ):
    return int(round(num))

class Map:
    """
    MCT Currently does not support the creation of a map between two arbitrary 
    domains, rather one of the two domains must have a decomposition
    of MCT's choosing, thus the calling routine can specify which of the two
    domains should be altered and this new/required domain can be returned.
    """
    def __init__(self, dom_src=None, dom_dst=None, mapName=None,
                 fileName=None, newdom=None, adj_areas=None ):
        self.sMat = SparseMatrix.SparseMatrix()
        self.sMatGlobal = SparseMatrix.SparseMatrix()
        self.new = domain.Domain()
        # Rearranger:
        self.rearranger = rearranger.Rearranger()
        self.lrearranger = rearranger.Rearranger()
        # Area of source grid from mapping file
        self.areasrc = attributevector.AttributeVector([],[AREA_FIELD],10)
        self.areasrcGlobal = attributevector.AttributeVector([],[AREA_FIELD],10)
        # Area of destination grid from mapping file
        self.areadst = attributevector.AttributeVector([],[AREA_FIELD],10)
        self.areadstGlobal = attributevector.AttributeVector([],[AREA_FIELD],10)
        self.tag = const.default_tag
        self.sum = False
        self.firstrun = True
        return

    def initialize(self, dom_src=None, dom_dst=None, mapName=None,
                   fileName=None, newdom=None, adj_areas=None):
        ## Set map's src & dst domains, name, "new mct domain" choice
        self.mapname = mapName
        self.filename = fileName
        # Associated source and destination domains:
        self.source = dom_src
        self.destination = dom_dst
        self.newtype = newdom # Can be "SRC" or "DST"
        self.idtype = 0 # 0 = Normal, 1 = Identity(ID)
        if( control.bfbflag ):
            self.newtype = "src"
            initPrint( "bfbflag = "+str(control.bfbflag) )
            initPrint( ": overwriting newtype from "+str(newdom)+'to src' )
        
        ## read & test the map data on root processor(only)
        if( comm.component_pid == 0 ):
            self.read( self.filename )
            if( True ):#control.infodbug > 1 ):
                initPrint( "running global map test, dbug level = "+str(control.infodbug))
                self.testGlobal()
            else:
                initPrint( "skipping map test, dbug level = "+str(control.infodbug))
        # Force all coupler processors to wait for the read
        mpi.barrier( comm.local_comm )
        initPrint("Allocate local AttrVects...")
        src_size = mpi.bcast( self.areasrc.lsize(), 1, mpi.MPI_INT,
                              0, comm.local_comm )
        dst_size = mpi.bcast( self.areadst.lsize(), 1, mpi.MPI_INT,
                              0, comm.local_comm )
        self.areasrc = attributevector.AttributeVector([],[AREA_FIELD],src_size[0])
        self.areadst = attributevector.AttributeVector([],[AREA_FIELD],dst_size[0])
        
        ## scatter map_X data and create new/intermediate domain as required by MCT
        initPrint( "scattering map_X data" )
        initPrint( "\tareasrc.scatter:..." )

        initPrint( "\tsource.gsMap.gsize: "+str(self.source.gsMap.globalsize()))
        if( comm.component_pid == 0 ):
            initPrint( "\tareaSRCGlobal.lsize: "+str(self.areasrcGlobal.size()))
        initPrint( "\tcomm.component_pid: "+str(comm.component_pid))
        initPrint( "\tcomm.local_comm: "+str(comm.local_comm))
        self.areasrc.scatter( self.areasrcGlobal, self.source.gsMap, 0,
                                    comm.component_pid, comm.local_comm )
        initPrint( "\tareasrc.scatter:OK!" )
        initPrint( "\tareadst.scatter:..." )
        if( comm.component_pid == 0 ):
            initPrint( "\tareaDSTGlobal.lsize: "+str(self.areadstGlobal.size()))
        self.areadst.scatter( self.areadstGlobal, self.destination.gsMap, 0,
                                    comm.component_pid, comm.local_comm )
        initPrint( "\tareadst.scatter:OK!" )
        initPrint( "scattered map_X data" )
        
        if (control.infodbug > 2 ):
            pass
        
        ## create new/intermediate destination domain: create it:
        if( self.newtype == "dst" ):
            self.initializeScatterByColumns(dom_src, dom_dst, mapName,
                                         fileName, newdom, adj_areas)
        elif( self.newtype =="src" ):
            self.initializeScatterByRows(dom_src, dom_dst, mapName,
                                    fileName, newdom, adj_areas)
        else:
            initPrint( "Invalid newdom value = "+self.newtype)

        ## Adjust mapping weights based on areas of domains
        ## vs areas in mappings:
        if (adj_areas!=None):
            initPrint( "WARNING: Do not use adj_areas option right now!")

        initPrint( "done initializing map: "+self.mapname )
        return

    
    def initializeScatterByColumns(self, dom_src=None, dom_dst=None,
                                   mapName=None, fileName=None,
                                   newdom=None, adj_areas=None):
        initPrint( "Scatter Matrix by Column...")
        if( control.infodbug == 1 ):
            if( comm.component_pid == 0 ):
                initPrint( "Global sMat rows x cols = "+str(self.sMatGlobal.nRows())+'x'+str(self.sMatGlobal.nCols()))
                rows,length = self.sMatGlobal.exportGlobalRowIndices()
                initPrint( "length of global row indices array:"+str(length))
                # initPrint( "rows:",rows

        self.sMat.ScatterByColumn( self.source.gsMap,
                                   self.sMatGlobal,
                                   0,
                                   comm.component_pid,
                                   comm.local_comm )

        if( comm.component_pid == 0 ):
            """
            Executes a number of routines on the SparseMatrix
            that hang on a non-root processor in this component.
            """
            self.testLocal()
        
        if( control.infodbug == 1 ):
            if( comm.component_pid == 0 ):
                initPrint("Global sMat rows x cols = "+str(self.sMatGlobal.nRows())+'x'+str(self.sMatGlobal.nCols()))
                initPrint("scattered sMat rows x cols = "+str(self.sMat.nRows())+'x'+str(self.sMat.nCols()))
                initPrint("scattered sMat lSize "+str(self.sMat.lsize()))
                
            gnumels = self.sMat.GlobalNumElements(comm.local_comm)
            if(comm.component_pid == 0 ):
                initPrint( "scattered sMat GnumEl "+str(gnumels) )
                rows,length = self.sMat.exportGlobalRowIndices()
                initPrint( "length of global row indices array:"+str(length) )
                # initPrint( "rows:",rows

        # create new dst intermediate map
        initPrint( "Creating new destinaton intermediate map..." )
        if (comm.component_pid):
            initPrint( "self.destination.lsize(): "+str(self.destination.gsMap.lsize(comm.local_comm)))
            initPrint( "self.destination.gsMap.ProcessStorage(0): "+str(self.destination.gsMap.ProcessStorage(0)))
            initPrint( "self.source.gsMap.ProcessStorage(0): "+str(self.source.gsMap.ProcessStorage(0)))
        # Copy destination to new:
        self.new.name = "[%s]->[%s]"%(self.source.name,
                                  self.destination.name)
        self.new.suffix = self.destination.suffix
        self.new.n  = self.destination.n
        self.new.ni = self.destination.ni
        self.new.nj = self.destination.nj
        # This modified self.new.gsMap in place
        initPrint( "Before smat.sparsematrixtoyglobalsegmap:" )
        initPrint( "comm.mph_component_id = "+str(comm.mph_component_id))
        self.sMat.SparseMatrixToYGlobalSegMap(self.new.gsMap,
                                              0,
                                              comm.local_comm,
                                              comm.mph_component_id
                                              )
        initPrint( "self.new.gsMap" )
        #initPrint( "\ttype of:",type(self.new.gsMap)
        #initPrint( "\tdir of :",dir(self.new.gsMap)
        initPrint( "\tgsize of new: "+str(self.new.gsMap.globalsize()))
        initPrint( "\tlocal size of new: "+str(self.new.gsMap.lsize(comm.local_comm)))
        initPrint( "self.destination.gsMap" )
        #initPrint( "\ttype of:",type(self.destination.gsMap)
        #initPrint( "\tdir of:",dir(self.destination.gsMap)
        initPrint( "\tgsize of destination: "+str(self.destination.gsMap.globalsize()))
        initPrint( "\tlocal size of destination: "+str(self.destination.gsMap.lsize(comm.local_comm)))

        if (control.infodbug > 1):
            initPrint( "new gsMap gSize = "+str(self.new.gsMap.globalsize()))
            initPrint( "new gsMap lSize = "+str(self.new.gsMap.lsize(comm.local_comm)))

        initPrint( "Created new destination intermediate map..." )
        ##
        # these modify self.sMat in place.
        initPrint( "GlobalToLocalMatrix: ROW" )
        self.new.gsMap.GlobalToLocalMatrix( self.sMat, "ROW", comm.local_comm )
        initPrint( "self.new.gsMap.ProcessStorage(0): "+str(self.new.gsMap.ProcessStorage(0)))
        initPrint( "GlobalToLocalMatrix: COLUMN")
        self.source.gsMap.GlobalToLocalMatrix( self.sMat, "COLUMN", comm.local_comm )
        initPrint( "self.new.gsMap.ProcessStorage(0): "+str(self.new.gsMap.ProcessStorage(0)))

        initPrint( 'Rearranger Init')
        initPrint( '\tdestination gsmap compid: '+str(self.destination.gsMap.compid()))
        initPrint( '\tnew gsmap compid: '+str(self.new.gsMap.compid()) )
        initPrint( "\tdestination lgrid lsize: "+str(self.destination.lgrid.size()))
        self.rearranger.init( self.new.gsMap, self.destination.gsMap, comm.local_comm )
        initPrint( 'Local Rearranger Init')
        initPrint( "\tdestination lgrid lsize: "+str(self.destination.lgrid.size()))
        self.lrearranger.init(self.destination.gsMap, self.new.gsMap, comm.local_comm )
        initPrint( 'Setup New\'s Lgrid')
        initPrint( "\tdestination lgrid lsize: "+str(self.destination.lgrid.size()))
        
        self.new.lgrid.initv(self.destination.lgrid, self.new.gsMap.lsize(comm.local_comm) )
        ifields,rfields = self.new.lgrid.getFields()
        size = self.new.lgrid.lsize()
        initPrint( "\t*** new lgrid ***\n\t",ifields,"\n\t",rfields,"\n\t","size=",size)
        dbg,size = self.destination.lgrid.exportRAttr("aream")
        initPrint( "\tslice of destination.lgrid's aream field:",dbg[0:200])
        dbg,size = self.new.lgrid.exportRAttr("aream")
        initPrint( "\tslice of new.lgrid's aream field:",dbg[0:200] )
         
        initPrint( 'Local Rearranger Rearranging things...' )
        initPrint( "\tdestination lgrid lsize: "+str(self.destination.lgrid.size()))
       
        #initPrint( "\t**this is currently being skipped!**"
        # Expects real MCT AV type,
        self.lrearranger.rearrange( self.destination.lgrid, self.new.lgrid, self.tag, self.sum)
        initPrint( 'Local Rearranging done!' )
        initPrint( "\tdestination lgrid lsize: "+str(self.destination.lgrid.size()) )

        # call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%new%gsMap,"ROW"   ,cpl_comm_comp)
        # call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%src%gsMap,"COLUMN",cpl_comm_comp)
        # call cpl_mct_rearr_init (map_X%new%gsMap,map_X%dst%gsMap,cpl_comm_comp,map_X%rearr)
        # call cpl_mct_rearr_init (map_X%dst%gsMap,map_X%new%gsMap,cpl_comm_comp,lrearr)
        # call cpl_mct_aVect_init (map_X%new%lGrid,map_X%dst%lGrid,cpl_mct_gsMap_lSize(map_X%new%gsMap,cpl_comm_comp))
        # call cpl_mct_rearr_rearrange(map_X%dst%lGrid,map_X%new%lGrid,lrearr)
        # call cpl_mct_rearr_clean(lrearr)
        return # DST Init END
    
    def initializeScatterByRows(self, dom_src=None, dom_dst=None, mapName=None,
                   fileName=None, newdom=None, adj_areas=None):
        initPrint( "Scatter Matrix by Row...")
        if( control.infodbug == 1 ):
            if( comm.component_pid == 0 ):
                initPrint( "Global sMat rows x cols ="+str(self.sMatGlobal.nRows())+'x'+str(self.sMatGlobal.nCols()) )
                rows,length = self.sMatGlobal.exportGlobalRowIndices()
                initPrint( "length of global row indices array:"+str(length) )
                # initPrint( "rows:",rows
            
        self.sMat.ScatterByRow( self.source.gsMap,
                                   self.sMatGlobal,
                                   0,
                                   comm.component_pid,
                                   comm.local_comm )
        if( comm.component_pid == 0 ):
            """
            Executes a number of routines on the SparseMatrix
            that hang on a non-root processor in this component.
            """
            self.testLocal()

        if( control.infodbug == 1 ):
            if( comm.component_pid == 0 ):
                initPrint( "Global sMat rows x cols = "+str(self.sMatGlobal.nRows())+'x'+str(self.sMatGlobal.nCols() ))
                initPrint( "scattered sMat rows x cols = "+str(self.sMat.nRows())+'x'+str(self.sMat.nCols()))
                initPrint( "scattered sMat lSize "+str(self.sMat.lsize()))
                
            gnumel = self.sMat.GlobalNumElements(comm.local_comm)
            if( comm.component_pid == 0 ):
                initPrint( "scattered sMat GnumEl "+str( gnumel ))
                #rows,length = self.sMat.exportGlobalRowIndices()
                #initPrint( "length of global row indices array:",length
                # initPrint( "rows:",rows
                
        # create new dst intermediate map
        initPrint( "Creating new destinaton intermediate map..." )
        if( comm.component_pid == 0 ):
            initPrint( "self.destination.lsize(): "+str(self.destination.gsMap.lsize(comm.local_comm)))
            initPrint( "self.destination.gsMap.ProcessStorage(0): "+str(self.destination.gsMap.ProcessStorage(0)))
            initPrint( "self.source.gsMap.ProcessStorage(0): "+str(self.source.gsMap.ProcessStorage(0)))
        # Copy source to new:
        self.new.name = "[%s]->[%s]"%(self.source.name,
                                  self.destination.name)
        self.new.suffix = self.source.suffix
        self.new.n  = self.source.n
        self.new.ni = self.source.ni
        self.new.nj = self.source.nj
        # This modified self.new.gsMap in place
        initPrint( "Before smat.sparsematrixtoyglobalsegmap:" )
        initPrint( "comm.mph_component_id = "+str(comm.mph_component_id))
        self.sMat.SparseMatrixToXGlobalSegMap(self.new.gsMap,
                                                           0,
                                                           comm.local_comm,
                                                           comm.mph_component_id
                                                           )
        initPrint( "self.new.gsMap" )
        #initPrint( "\ttype of:",type(self.new.gsMap)
        #initPrint( "\tdir of :",dir(self.new.gsMap)
        initPrint( "\tgsize of new:"+str(self.new.gsMap.globalsize()))
        initPrint( "\tlocal size of new:"+str(self.new.gsMap.lsize(comm.local_comm)))
        initPrint( "self.destination.gsMap")
        #initPrint( "\ttype of:",type(self.destination.gsMap)
        #initPrint( "\tdir of:",dir(self.destination.gsMap)
        initPrint( "\tgsize of destination:"+str(self.destination.gsMap.globalsize()))
        initPrint( "\tlocal size of destination:"+str(self.destination.gsMap.lsize(comm.local_comm)))
        
        if (control.infodbug > 1):
            initPrint( "new gsMap gSize ="+str(self.new.gsMap.globalsize()))
            initPrint( "new gsMap lSize ="+str(self.new.gsMap.lsize(comm.local_comm)))
            
        initPrint( "Created new destination intermediate map..." )
        ##
        # these modify self.sMat in place.
        initPrint( "GlobalToLocalMatrix: ROW")
        self.destination.gsMap.GlobalToLocalMatrix( self.sMat, "ROW", comm.local_comm )
        initPrint( "self.new.gsMap.ProcessStorage(0):"+str(self.new.gsMap.ProcessStorage(0)))
        initPrint( "GlobalToLocalMatrix: COLUMN")
        self.new.gsMap.GlobalToLocalMatrix( self.sMat, "COLUMN", comm.local_comm )
        initPrint( "self.new.gsMap.ProcessStorage(0):"+str(self.new.gsMap.ProcessStorage(0)))
        
        initPrint( 'Rearranger Init')
        initPrint( '\tdestination gsmap compid:'+str(self.destination.gsMap.compid()))
        initPrint( '\tnew gsmap compid:'+str(self.new.gsMap.compid()))
        self.rearranger.init( self.source.gsMap, self.new.gsMap, comm.local_comm )
        initPrint( 'Setup New\'s Lgrid')
        self.new.lgrid.initv(self.destination.lgrid, self.new.gsMap.lsize(comm.local_comm) )
        initPrint( 'Rearranger Rearranging things...')
        self.rearranger.rearrange( self.source.lgrid, self.new.lgrid, self.tag, self.sum)
        initPrint( 'Rearranging done!')
        """
        call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%dst%gsMap,"ROW"   ,cpl_comm_comp)
        call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%new%gsMap,"COLUMN",cpl_comm_comp)
        call cpl_mct_rearr_init (map_X%src%gsMap,map_X%new%gsMap,cpl_comm_comp,map_X%rearr)
        call cpl_mct_aVect_init (map_X%new%lGrid,map_X%dst%lGrid,cpl_mct_gsMap_lSize(map_X%new%gsMap,cpl_comm_comp))
        call cpl_mct_rearr_rearrange(map_X%src%lGrid,map_X%new%lGrid,map_X%rearr)
        """
        return # SRC Init END
        
    def __del__( self ):
        """
       Deallocates local mct types.
        """
        #debugPrint( "Warning:  this destructor is dangerous!"
        # Domain and Map objects:
        
        # These python objects don't need to be clean()ed
        #self.smat.clean()
        #self.rearr.clean()
        #self.new.clean() 
        return

    def testLocal(self):
        testPrint( "consistancy check on LOCAL mapping matrix data..." )
        testPrint( " map name : ", self.mapname, 
                   " map src domain name : ", self.source.name,
                   " map dst domain name : ", self.destination.name,
                   " map IDtype : ", self.newtype )
        
        ns = self.sMat.lsize()
        igrow,igrow_size = self.sMat.exportGlobalRowIndices()
        igcol,igcol_size = self.sMat.exportGlobalColumnIndices()
        iwgt,iwgt_size = self.sMat.exportMatrixElements()

        testPrint("Rows x Columns = %s x %s"%(self.sMatGlobal.nRows(),self.sMatGlobal.nCols()))
        testPrint("* number of non-zero elements:",ns)
        
        minwgt,maxwgt = iwgt[0],iwgt[0]
        mini,maxi = 0,0
        for i,v in enumerate(iwgt):
            if( v <= minwgt ):
                minwgt = v
                mini = i
            if( v >= maxwgt ):
                maxwgt = v
                maxi = i

        testPrint( "* min weight value: iwgt[", mini , "] =", minwgt )
        testPrint( "* max weight value: iwgt[", maxi , "] =", maxwgt )

        minrow,maxrow = igrow[0],igrow[0]
        mini,maxi = 0,0
        for i,v in enumerate(igrow):
            if( v <= minrow ):
                minrow = v
                mini = i
            if( v >= maxrow ):
                maxrow = v
                maxi = i
        
        testPrint( "* first/last row values: %s / %s "%(igrow[0],igrow[ns-1]) )
        testPrint( "* min row value: igRow[%s] = %s"%(mini,minrow))
        testPrint( "* max row value: igRow[%s] = %s"%(maxi,maxrow))

        mincol,maxcol = igcol[0],igcol[0]
        mini,maxi = 0,0
        for i,v in enumerate(igcol):
            if( v <= mincol ):
                mincol = v
                mini = i
            if( v >= maxcol ):
                maxcol = v
                maxi = i
        
        testPrint( "* first/last column values: %s / %s"%(igcol[0],igcol[ns-1]))
        testPrint( "* min column value: igCol[%s] = %s"%(mini,mincol))
        testPrint( "* max column value: igCol[%s] = %s"%(maxi,maxcol))
        return
    # end of testLocal
    
    def testGlobal( self ):
        testPrint( "consistancy check on GLOBAL mapping matrix data..." )
        testPrint( " map name : ", self.mapname, 
                   " map src domain name : ", self.source.name,
                   " map dst domain name : ", self.destination.name,
                   " map IDtype : ", self.newtype )
        ns = self.sMatGlobal.lsize()
        igrow,igrow_size = self.sMatGlobal.exportGlobalRowIndices()
        igcol,igcol_size = self.sMatGlobal.exportGlobalColumnIndices()
        iwgt,iwgt_size = self.sMatGlobal.exportMatrixElements()

        testPrint("Rows x Columns = %s x %s"%(self.sMatGlobal.nRows(),self.sMatGlobal.nCols()))
        testPrint("* number of non-zero elements:",ns)
        
        minwgt,maxwgt = iwgt[0],iwgt[0]
        mini,maxi = 0,0
        for i,v in enumerate(iwgt):
            if( v <= minwgt ):
                minwgt = v
                mini = i
            if( v >= maxwgt ):
                maxwgt = v
                maxi = i

        testPrint( "* min weight value: iwgt[", mini , "] =", minwgt )
        testPrint( "* max weight value: iwgt[", maxi , "] =", maxwgt )

        minrow,maxrow = igrow[0],igrow[0]
        mini,maxi = 0,0
        for i,v in enumerate(igrow):
            if( v <= minrow ):
                minrow = v
                mini = i
            if( v >= maxrow ):
                maxrow = v
                maxi = i
        
        testPrint( "* first/last row values: %s / %s "%(igrow[0],igrow[ns-1]) )
        testPrint( "* min row value: igRow[%s] = %s"%(mini,minrow))
        testPrint( "* max row value: igRow[%s] = %s"%(maxi,maxrow))

        mincol,maxcol = igcol[0],igcol[0]
        mini,maxi = 0,0
        for i,v in enumerate(igcol):
            if( v <= mincol ):
                mincol = v
                mini = i
            if( v >= maxcol ):
                maxcol = v
                maxi = i
        
        testPrint( "* first/last column values: %s / %s"%(igcol[0],igcol[ns-1]))
        testPrint( "* min column value: igCol[%s] = %s"%(mini,mincol))
        testPrint( "* max column value: igCol[%s] = %s"%(maxi,maxcol))
        
        """ 
        !  NOTES: 
        !  o stuff below is removed because it alters the mapping data.
        !  o any testing done should alter the data being tested.
        !  o if any sorting is required, do it on file input
        !
        !  type(cpl_mct_list) :: sortKey       ! SparseMatrix sorting key list:
        !  logical        :: descending(2) ! Descending order flags for sort test 2a
        !
        !  !----------------------------------------------------------------------------
        !  ! sort and re-test
        !  !----------------------------------------------------------------------------
        !
        !  write(6,F00) "Descending sort test..."
        !
        !  call cpl_mct_list_init(sortKey,"grow:gcol")
        !  descending = .true.
        !  call cpl_mct_sMat_SortPermute(sMat, sortKey, descending)
        !
        !  write(6,*) subName, "sMat%data%iAttr(igrow, 1) = ",sMat%data%iAttr(igrow,1)
        !  write(6,*) subName, "sMat%data%iAttr(igcol, 1) = ",sMat%data%iAttr(igcol,1)
        !  write(6,*) subName, "sMat%data%iAttr(igrow,ns) = ",sMat%data%iAttr(igrow,ns)
        !  write(6,*) subName, "sMat%data%iAttr(igcol,ns) = ",sMat%data%iAttr(igcol,ns)
        !
        !  !----------------------------------------------------------------------------
        !  ! sort and re-test
        !  !----------------------------------------------------------------------------
        !
        !  write(6,F00) "Ascending sort test..."
        !
        !  call cpl_mct_sMat_SortPermute(sMat,sortKey)
        !
        !  write(6,*) subName, "sMat%data%iAttr(igrow, 1) = ",sMat%data%iAttr(igrow,1)
        !  write(6,*) subName, "sMat%data%iAttr(igcol, 1) = ",sMat%data%iAttr(igcol,1)
        !  write(6,*) subName, "sMat%data%iAttr(igrow,ns) = ",sMat%data%iAttr(igrow,ns)
        !  write(6,*) subName, "sMat%data%iAttr(igcol,ns) = ",sMat%data%iAttr(igcol,ns)
        """
        return
    # End testGlobal

    # Info (with help from __str__ and __repr__
    def info( self ):
        print "mapping name:",self.mapname
        print "mapping type:",self.newtype
        print "mapping ID type:",self.idtype
        print "mapping src domain:"
        self.source.info()
        # insert domain.info call
        print "mapping dst domain:"
        self.destination.info()
        # insert domain info call
        print "mapping new domain:"
        self.new.info()
        print "mapping SparseMatrix:",self.sMat
        print "SparseMatrix nrows:",self.sMat.nRows()
        print "SparseMatrix ncols:",self.sMat.nCols()
        # call cpl_mct_aVect_info(2,mapping%sMat%data,cpl_comm_comp,cpl_comm_comp_pid)
        return

    def __str__(self):
        return "PyMCT.Map: %s - %s"%(self.mapname, self.filename)

    def __repr__(self):
        return "PyMCT.Map: %s - %s"%(self.mapname, self.filename)

    def read(self,filename ):
        """
        Modifies the internal state of this map object to contain 
        all the appropriate data gathered from the netcdf mapping 
        file whose path is provided above as the filename argument.
        """
        # Example of reading and writing a simple field:
        
        debugPrint('reading file...')
        file = pycdf.CDF( filename )
        # Read dimensions of file:
        self.dimensionsDict = file.dimensions()
        self.sizeSparseMatrix = self.dimensionsDict['n_s']
        self.sizeInputVector  = self.dimensionsDict['n_a']
        self.sizeOutputVector = self.dimensionsDict['n_b']
        debugPrint('* matrix dimensions rows x cols :'+str(self.sizeInputVector)+'x'+str(self.sizeOutputVector))
        debugPrint('* number of non-zero elements:'+str(self.sizeSparseMatrix))

        
        # Read Variables of file:
        # Read row, col, s, areaInput('area_a'), and areOutput('area_b')
        # file.var("name") returns a pyCDF Var that points to "row"
        # Get returns the actual data, in this case correctly typed numeric arrays.
        self.variablesDict = file.variables()
        self.row = file.var("row").get()
        self.col = file.var("col").get()
        self.s = file.var("S").get()
        self.areaInput = file.var("area_a").get()
        self.areaOutput = file.var("area_b").get()
        file.close()
        
        # Init & load the mct sMat data type:
        self.sMatGlobal.init( self.sizeOutputVector, self.sizeInputVector, self.sizeSparseMatrix )
        # cpl_mct_sMat_init must be given the number of rows and
        # columns that would be in the full matrix.
        # nRows = size of output vector
        # nCols = size of input vector
        # index(I|R)A( 'attribute vector field', 'error string', 'die string' )
        
        """
        'exportGlobalColumnIndices', 'exportGlobalRowIndices',
        'exportLocalColumnIndices', 'exportLocalRowIndices',
        'exportMatrixElements', 'getClassInfo', 'global_col_range',
        'global_row_range', 'importGlobalColumnIndices', 'importGlobalRowIndices',
        'importLocalColumnIndices', 'importLocalRowIndices', 'importMatrixElements'
        """
        # Can we do this? We need to import export right?
        #self.sMatGlobal.iAttr(igrow,:) = row
        debugPrint("rows("+str(len(self.row))+"):"+str(self.row[:10])+"\n"+str(self.row[-10:]))
        self.sMatGlobal.importGlobalRowIndices( self.row, self.sizeSparseMatrix)
        junkrows,length = self.sMatGlobal.exportGlobalRowIndices( )
        debugPrint("junkrows("+str(length)+
                   ")"+str(junkrows[:10])+"\n"+str(junkrows[-10:]))
        #debugPrint(self.row
        #self.sMatGlobal.iAttr(igcol,:) = col
        self.sMatGlobal.importGlobalColumnIndices( self.col, self.sizeSparseMatrix)
        #debugPrint(self.col
        #self.sMatGlobal.iRttr(iwgt, :) = s
        self.sMatGlobal.importMatrixElements( self.s,self.sizeSparseMatrix )
        
        self.areasrcGlobal = attributevector.AttributeVector( [],
                                                              [AREA_FIELD],
                                                              self.sizeInputVector )        
        self.areadstGlobal = attributevector.AttributeVector( [],
                                                              [AREA_FIELD],
                                                              self.sizeOutputVector )
        self.areasrc = attributevector.AttributeVector([],[AREA_FIELD],self.sizeInputVector)
        self.areadst = attributevector.AttributeVector([],[AREA_FIELD],self.sizeOutputVector)
        #iarea = self.areasrcGlobal.indexRA(AREA_FIELD)
        #self.areasrcGlobal.rAttr[ iarea, 1:sizeInputVector ] = \
        #                        areaInput[1:sizeInputVector]
        self.areasrcGlobal.importRAttr( AREA_FIELD, self.areaInput )                        
        #iarea = self.areadstGlobal.indexRA(AREA_FIELD)
        #self.areadstGlobal.rAttr[ iarea, 1:sizeOutputVector ] = \
        #                        areaOutput[1:sizeOutputVector]
        self.areadstGlobal.importRAttr( AREA_FIELD, self.areaOutput )
        debugPrint('... done reading file')
        return
        
    def mapBundle( self, inbun, outbun ):
        """
        Map from one bundle to another.
        """
        return

    def mapAV( self, inav, avfs=None, fsname=None, avfd=None, fdname=None ):
        """
        Map a attribute vector of fields from one domain
        to another*.

        outav = mapAV( self, inav, avfs=None, fsname=None, avfd=None, fdname=None )

        ARGUMENTS:
         inav,        ! input bundle (attrvect)
        OPTIONAL:
         avfs=None,  ! src fraction input bundle
         fsname=None, ! Name of a field in avfs
         avfd=None,  ! destination fraction input bundle
         fdname=None  ! Name of a field in avfd

        RETURNS:
         outav        ! output bundle(attrvect)
         
        *It is assumed tht both vectors have an identical
        list of fields.
        """
        ifields,rfields = inav.getFields()
        outav = None # Placeholder
        ## identity map -> mapping is simply a copy of the bundle:
        if( self.idtype == 1):
            mapPrint("mapAV: Identity Mapping")
            outav.copy( inav )
            return
        
        normalize = False
        if ( avfs and avfd and fsname and fdname ):
            ## If Normalization is requested, create a local copy of the vector:
            inav_local = attributevector.AttributeVector( ifields, rfields, inav.lsize() )
            normalize = True
            #self.info()
            mapPrint("normalization vectors:")
            mapPrint("avfs(size):"+str(avfs))
            mapPrint("avfd(size):"+str(avfd))
            #ifields, rfields = inav.getFields()
            #inav_local.initialize( ifields, rfields, inav.size() )
            inav_local.initv( inav, inav.lsize() )
            #mapPrint("*** Warning:  MUST call bundle mult! ***"
            # call cpl_bundle_mult( inav_local, bunfs, fsname )
            # How we'll do bundle mult:
            frac,size = avfs.exportRAttr( fsname )
            for f in rfields:
                if( f == "index" ):
                    debugPrint("*** Not normalizing 'index' field! ***")
                    continue
                temp,size = inav.exportRAttr( f )
                temp *= frac
                inav_local.importRAttr( f, temp )
        elif( avfs == avfd == fsname == fdname == None ):
            pass
        else:
            debugPrint("mapAV: ERROR: optional arguments inconsistent")
            
        ##
        ## Do the mapping:
        ##
        if (self.newtype =="src" ):    
            av = attributevector.AttributeVector(ifields,rfields,self.new.lgrid.lsize())
            outav = attributevector.AttributeVector( ifields, rfields,
                                                     self.destination.lgrid.lsize() )
        
            # This block is 2 if statements in the fortran version
            mapPrint("Rearranging...  SRC Mapping:")
            # initialize temporary bundle
            if normalize:
                mapPrint("\twith normalization")
                # In fortran, just this line:
                # call cpl_bundle_initv( bunn, name, buni_local, mapx%new )
                # now in python:
                av.initv( inav_local, self.new.lgrid.lsize() )
                # Redistribute the data:
                # call cpl_mct_rearr_rearrange( buni_local%data, bunn%data, mapx%rearr )
                mapPrint("SourceAV:"+str(inav_local))
                mapPrint("TargetAV:"+str(av))
                self.rearranger.rearrange( inav_local, av, self.tag, self.sum )
            else:
                mapPrint("without normalization")
                # in fortran, just this line:
                # call cpl_bundle_initv( bunn, name, buni, mapx%new )
                av.initv( inav, self.new.lgrid.lsize() )
                # Redistribute the data:
                # call cpl_mct_rearr_rearrange( buni%data, bunn%data, mapx%rearr)
                mapPrint("SourceAV:"+str(inav))
                mapPrint("TargetAV:"+str(av))
                self.rearranger.rearrange( inav, av, self.tag, self.sum )            
            # map the data
            # CALL cpl_mct_Smat_AvMult( bunn%data, mapx%smat, buno%data )
            mapPrint("Mapping Data - Before sMatAvMult_DataLocal...")
            mapPrint("AV:"+str(av))
            mapPrint("OutAV:"+str(outav))
            mapPrint("self.sMat:"+str(self.sMat))
            mapPrint("Mapping Data - Calling sMatAvMult_DataLocal...")
            outav.sMatAvMult_DataLocal( av, self.sMat )
            mapPrint("\tDone!")
        elif (self.newtype == "dst" ):
            if( "Faii_taux" in rfields and self.firstrun ):
                inav_taux,size = inav.exportRAttr("Faii_taux")
                mapPrint("Faii_taux within inav - map.py(slice)(during mapping):",inav_taux[0:200])
                if(normalize):
                    inav_local_taux,size = inav_local.exportRAttr("Faii_taux")
                    mapPrint("Faii_taux within inav_local - map.py(slice)(during mapping):",inav_local_taux[0:200])
                
            outav = attributevector.AttributeVector( ifields, rfields,
                                                     self.destination.lgrid.lsize() )
        
            # initialize temporary bundle
            # In fortran, just this line:
            # call cpl_bundle_initv( bunn, name, buni_local, mapx%new )
            # now in python:
            mapPrint("Mapping data...  DST Mapping:")
            av = attributevector.AttributeVector( ifields, rfields,
                                                  self.new.lgrid.lsize() )
            # map the data
            if normalize:
                mapPrint("\twith normalization")
                # CALL cpl_mct_Smat_AvMult( bunn%data, mapx%smat, buno%data )
                av.sMatAvMult_DataLocal( inav_local, self.sMat )
            else:
                mapPrint("\twithout normalization")
                av.sMatAvMult_DataLocal( inav, self.sMat )
            if( "Faii_taux" in rfields and self.firstrun ):
                inav_taux,size = inav.exportRAttr("Faii_taux")
                mapPrint("Faii_taux within inav - map.py(slice)(during mapping after av.sMatAvMult_dataLocal):",inav_taux[0:200])
                av_taux,size = av.exportRAttr("Faii_taux")
                mapPrint("Faii_taux within av - map.py(slice)(during mapping after av.sMatAvMult_dataLocal):",av_taux[0:200])
                matrixelements,size = self.sMat.exportMatrixElements()
                inavfile = open("inav.txt","w")
                avfile = open("av.txt","w")
                smatfile = open("smat.txt","w")
                
                for i in xrange( inav.lsize() ):
                    inavfile.write( str(inav_taux[i])+"\n" )
                for i in xrange( av.lsize() ):
                    avfile.write( str(av_taux[i]) + "\n")
                for i in xrange( self.sMat.lsize() ):
                    smatfile.write( str(matrixelements[i]) + "\n")
                    
                inavfile.close()
                avfile.close()
            # redistribute the data
            # call cpl_mct_rearr_rearrange( buni_local%data, bunn%data, mapx%rearr )
            mapPrint("Rearranging data...")
            sum = True # this may not be necessary (RLJ)
            mapPrint("SourceAV lsize:",av.lsize() )
            mapPrint("TargetAV lsize:",outav.lsize() )
            
            self.rearranger.rearrange( av, outav, self.tag, sum )

            #if( "Faii_taux" in rfields and self.firstrun ):
            #    outav_taux,size = outav.exportRAttr("Faii_taux")
            #    mapPrint("Faii_taux within map.py(slice)(during mapping after rearranging):",outav_taux[0:200])
            #newops= self.new.gsMap.OrderedPoints( 0 )
            #debugPrint("\tnew global seg map ordered points:",newops)
            #destops = self.destination.gsMap.OrderedPoints(0)
            #debugPrint("\tdest global seg map ordered points:",destops)
            self.firstrun = False
            mapPrint("\tDone!")
        else:
            mapPrint("error, invalid newtype:"+str(self.newtype))

        ## Finish Normalization:
        if normalize:#let's skip this for now, /=frac causes an overflow exception
            # call cpl_bundle_divide(buno, bunfd, fdname)
            outav.divide( avfd, fdname )
                
        return outav
    # End of mapAV

"""
6+ npFix* methods:
"""
def npFix( av_in, av_out_tmp, domain_i, domain_o, field1, field2 ):
    """
    Fix NP Values with respect to mapping vector fields
    ! !IROUTINE: cpl_map_npFixNew3 - correct the north pole mapping of velocity fields
    !
    ! !DESCRIPTION:
    !    Correct the north pole mapping of velocity fields from the atm to ocn
    !    grids.  This assumes the input grid is a regular lat/lon with the north
    !    pole surrounded by the last latitude line in the input array.  The
    !    longitudes in the last latitude must be ordered and equally spaced.
    !
    !    This version (New3) is the same as New except it saves the gGrid data
    !    type from the first call.  This assumes the gGrid used in all calls
    !    to npfix is the same.  It is different from New2 in that it doesn't
    !    use gather to compute npu and npv.  This is bfb with New2 and New on 
    !    9/1/2003.        
    """
    root_pid = 0
    ###
    global NPFIX_FIRSTCALL
    if( NPFIX_FIRSTCALL ):
        global global_gathered_grid 
        global_gathered_grid.initv( domain_i.lgrid, domain_i.gsMap.lsize(comm.local_comm) )
        print "cpl.map.npFix: computing bilinear weights & indices for NP region."
        global_gathered_grid.gather( domain_i.lgrid, domain_i.gsMap, root_pid, comm.component_pid, comm.local_comm )
        global_gathered_grid.bcast( root_pid, comm.component_pid, comm.local_comm )
        NPFIX_FIRSTCALL=False
    ###
    # Create output attrvect:
    av_out = attributevector.AttributeVector()
    av_out.initv( av_out_tmp, av_out_tmp.lsize() )
    kui,kui_size = av_in.exportRAttr(field1)
    kvi,kvi_size = av_in.exportRAttr(field2)
    kuo,kuo_size = av_out.exportRAttr(field1)
    kvo,kvo_size = av_out.exportRAttr(field2)
    kin,kin_size = domain_i.lgrid.exportRAttr("index")
    klati,klati_size = domain_i.lgrid.exportRAttr("lat")
    kloni,kloni_size = domain_i.lgrid.exportRAttr("lon")
    klato,klato_size = domain_o.lgrid.exportRAttr("lat")
    klono,klono_size = domain_o.lgrid.exportRAttr("lon")

    nmin = (domain_i.ni) * (domain_i.nj - 1) + 1
    nmax = domain_i.n
    num = domain_i.ni
    ###
    # Pack the local data and gather it:
    local_data = pyarray.zeros( (3, num), pyarray.Float64 )
    npts = av_in.lsize()
    m = 0
    for n in xrange(npts):
        if( kin[n] >= nmin ):
            local_data[0][m] = kin[n]
            local_data[1][m] = kui[n]
            local_data[2][m] = kvi[n]
            m+=1
    root = 0
    global_data = mpi.allgather( local_data,
                                 num*3,
                                 mpi.MPI_DOUBLE,
                                 num*3,
                                 mpi.MPI_DOUBLE,
                                 comm.local_comm )
    #print "len(global_data)=",len(global_data)
    global_data = global_data.resize((3, num, comm.component_npe))
    # Unpack the gathered data on each node:
    m = 0
    local_data = pyarray.zeros( (3, num), pyarray.Float64 )
    for n2 in xrange(num):
        for n3 in xrange( comm.component_npe ):
            if( global_data[0][n2][n3] > 0.1 ):
                try:
                    index = nint( global_data[0][n2][n3] ) - nmin
                    local_data[0][index] = global_data[0][n2][n3]
                    local_data[1][index] = global_data[1][n2][n3]
                    local_data[2][index] = global_data[2][n2][n3]
                    m+=1
                except IndexError:
                    print index, n2, n3
    # Check for errors in the gathered/received data:
    if( m != num ):
        print "\terror in allgather!"
    for n2 in xrange(num):
        if( local_data[0][n2] < 0.1 ):
            print "\terror in allgather2!"
    ###
    ilon1 = pyarray.zeros( (num), pyarray.Float64 )
    ilon2 = pyarray.zeros( (num), pyarray.Float64 )
    ilat1 = pyarray.zeros( (num), pyarray.Float64 )
    ilat2 = pyarray.zeros( (num), pyarray.Float64 )

    klati,klati_size = global_gathered_grid.exportRAttr("lat")
    kloni,kloni_size = global_gathered_grid.exportRAttr("lon")
    latmax = klati[nmin]
    npu = 0.0
    npv = 0.0
    for n in xrange(num):
        try:
            np1 = (n%num) 
            nn1 = nmin + n - 1
            nn2 = nmin + np1 - 1
            rtmp = kloni[nn1]
            ilon1[n] = (rtmp + 360) % (360.0)
            rtmp = kloni[nn2]
            ilon2[n] = (rtmp + 360) % (360.0)
            ilat1[n] = klati[nn1]
            ilat2[n] = klati[nn2]
            if( ilon2[n] < ilon1[n] ):
                ilon2[n] = ilon2[n] + 360.0

            latmax = max( latmax, ilat1[n] )
        
            theta1 = ilon1[n] * const.deg2rad
            npu = npu + math.cos( theta1 )*local_data[1][n] - math.sin(theta1)*local_data[2][n]
            npv = npv + math.sin( theta1 )*local_data[1][n] + math.cos(theta1)*local_data[2][n]
        except IndexError:
            print n, nn1, nn2
            
    npu = npu / float(num)
    npv = npv / float(num)
    ###
    found = False
    olon, olat = 0,0
    npts = av_out.lsize()
    for m in xrange(npts):
        olat = klato[m]
        if( olat >= latmax ):
            rtmp = klono[m]
            olon = rtmp % 360.0
            n = 1
            found = False
            while ( ( n <= num) and (not found) ):
                if( ( olon >= ilon1[n] and olon < ilon2[n] ) or
                    ( olon+360.0 >= ilon1[n] and olon < ilon2[n] )):
                    np1 = n % num
                    ilat = (ilat1[n] + ilat2[n]) * 0.5
                    if( ilon2[n] == ilon1[n] ):
                        alpha = 0.5
                    elif( olon >= ilon1[n] and olon < ilon2[n] ):
                        alpha = (olon - ilon1[n]) / (ilon2[n] - ilon1[n])
                    elif( olon+360.0 >= ilon1[n] and olon < ilon2[n] ):
                        alpha = (olon+360.0 - ilon1[n]) / (ilon2[n] - ilon1[n])
                    else:
                        print "cpl.map.npFix: ERROR: olon ",olon, ilon1[n], ilon2[n]
                if(ilat >= 90.0):
                    beta = 1.0
                else:
                    beta = (olat-ilat)/(90. - ilat)
                w1 = (1.0-alpha)*(1.0-beta)
                w2 = (alpha)*(1.0-beta)
                w3 = (alpha)*(beta)
                w4 = (1.0-alpha)*(beta)

                theta1 = ilon1[n] * const.deg2rad
                theta2 = ilon2[n] * const.deg2rad

                f1 = local_data[1][n]
                f2 = local_data[1][np1]
                f3 = math.cos(theta1)*npu + math.sin(theta1)*npv
                f4 = math.cos(theta2)*npu + math.sin(theta2)*npv
                rtmp = w1 * f1 + w2*f2 + w3*f3 + w4*f4
                kuo[m] = rtmp

                f1 = local_data[2][n]
                f2 = local_data[2][np1]
                f3 = -1*math.sin(theta1)*npu + math.cos(theta1)*npv
                f4 = -1*math.sin(theta2)*npu + math.cos(theta2)*npv
                rtmp = w1 * f1 + w2*f2 + w3*f3 + w4*f4
                kvo[m] = rtmp
                found = True
            n = n+1
        if( not found ):
            print "cpl.map.npFix: ERROR: found = false ",found,m,olon,olat                        
    ###
    av_out.importRAttr(field1,kuo)
    av_out.importRAttr(field2,kvo)
    return av_out

    
