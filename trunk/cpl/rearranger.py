import MCT
import MCT.Rearranger

import attributevector
import debug
debugPrint = debug.newDebugPrint(False)

avtype = attributevector.AttributeVector

class Rearranger:
    def __init__(self):
        self.rearr = MCT.Rearranger.Rearranger()
        self._initialized = False
        return

    def init( self, domain0, domain1, comm ):
        self.rearr.init( domain0, domain1, comm )
        self._initialized = True
        return
    
    def rearrange( self, ingrid, outgrid, tag, sum_flag ):
        debugPrint( "AVs: (id, lsize, nIAttr, nRAttr)" )
        debugPrint( "SourceAV (%s,%s,%s,%s):"%(ingrid, ingrid.lsize(),
                                               ingrid.nIAttr(),ingrid.nRAttr()) )
        debugPrint( "TargetAV (%s,%s,%s,%s):"%(outgrid, outgrid.lsize(),
                                               outgrid.nIAttr(),outgrid.nRAttr()) )
        debugPrint( "tag :",tag)
        debugPrint( "sum :", sum_flag)
        if( isinstance( ingrid, avtype ) ):
            if( isinstance( outgrid, avtype )):
                debugPrint("both in and out grids are AttributeVectors")
                return self.rearr.rearrange( ingrid._av, outgrid._av, tag, sum_flag )
            else:
                debugPrint("only the in grid is an AttributeVector")
                return self.rearr.rearrange( ingrid._av, outgrid, tag, sum_flag )
        elif( isinstance(outgrid,avtype) ):
            debugPrint("only the out grid is an AttributeVector")
            return self.rearr.rearrange( ingrid, outgrid.av, tag, sum_flag )
        else:
            debugPrint("both inputs are simple AttrVects")
            return self.rearr.rearrange( ingrid, outgrid, tag, sum_flag )
            
    def __del__(self):
        if( self._initialized ):
            self.rearr.clean()
        return
    
