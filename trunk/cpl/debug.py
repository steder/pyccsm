"""
"""

def newDebugPrint( flag ):
    """
    Factory function that 
    defines either a null method or a method that simply
    prints it's arguments.
    """
    if( flag ):
        def debugPrint( *msgs ):
            msg = ""
            for m in msgs:
                msg += str(m)
            print msg
    else:
        def debugPrint( *msg ):
            pass
    return debugPrint

debugPrint = newDebugPrint(True)
statusPrint = newDebugPrint(True)
mapPrint = newDebugPrint(True)
