import comm

def init( name ):
    """
    This routine establishes coupled component communicator
    groups, coupled component id's, and pid's relative to world
    and component communicator groups
    
    local_communicator = init( "my_component_name_in_mph" )
    """
    return comm.init( name )

def finalize( ):
    """
    Terminates the coupling/mpi environment
    """
    return

def contractInit( my_name, other_name, fields, ibuf=None, buf=None, ibufr=None, bunname=None, decomp=None ):
    """
    Initializes a contract and returns a contract object.
    """
    decomp_type = 1
    if decomp:
        decomp_type = decomp
    if ibufi:
        
