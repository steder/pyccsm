"""
Infobuffer.py

Python version of cpl_infobuf_mod.F90
"""
# Python Modules

# External Modules
import Numeric
import mpi

# Internal Modules
import fields

import debug
debugPrint = debug.newDebugPrint(False)

class Infobuffer:
    """
    Defines an Infobuffer class that stores configuration information
    and parameters for accessing the states and fluxes of attribute 
    vectors.
    
    Please note that all the send a receive operations are really 
    2 seperate operations, one for sending each field of the infobouffer.
    """
    integer_size = fields.ibuf["total"]
    real_size = fields.rbuf["total"]
    ilist = fields.ibuf.keys()
    rlist = fields.rbuf.keys()
    def __init__(self):
        # Create Buffers and initialize to 0
        #print "integer_size=%s,real_size=%s"%(self.integer_size,self.real_size)
        self.ibuf = Numeric.zeros(self.integer_size, Numeric.Int32)
        self.rbuf = Numeric.zeros(self.real_size, Numeric.Float32)
        #print self.ibuf
        #print self.rbuf
        return
    
    def getFields(self):
        """
        Returns a tuple containing two lists:  the first
        list represents all the integer buffer fields.

        The second list represents all the real (or floating point)
        fields.
        """
        ifields = fields.ibuf.keys()
        rfields = fields.rbuf.keys()
        ifields.remove("total")
        rfields.remove("total")
        return (ifields,rfields)



    def ibufGet( self, key ):
        """
        Takes a single argument, a string corresponding
        to one of the values defined in the
        cpl.fields module.

        For example:
        total size of(int)infobuffers = myinfobuf.ibufGet("total")
        """
        if( type(key) == type("") ):
            try:
                return self.ibuf[ fields.ibuf[key] ]
            except KeyError:
                raise KeyError,"Unknown Integer Buffer Key %s"%(key)
        else:
            raise AttributeError,"Key must be a string.  Look at cpl.fields for defined values."

    def ibufGetList( self, keys ):
        if( type(keys) == type([]) ):
            values  = []
            for key in keys:
                values.append( self.ibufGet( key ) )
            return values
        else:
            raise AttributeError,"Key must be a list of strings.  Look at cpl.fields for defined values."

    def ibufGetDict( self, keys ):
        if( type(keys) == type([]) ):
            values  = {}
            for key in keys:
                values[key] = self.ibufGet( key )
            return values
        else:
            raise AttributeError,"Key must be a list of strings.  Look at cpl.fields for defined values."
    
    def rbufGet( self, key ):
        """
        Takes a single argument, a string corresponding
        to one of the values defined in the
        cpl.fields module.

        For example:
        total size of(int)infobuffers = myinfobuf.ibufGet("total")
        """       
        if( type(key) == type("") ):
            try:
                return self.rbuf[ fields.rbuf[key] ]
            except KeyError:
                raise KeyError,"Unknown Real Buffer Key %s"%(key)
        else:
            raise AttributeError,"Key must be a string.  Look at cpl.fields for defined values."

    def rbufGetList( self, keys ):
        if( type(keys) == type([]) ):
            values  = []
            for key in keys:
                values.append( self.rbufGet( key ) )
            return values
        else:
            raise AttributeError,"Key must be a list of strings.  Look at cpl.fields for defined values."

    def rbufGetDict( self, keys ):
        if( type(keys) == type([]) ):
            values  = {}
            for key in keys:
                values[key] = self.rbufGet( key )
            return values
        else:
            raise AttributeError,"Key must be a list of strings.  Look at cpl.fields for defined values."

    def ibufSet( self, key, value ):
        """
        Takes 2 arguments:

        A key string that determines what index to set in the infobuffer,
        and an integer(or real) value to set at that infobuffer index.

        myinfobuffer.ibufSet( "gisize", 300 )
        """
        if( type(key) == type("") ):
            if( type(value) == type(0) ):
                try:
                    self.ibuf[ fields.ibuf[key] ] = value
                except KeyError:
                    raise KeyError,"Unknown Integer Buffer Key %s"%(key)
            else:
                raise AttributeError,"Value must be an integer(%s)"%(value)
        else:
            raise AttributeError,"Key must be a string.  Look at cpl.fields for defined values."

    def ibufSetList( self, keys, values ):
        if( type(keys) == type([]) ):
            if( type(values) == type([]) ):
                for key,value in zip(keys,values):
                    self.ibufSet( key, value )
            else:
                raise AttributeVector,"Values msut be a list of integers"
        else:
            raise AttributeError,"Keys must be a list of strings.  Look in cpl.fields for defined values."

    def ibufSetDict( self, dictionary ):
        if( type(dictionary) == type({}) ):
            for key,value in dictionary.iteritems():
                self.ibufSet( key, value )
        else:
            raise AttributeError,"Input must be a dictionary with string keys and integer values.  Look in cpl.fields for defined values."
    
    def rbufSet( self, key, value ):
        """
        Takes 2 arguments:

        A key string that determines what index to set in the infobuffer,
        and a real value to set at that infobuffer index.

        myinfobuffer.rbufSet( "gisize", 300.0 )
        """
        if( type(key) == type("") ):
            if( ( type(value) == type(0) ) or ( type(value)==type(1.0) ) ):
                try:
                    self.rbuf[ fields.rbuf[key] ] = value
                except KeyError:
                    raise KeyError,"Unknown Real Buffer Key %s"%(key)
            else:
                raise AttributeError,"Value must be an Real or Integer(%s)"%(value)
        else:
            raise AttributeError,"Key must be a string.  Look at cpl.fields for defined values."

    def rbufSetList( self, keys, values ):
        if( type(keys) == type([]) ):
            if( type(values) == type([]) ):
                for key,value in zip(keys,values):
                    self.rbufSet( key, value )
            else:
                raise AttributeVector,"Values msut be a list of real numbers"
        else:
            raise AttributeError,"Keys must be a list of strings.  Look in cpl.fields for defined values."

    def rbufSetDict( self, dictionary ):
        if( type(dictionary) == type({}) ):
            keys = dictionary.keys()
            values = dictionary.values()
            for key,value in dictionary.iteritems():
                self.rbufSet( key, value )
        else:
            raise AttributeError,"Input must be a dictionary with string keys and real number values.  Look in cpl.fields for defined values."
    
    def send( self, pid, tag, comm=mpi.MPI_COMM_WORLD ):
        """
        Thin wrapper around mpi.send:
        
        Takes 3 arguments, a process id(rank) to send to,
        a tag to use when sending the message, and finally the
        communicator in which to send it.
        
        The communicator defaults to mpi.WORLD, both the process id and
        tag always need to be supplied by the user.
        """
        status = mpi.send( self.ibuf, self.integer_size, mpi.MPI_INT, pid, tag, comm )
        status = mpi.send( self.rbuf, self.real_size, mpi.MPI_DOUBLE, pid, tag, comm )
        return

    def recv( self, pid=mpi.MPI_ANY_SOURCE, tag=mpi.MPI_ANY_TAG, comm=mpi.MPI_COMM_WORLD ):
        """
        Thin wrapper around mpi.recv:
        
        Takes 3 arguments, a process id(rank) to recv from,
        a tag to use when sending the message, and finally the 
        communicator
 
        recv returns a tuple containing the ibuf and rbuf values of 
        as a tuple
        
        comm is optional, it is defined as mpi.WORLD by default.
        """
        debugPrint("receiving integer infobuffer from node %s..."%(pid))
        self.ibuf = mpi.recv( self.integer_size, mpi.MPI_INT, pid, tag, comm )
        debugPrint("receiving real infobuffer from node %s..."%(pid))
        self.rbuf = mpi.recv( self.real_size, mpi.MPI_DOUBLE, pid, tag, comm )
        return (self.ibuf,self.rbuf)
        
    def bcast( self, comm=mpi.MPI_COMM_WORLD, root=0 ):
        """
        Thin wrapper around mpi.bcast
        
        Takes 2 arguments, a communicator and a root node to accumulate
        the broadcast on.
        
        The communicator defaults to mpi.WORLD, and the root 
        processor defaults to process 0.
        
        The communicator object is doublechecked for type safety,
        and then both the integer and real buffers are sent as two seperate 
        broadcasts.
        """
        #comm = mpi.communicator(comm)
        debugPrint("broadcasting integer infobuffer...")
        self.ibuf = mpi.bcast( self.ibuf, self.integer_size, mpi.MPI_INT, root, comm )
        debugPrint("broadcasting real infobuffer...")
        self.rbuf = mpi.bcast( self.rbuf, self.real_size, mpi.MPI_DOUBLE, root, comm )
        #self.ibuf = comm.bcast( self.ibuf, root )
        #self.rbuf = comm.bcast( self.rbuf, root )
        return

if __name__=="__main__":
    print "Creating an infobuffer:"
    buf = Infobuffer()
    print "Displaying Fields:"
    ifields, rfields = buf.getFields()
    buf.ibufSet("gsize",100)
    print "Integer Fields:",ifields
    print "Real Fields:",rfields

    print "Testing Broadcast..."
    import sys
    rank,size = mpi.init(len(sys.argv),sys.argv)
    buf.bcast()
    print "gsize=%s"%(buf.ibufGet("gsize"))
    mpi.finalize()
