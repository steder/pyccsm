import sys
import mpi
import cpl
import cpl.infobuffer

rank,size = mpi.init( len(sys.argv), sys.argv )

if ( rank == 0 ):
    ib = cpl.infobuffer.Infobuffer()
    ib.ibufSet("gsize",100)
    ib.ibufSet("gisize",10)
    ib.ibufSet("gjsize",10)
    ib.send( 1, 7, mpi.MPI_COMM_WORLD )
else:
    ic = cpl.infobuffer.Infobuffer()
    ic.recv( 0, 7, mpi.MPI_COMM_WORLD )
    print "Received gsize=%s\tgisize=%s\tgjsize=%s"%(ic.ibufGet("gsize"),ic.ibufGet("gisize"),ic.ibufGet("gjsize"))
mpi.finalize()
