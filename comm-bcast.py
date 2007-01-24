import sys
import mpi

def main( local_comm ):
    name = "test"
    local_rank = mpi.comm_rank( local_comm )
    local_size = mpi.comm_size( local_comm )
    
    print "%s (%s,%s): creating root communicator!"%(name,local_rank,local_size)
    sys.stdout.flush()
    if local_rank == 0:
        tmp_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 5, 0 )
        print "%s (%s,%s): joined root communicator %s"%(name,local_rank,local_size,tmp_comm)
        sys.stdout.flush()
    else:
        tmp_comm = mpi.comm_split( mpi.MPI_COMM_WORLD, 6, 0 )
        print "%s (%s,%s): Joined non-root communicator %s"%(name,local_rank,local_size,tmp_comm)
        sys.stdout.flush()
    print "%s (%s,%s): Distributing root communicator!"%(name,local_rank,local_size)
    root_comm = mpi.bcast( tmp_comm, 1, mpi.MPI_INT, 0, local_comm )
    root_comm = root_comm[0]
    print "%s (%s,%s): Distributed root communicator!"%(name,local_rank,local_size)
    # Get total number of components and distribute to every node:
    # ncomponents = mpi.allreduce( ncomponents, 1, mpi.MPI_INT, mpi.MPI_SUM, root_comm )
    print "%s (%s,%s): Root Comm = %s"%(name,local_rank,local_size,root_comm)
    ncomponents = mpi.comm_size( root_comm )
    print "%s(%s,%s): ncomponents = %s"%(name, local_rank, local_size, ncomponents )

if __name__=="__main__":
    rank,size = mpi.init( len(sys.argv), sys.argv )
    main( mpi.MPI_COMM_WORLD )
    mpi.finalize()
