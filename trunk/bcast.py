import sys,mpi
rank,size = mpi.init(len(sys.argv),sys.argv)
sigma = mpi.bcast( 7, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD )
print "broadcasted result:",sigma 
mpi.finalize()
