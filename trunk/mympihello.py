#!/home/steder/software/bin/python2.3
import sys
import os
import mpi
mpi.mpi_init(len(sys.argv),sys.argv)
rank = mpi.mpi_comm_rank( mpi.MPI_COMM_WORLD )
size = mpi.mpi_comm_size( mpi.MPI_COMM_WORLD )
print "Working Directory = ",os.getcwd()
print "Processor",rank,"thinks LD_ASSUME_KERNEL=",os.getenv("LD_ASSUME_KERNEL","UNDEFINED")
print "Processor",rank,"thinks SIDL_DLL_PATH=",os.getenv("SIDL_DLL_PATH","UNDEFINED")
print "Hello MPI World from Node",rank,"of",size
mpi.mpi_finalize()
