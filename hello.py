#!/home/steder/software/bin/pyMPI
import sys
import os
import mpi
print "Working Directory = ",os.getcwd()
rank,size = mpi.init( len(sys.argv), sys.argv )
print "Processor",rank,"thinks LD_ASSUME_KERNEL=",os.getenv("LD_ASSUME_KERNEL","UNDEFINED")
print "Processor",rank,"thinks SIDL_DLL_PATH=",os.getenv("SIDL_DLL_PATH","UNDEFINED")
print "Hello MPI World from Node",rank,"of",size
mpi.finalize()
