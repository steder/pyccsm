"""
This script starts a MPI2 program properly so that it will
work with PyMCT.
"""

import sys, os, string

PYTHON="python2.4"
mpis=os.popen("which mpiexec").read()
mpis=mpis.strip(" ")
mpis=mpis.strip("\n")
MPIEXEC_SCRIPT=mpis+".py"

if(len(sys.argv)>=3):
    NUMPROCS=sys.argv[1]
    PROGRAM=sys.argv[2]
    EXEC=PYTHON+" "+MPIEXEC_SCRIPT+" -l -n "+NUMPROCS+" "+PYTHON+" "+PROGRAM
    output=os.popen( EXEC ).read()
    print output
else:
    print "Provide the script you want to run as an argument to run.py"
    print "i.e.:"
    print "$ python2.4 run.py twocmp-con.py"
    print "\n\nMPIEXEC_SCRIPT:",MPIEXEC_SCRIPT
    sys.exit()

