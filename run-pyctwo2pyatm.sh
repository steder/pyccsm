#python2.4 /usr/local/pkgs/mpich2-eth-pgi-4.1/bin/mpiexec.py -l -envall -n 2 /home/steder/software/bin/python2.4 /home/steder/pyCPL/coupler.py

MPIEXEC=/usr/local/pkgs/mpich2-eth-pgi-4.1/bin/mpiexec.py 
PYTHON=/home/steder/software/bin/python2.4 
$PYTHON $MPIEXEC -l -envall -n 2 $PYTHON /home/steder/pyCPL/coupler.py : -n 1 $PYTHON /home/steder/pyCPL/xatm.py