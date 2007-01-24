#! /usr/bin/env sh
python=$(which python)
mpiexec=$(which mpiexec)
mpiscript="$mpiexec.py"

if [[ $# -ne 2 ]]; then
    echo "usage:"
    echo "\t./run.sh <Number of Processors> <(parallel)python script>"
    echo "ex:"
    echo "\t./run.sh 2 send-recv.py"
else
#echo $1 $2 $mpiscript
    echo "run.sh: starting $2 on $1 processors..."
    echo "run.sh: using $python"
    echo "run.sh: using $mpiscript"
    $python $mpiscript -l -n $1 $python $PWD/$2
fi