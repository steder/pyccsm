#!/usr/bin/env python2.4
"""
pyCPL6 Main

---


"""
import os.path
import time
import traceback
abspath = os.path.abspath(__file__)
basedir = os.path.dirname(abspath)
ROOT=basedir

import coupler

def main():
    print "pyCPL2: Start!"
    # Coupler Path, doComputations, doMapping, initializeRiverMapping,Logging)
    coupler.TIMER.start("total")
    coupler.TIMER.start("init")
    print "pyCPL2: Coupler root is",ROOT
    """
    >>> str(time.ctime()).replace(" ","")
    'WedJan1012:10:102007'
    """
    LOGROOT= os.path.join("/home/steder/exe/pyccsm/nc2",str(time.ctime()).replace(" ",""))
    try:
        if( coupler.comm.world_pid == 0 ):
            os.mkdir( LOGROOT )
    except:
        traceback.print_exc()
    print "pyCPL: NC Output path is",LOGROOT
    
    c = coupler.Coupler(os.path.join(ROOT,"cpl","cpl.nml"),True,True,True,True,LOGROOT)
    coupler.TIMER.stop("init")
    c.run(5)
    coupler.TIMER.stop("total")
    #print "StartDate:",cpl.date
    print "pyCPL2: Done!" 
    print coupler.TIMER
    
if __name__=="__main__":
    main()
