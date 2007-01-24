"""
control.py

Python version of mod.F90

Defines the following routines:

readNameList
init
update

This module represents a mmajor  subsystem of pyCPL.  It contains
data type definitions and associated methods used for
controlling coupler integration.  Here "controlling" referes to issues
such as:
    * Selecting integration start and stop dates
    * Determinging when history and restart data files should be made
    * Determinging when diagnostic data should be generated
    * Verifying that all coupled system component models( eg. atm, ocn, etc)
      are synchronized in time.
    * Reads and parses namelist variables
    * makes simulation control variables available to other modules.
"""

########## Public Interface

# Do we need python namelist support?  A method for reading / parsing
# namelist files?  Should be easy enough to implement, just need a 
# namelist file to work with.
#public :: readNList  # read & parse namelist input

# Ignore this because we don't need to concern ourselves with timings 
# until we get everything else working
#public :: init       # initialize alarms, etc.

# Unnecessary, we can access and update control flags directly
# by just modifying the modules public members.
#public :: update     # update control flags

########## Control Flags

def init( date ):
    """
    I believe that the defaults are currently okay.  If we need to have
    an initialization similar to the comm module this is where that 
    belongs.
    """
    return

stopnow = False   # T => stop model now
stopeod = False   # T => stop model at end of day
restnow = False   # T => create restart data now
resteod = False   # T => create restart data at EOD

histnow = False   # T => create history data now
histeod = False   # T => create history data at EOD
#  histSave = False  # T => archive history data now
hist64bit = False # T => use 64 bit netCDFfiles 
avhistnow = False # T => create history data now
avhisteod = False # T => create history data at EOD

diagnow = False  # T => print diagnostic data now
diageod = False   # T => print diagnostic data at eod
avdiagnow = False # T => print tavg diag data now
avdiageod = False # T => print tavg diag data at eod
bfbflag = False   # T => bfb with different pes

casename = ""  # case name
casedesc = ""  # case description

resttype = ""  # restart type: init,cont,branch
restcdate = 0  # restart cDate from namelist
restdate = 0  # restart date
restpfn = ""   # restart pointer file name
restbfn = ""   # restart branch  file name

lagocn = False    # T => lag the ocn at startup
sendatmalb = False # T => send albedo ICs to atm
sendlnddom = False # T => send lnd domain to lnd
icdata_a = False   # T => use IC data provided by atm
icdata_i = False   # T => use IC data provided by ice
icdata_l = False   # T => use IC data provided by lnd
icdata_o = False   # T => use IC data provided by ocn
icdata_r = False   # T => use IC data provided by roff
avhisttype = "" # tavg history file type

mapfn_a2of = "" # map data file: a->o fluxes
mapfn_a2os = "" # map data file: a->o states
mapfn_o2af = "" # map data file: o->a fluxes
mapfn_o2as = "" # map data file: o->a states
mapfn_r2o  = "" # map data file: r->o runoff

fluxalbav  = False # T => NO diurnal cycle in albedos
fluxepbal  = "Unset!" # selects E,P,R adjustment technique
fluxepfac  = 0.0 # E,P,R adjust factor recv'd from ocn
fluxashift = 0 # albedo calc time-shift (seconds)

orbeccen   = 1.670772e-2 # eccen of earth orbit (unitless)
orbobliqr  = 4.091238e-01 # earth's obliquity (rad)
orblambm0  = -3.250364e-02 # mean lon perihelion @ vernal eq (rad)
orbmvelpp  = 4.935568e0 # moving vernal equinox longitude

dead_a     = False # T => atm component is dead comp
dead_i     = False # T => ice component is dead comp
dead_l     = False # T => lnd component is dead comp
dead_o     = False # T => ocn component is dead comp
dead_ao    = False # T => atm and/or ocn are dead comp

ncpl_a     = 0 # atm/cpl communications per day
ncpl_i     = 0 # ice/cpl communications per day
ncpl_l     = 0 # lnd/cpl communications per day
ncpl_o     = 0 # ocn/cpl communications per day
ncpl_r     = 0 # rof/cpl communications per day

cdate_a    = 0 # atm coded date
cdate_i    = 0 # ice coded date
cdate_l    = 0 # lnd coded date
cdate_o    = 0 # ocn coded date

sec_a      = 0 # atm secs on date
sec_i      = 0 # ice secs on date
sec_l      = 0 # lnd secs on date
sec_o      = 0 # ocn secs on date

# Default = 1, 2 is a special decomp, and 903 is a testing invalid decomp
decomp_a   = 1 # atm decomposition type
decomp_l   = 1 # lnd decomposition type
decomp_o   = 1 # ocn decomposition type
decomp_i   = 1 # ice decomposition type
decomp_r   = 1 # rof decomposition type

infodbug = 1 # user specified dbug level
infobcheck = False # T => do bit-check now

def readNamelist( path ):
    """
    Namelist files are simple files with a variable name on the left and 
    the value of that variable on the right.  

    mybool = .false.
    myreal = 999.0
    mystring = 'etc etc , you get the picture'
    """
    file = open(path,"r")
    lines = file.readlines()
    file.close()
    # Replace literal " and ' in input
    lines = [ l.replace("\'","") for l in lines ]
    lines = [ l.replace("\"","") for l in lines ]
    # Drop the first and last lines
    lines = [ l for l in lines[1:len(lines)-1] ]
    # Get the assignments
    pairs = [ ( line.split("=")[0].strip(), line.split("=")[1].strip() )  for line in lines ] 
    vars = {}
    for x in pairs:
        # Convert to the appropriate type
        value = x[1]
        try:
            value = float(x[1])
        except:
            pass
        try:
            value = int(x[1])
        except ValueError:
            pass
        if( x[1].upper() == ".FALSE." ):
            value = False
        if( x[1].upper() == ".TRUE." ):
            value = True
        vars[x[0]] = value
    return vars
    
def update( date ):
    """
    This method updates the date that the model thinks it is and causes the 
    model to check if it is time to finish it's run or keep going.

    If an alarm is set to go off on a specific date or after a specific number
    of days, this is the same code that handles turning on and shutting off
    alarms.
    """
    return
