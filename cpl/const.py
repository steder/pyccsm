"""
cpl6.const

Translation from fortran of cpl_const_mod.F90
"""

import shr

default_tag = 600

pi       =  shr.const.PI     
rearth   =  shr.const.REARTH 
rearth2  =  shr.const.REARTH*shr.const.REARTH 
g        =  shr.const.G      
deg2rad  =  pi/180.0#_R8  
rad2deg  =  180.0/pi  
cpdair   =  shr.const.CPDAIR  
cpwv     =  shr.const.CPWV    
cpvir    =  cpwv/cpdair - 1.0#_R8 

zvir     =  shr.const.ZVIR    
latvap   =  shr.const.LATVAP  
latice   =  shr.const.LATICE  
stebol   =  shr.const.STEBOL  
karman   =  shr.const.KARMAN  
ocn_ref_sal  =  shr.const.OCN_REF_SAL 
ice_ref_sal  =  shr.const.ICE_REF_SAL 
spval        =  shr.const.SPVAL       
