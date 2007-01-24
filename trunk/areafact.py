"""
! !MODULE: areafact_mod -- handles surface fractions.
!
! !DESCRIPTION:
!    Defines, declares, initializes, and updates area correction bundles.
!
! !REMARKS:
!    These are used to correct flux fields received and sent to components
!    based on their areas and the areas in the weights files.
!
! !REVISION HISTORY:
!     2003-Jan-02 - T. Craig, 1st version.
!     2003-Jan-06 - T. Craig, moved work to areafact_set, removed use of cpl_map
!
! !INTERFACE: --
"""

#
import math

#
import Numeric

#
import cpl.comm
import cpl.domain
import cpl.attributevector
import cpl.control

#
from cpl.debug import *

debugPrint = newDebugPrint(False)

av_areafact_a = 0
av_areafact_l = 1
av_areafact_o = 2
av_areafact_i = 3
av_areafact_r = 4

areafact_fields = ["cpl2comp","comp2cpl"]

def init( domain_a, domain_i, domain_l, domain_o, domain_r ):
    """
    ! !IROUTINE: areafact_init - initialize the area factor bundle
    !
    ! !DESCRIPTION:
    !    Initialize the area fraction bundles.  All fractions are derived from the
    !    (time-invariant) component areas and the time-invariant weights areas.
    !
    ! !REVISION HISTORY:
    !     2003-Jan-02 - T. Craig, 1st version.
    !
    ! !INTERFACE: --
     
    ! !USES:
    
    ! !INPUT/OUTPUT PARAMETERS:
    
    type(cpl_domain),intent(in) :: domain_a ! domain of atm bundle
    type(cpl_domain),intent(in) :: domain_i ! domain of ice bundle
    type(cpl_domain),intent(in) :: domain_l ! domain of lnd bundle
    type(cpl_domain),intent(in) :: domain_o ! domain of ocn bundle
    type(cpl_domain),intent(in) :: domain_r ! domain of runoff bundle
    
    !EOP
    
    character(*),parameter :: subName = '(areafact_init) '
    """
    global av_areafact_a,av_areafact_i,av_areafact_l
    global av_areafact_o,av_areafact_r
    debugPrint("areafact.py: Initializing areafactor atm:")
    av_areafact_a = set(av_areafact_a,domain_a)
    debugPrint("areafact.py: Initializing areafactor ice:")
    av_areafact_i = set(av_areafact_i,domain_i)
    debugPrint("areafact.py: Initializing areafactor lnd:")
    av_areafact_l = set(av_areafact_l,domain_l)
    debugPrint("areafact.py: Initializing areafactor ocn:")
    av_areafact_o = set(av_areafact_o,domain_o)
    debugPrint("areafact.py: Initializing areafactor river:")
    av_areafact_r = set(av_areafact_r,domain_r)
    debugPrint("areafact.py: Initialized!")
    #   call cpl_mct_aVect_info(4,bun_areafact_a%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2a area factor')
    #   call cpl_mct_aVect_info(4,bun_areafact_i%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2i area factor')
    #   call cpl_mct_aVect_info(4,bun_areafact_l%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2l area factor')
    #   call cpl_mct_aVect_info(4,bun_areafact_o%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2o area factor')
    #   call cpl_mct_aVect_info(4,bun_areafact_r%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2r area factor')
    return


def set(av, dom):
    """
    ! !IROUTINE: areafact_set - initialize the area factor bundle
    !
    ! !DESCRIPTION:
    !    Initialize an area fraction bundle.  All fractions are derived from the
    !    (time-invariant) component areas and the time-invariant weights areas.
    !
    ! !REVISION HISTORY:
    !     2003-Jan-06 - T. Craig, 1st version.
    !
    ! !INTERFACE: 
    
    subroutine areafact_set(bun,bName,dom)

    ! !USES:
    
    ! !INPUT/OUTPUT PARAMETERS:
    
    type(cpl_bundle),intent(inout) :: bun    ! bundle of area factors to setup
    character(*)    ,intent(in)    :: bName  ! name assigned to bundle
    type(cpl_domain),intent(in)    :: dom    ! domain of atm bundle

    """
    #character(*),parameter :: subName = '(areafact_set) '
    #integer(IN) :: i1,i2,j1,j2,n,m1
    #integer(IN) :: isiz1

    #call cpl_bundle_init(bun,bName,bun_areafact_fields,dom)
    av = cpl.attributevector.AttributeVector([],areafact_fields,dom.lsize())
    #bun%data%rAttr = 1.0_R8
    # initialize av to 1.0:
    for f in areafact_fields:
        debugPrint("areafact.py: importing field %s"%(f))
        av.importRAttr(f,Numeric.ones((dom.lsize()),Numeric.Float64))
    
    #isiz1 = cpl_mct_aVect_lsize(av%data)
    isize = av.size()
    
    #i1 = cpl_mct_aVect_indexRA(av%data     ,"comp2cpl",perrWith=subName)
    #i2 = cpl_mct_aVect_indexRA(av%data     ,"cpl2comp",perrWith=subName)
    debugPrint("areafact.py: exporting comp2cpl")
    i1,i1_size = av.exportRAttr("comp2cpl")
    debugPrint("areafact.py: exporting cpl2comp")
    i2,i2_size = av.exportRAttr("cpl2comp")
    #j1 = cpl_mct_aVect_indexRA(av%dom%lGrid,"area"    ,perrWith=subName)
    #j2 = cpl_mct_aVect_indexRA(av%dom%lGrid,"aream"   ,perrWith=subName)
    #m1 = cpl_mct_aVect_indexRA(av%dom%lGrid,"mask"    ,perrWith=subName)
    debugPrint("areafact.py: exporting area")
    j1,j1_size = dom.lgrid.exportRAttr("area")
    debugPrint("areafact.py: exporting aream")
    j2,j2_size = dom.lgrid.exportRAttr("aream")
    debugPrint("areafact.py: exporting mask")
    m1,m1_size = dom.lgrid.exportRAttr("mask")
    #do n=1,isiz1
    for n in xrange(isize):
        if (abs(m1[n]) >= 1.0e-06):
            i1[n] = j1[n]/j2[n]
            i2[n] = 1.0 / i1[n]

    debugPrint("areafact.py: importing comp2cpl")
    av.importRAttr("comp2cpl",i1)
    debugPrint("areafact.py: importing cpl2comp")
    av.importRAttr("cpl2comp",i2)

    return av

