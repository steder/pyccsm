!===============================================================================
! Cvs: $Id: areafact_mod.F90,v 1.1.1.1 2005/02/03 22:29:00 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/areafact_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
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
! !INTERFACE: ------------------------------------------------------------------

module areafact_mod

! !USES:

   use shr_sys_mod         ! share system routines
   use cpl_kind_mod        ! kinds
   use cpl_mct_mod         ! mct routines 
   use cpl_comm_mod        ! comms
   use cpl_domain_mod      ! defines domain
   use cpl_bundle_mod      ! defines bundle
   use cpl_control_mod     ! control paramters

   implicit none

   private ! except

! !PUBLIC TYPES:

  ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: areafact_init  

! !PUBLIC DATA MEMBERS:

   type(cpl_bundle),public :: bun_areafact_a   ! area corrections
   type(cpl_bundle),public :: bun_areafact_l   ! area corrections
   type(cpl_bundle),public :: bun_areafact_o   ! area corrections
   type(cpl_bundle),public :: bun_areafact_i   ! area corrections
   type(cpl_bundle),public :: bun_areafact_r   ! area corrections

   character(*),parameter :: bun_areafact_fields =  'cpl2comp:comp2cpl'

!EOP

    character(*),parameter :: modName = "areafact_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: areafact_init - initialize the area factor bundle
!
! !DESCRIPTION:
!    Initialize the area fraction bundles.  All fractions are derived from the
!    (time-invariant) component areas and the time-invariant weights areas.
!
! !REVISION HISTORY:
!     2003-Jan-02 - T. Craig, 1st version.
!
! !INTERFACE: ------------------------------------------------------------------

subroutine areafact_init(domain_a,domain_i,domain_l,domain_o,domain_r)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_domain),intent(in) :: domain_a ! domain of atm bundle
   type(cpl_domain),intent(in) :: domain_i ! domain of ice bundle
   type(cpl_domain),intent(in) :: domain_l ! domain of lnd bundle
   type(cpl_domain),intent(in) :: domain_o ! domain of ocn bundle
   type(cpl_domain),intent(in) :: domain_r ! domain of runoff bundle

!EOP

   character(*),parameter :: subName = '(areafact_init) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call areafact_set(bun_areafact_a,"areafact_a",domain_a)
   call areafact_set(bun_areafact_i,"areafact_i",domain_i)
   call areafact_set(bun_areafact_l,"areafact_l",domain_l)
   call areafact_set(bun_areafact_o,"areafact_o",domain_o)
   call areafact_set(bun_areafact_r,"areafact_r",domain_r)

   call cpl_mct_aVect_info(4,bun_areafact_a%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2a area factor')
   call cpl_mct_aVect_info(4,bun_areafact_i%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2i area factor')
   call cpl_mct_aVect_info(4,bun_areafact_l%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2l area factor')
   call cpl_mct_aVect_info(4,bun_areafact_o%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2o area factor')
   call cpl_mct_aVect_info(4,bun_areafact_r%data,cpl_comm_comp,cpl_comm_comp_pid,fld='cpl2comp',istr='c2r area factor')

end subroutine areafact_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: areafact_set - initialize the area factor bundle
!
! !DESCRIPTION:
!    Initialize an area fraction bundle.  All fractions are derived from the
!    (time-invariant) component areas and the time-invariant weights areas.
!
! !REVISION HISTORY:
!     2003-Jan-06 - T. Craig, 1st version.
!
! !INTERFACE: ------------------------------------------------------------------

subroutine areafact_set(bun,bName,dom)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout) :: bun    ! bundle of area factors to setup
   character(*)    ,intent(in)    :: bName  ! name assigned to bundle
   type(cpl_domain),intent(in)    :: dom    ! domain of atm bundle

!EOP

   character(*),parameter :: subName = '(areafact_set) '
   integer(IN) :: i1,i2,j1,j2,n,m1
   integer(IN) :: isiz1

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_bundle_init(bun,bName,bun_areafact_fields,dom)
   bun%data%rAttr = 1.0_R8
   bun%cnt = 1

   isiz1 = cpl_mct_aVect_lsize(bun%data)

   i1 = cpl_mct_aVect_indexRA(bun%data     ,"comp2cpl",perrWith=subName)
   i2 = cpl_mct_aVect_indexRA(bun%data     ,"cpl2comp",perrWith=subName)
   j1 = cpl_mct_aVect_indexRA(bun%dom%lGrid,"area"    ,perrWith=subName)
   j2 = cpl_mct_aVect_indexRA(bun%dom%lGrid,"aream"   ,perrWith=subName)
   m1 = cpl_mct_aVect_indexRA(bun%dom%lGrid,"mask"    ,perrWith=subName)

   do n=1,isiz1
      if (abs(bun%dom%lGrid%rAttr(m1,n)) >= 1.0e-06) then
         bun%data%rAttr(i1,n) = bun%dom%lGrid%rAttr(j1,n)/bun%dom%lGrid%rAttr(j2,n)
         bun%data%rAttr(i2,n) = 1.0_R8 / bun%data%rAttr(i1,n)
      endif
   enddo

end subroutine areafact_set

!===============================================================================
!===============================================================================

end module areafact_mod
