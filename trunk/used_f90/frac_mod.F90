!===============================================================================
! CVS: $Id: frac_mod.F90,v 1.1.1.1 2005/02/03 22:29:00 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/frac_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: frac_mod -- handles surface fractions.
!
! !DESCRIPTION:
!    Defines, declares, initializes, and updates surface fractions.
!
! !REMARKS:
!    These fractions are used for merging fields onto various domains.
!    This particular implementation of this module makes certain assumptions 
!    about which domains exist and the relationships between them.  These 
!    assumptions are hard-coded into this software implementation.
!
! ASSUMPTIONS: 
!    o atm & lnd grid cells and domain decompostion are identical
!    o ice & ocn domains are identical (includes decomposition)
!    o all atm cells are fully active
!    o all ocn cells are either fully active or fully inactive
!    o lnd cells can be partially active -- the fraction of a lnd cell that is 
!      active is the fraction that is not occupied by ocn
!    o ice cells can be partially active -- the fraction that is active is
!      determined by the ice component itself
!     
!    For each domain (atm,lnd,ice,ocn) there are four fractions: fa,fi,fl,fo,
!    three that could be used for merging, and one which indicatates the 
!    fraction of the cell which is active.
!    o merging on atm domain: Fa = fi*Fi + fl*Fl + fo*Fo   (fi + fl + fo = 1)
!    o merging on ice domain: Fi = fa*Fa + fo*Fo = Fa + Fo (fa = fo = 1, fl=0)
!    o merging on lnd domain: Fl = fa*Fa = Fa              (fa = 1, fo = fi = 0)
!    o merging on ocn domain: Fo = fa*Fa + fo*Fo           (fa + fi = 1, fl = 0)
!    o on the atm domain: fa = 1 (atm cells are fully active)
!    o on the ice domain: fi = is time-variant and determined by the ice model
!    o on the lnd domain: fl = 1 - fo and is time-invariant
!    o on the ocn domain: fo = 1 (ocn cells are fully active)
!
! !REVISION HISTORY:
!     2002-Aug-21 - B. Kauffman, 1st version.
!
! !INTERFACE: ------------------------------------------------------------------

module frac_mod

! !USES:

   use shr_sys_mod         ! shared system routines
   use shr_timer_mod       ! shared timer routines
   use shr_mpi_mod         ! shared mpi layer
   use cpl_kind_mod        ! kinds
   use cpl_comm_mod        ! mpi/mph communicator info
   use cpl_mct_mod         ! mct interface
   use cpl_const_mod       ! constants
   use cpl_domain_mod      ! defines domain
   use cpl_bundle_mod      ! defines bundle
   use cpl_map_mod         ! access to map data types and methods
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug

   implicit none

   private ! except

! !PUBLIC TYPES:

  ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: frac_init  ! one-time initialization of fraction values
   public :: frac_set   ! time-variant update of fraction values

! !PUBLIC DATA MEMBERS:

   !--- note: these could be declared in data_mod.F90            ---
   !--- or in a cpl/frac_mod.F90 & passed down from main program ---

   type(cpl_bundle),public :: bun_frac_a   ! surface fractions on atm domain
   type(cpl_bundle),public :: bun_frac_i   ! surface fractions on ice domain
   type(cpl_bundle),public :: bun_frac_l   ! surface fractions on lnd domain
   type(cpl_bundle),public :: bun_frac_o   ! surface fractions on ocn domain

   character(*),parameter :: frac_fields = 'afrac:ifrac:lfrac:ofrac'

!EOP

    character(*),parameter :: modName = "frac_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: frac_init - initialize the surface fraction bundles
!
! !DESCRIPTION:
!    Initialize the fraction bundles.  All fractions are derived from the
!    (time-invariant) ice/ocn domain masks plus the (time-variant) ice fraction.
!    This initialization routine sets the time-invariant values.
!
! !REVISION HISTORY:
!     2002-aug-21 - B. Kauffman, 1st version.
!
! !INTERFACE: ------------------------------------------------------------------

subroutine frac_init(map_o2a,domain_a,domain_i,domain_l,domain_o)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map   ),intent(inout) :: map_o2a  ! use to map ice-frac from ocn -> atm
   type(cpl_domain),intent(in   ) :: domain_a ! domain of atm fraction bundle
   type(cpl_domain),intent(in   ) :: domain_i ! domain of ice fraction bundle
   type(cpl_domain),intent(in   ) :: domain_l ! domain of lnd fraction bundle
   type(cpl_domain),intent(in   ) :: domain_o ! domain of ocn fraction bundle

!EOP

   !----- local -----
   integer(IN)            :: km          ! mask index
   integer(IN)            :: ka,ki,kl,ko ! atm,ice,lnd,ocn frac indicies

   !----- formats -----
   character(*),parameter :: subName =   '(frac_init) '
   character(*),parameter :: F00     = "('(frac_init) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! initialize fraction bundles with spvals
   !----------------------------------------------------------------------------
   call cpl_bundle_init(bun_frac_a,"frac_a",frac_fields,domain_a)
   call cpl_bundle_init(bun_frac_i,"frac_i",frac_fields,domain_i)
   call cpl_bundle_init(bun_frac_l,"frac_l",frac_fields,domain_l)
   call cpl_bundle_init(bun_frac_o,"frac_o",frac_fields,domain_o)

   bun_frac_a%data%rAttr = cpl_const_spval
   bun_frac_i%data%rAttr = cpl_const_spval
   bun_frac_l%data%rAttr = cpl_const_spval
   bun_frac_o%data%rAttr = cpl_const_spval

   bun_frac_a%cnt = 1
   bun_frac_i%cnt = 1
   bun_frac_l%cnt = 1
   bun_frac_o%cnt = 1

   !----------------------------------------------------------------------------
   ! initialize values on ice grid (based on zero ice fraction)
   !----------------------------------------------------------------------------
   bun_frac_i%data%rAttr(:,:) = 0.0_r8
   km = cpl_mct_aVect_indexRA(bun_frac_i%dom%lGrid,"mask" ,perrWith=subName)
   ka = cpl_mct_aVect_indexRA(bun_frac_i%data     ,"afrac",perrWith=subName)
   ko = cpl_mct_aVect_indexRA(bun_frac_i%data     ,"ofrac",perrWith=subName)
   where ( nint(bun_frac_i%dom%lGrid%rAttr(km,:)) /= 0 )
      bun_frac_i%data%rAttr(ka,:) = 1.0_r8
      bun_frac_i%data%rAttr(ko,:) = 1.0_r8
   end where

   !----------------------------------------------------------------------------
   ! initialize values on ocn grid (same as for ice grid)
   !----------------------------------------------------------------------------
   bun_frac_o%data%rAttr(:,:) = bun_frac_i%data%rAttr(:,:) 

   !----------------------------------------------------------------------------
   ! initialize values on atm grid (needs/assumes zero ice) 
   !----------------------------------------------------------------------------

   !--- map all fractions from ocn to atm grid ---
   call cpl_map_bun(bun_frac_o,bun_frac_a,map_o2a) 

   !--- clean up atm fraction: must be 1 everywhere ---
   ka = cpl_mct_aVect_indexRA(bun_frac_a%data,"afrac",perrWith=subName)
   bun_frac_a%data%rAttr(ka,:) = 1.0_r8

   !--- clean up ice fraction: must be 0 everywhere ---
   ki = cpl_mct_aVect_indexRA(bun_frac_a%data,"ifrac",perrWith=subName)
   bun_frac_a%data%rAttr(ki,:) = 0.0_r8 

   !--- clean up ocn fraction: must be in [0,1] ---
   ko = cpl_mct_aVect_indexRA(bun_frac_a%data,"ofrac",perrWith=subName)
   where (bun_frac_a%data%rAttr(ko,:) > 1.0_r8) bun_frac_a%data%rAttr(ko,:)=1.0_r8
   where (bun_frac_a%data%rAttr(ko,:) < 0.0_r8) bun_frac_a%data%rAttr(ko,:)=0.0_r8

   !--- compute lnd fraction: lnd = 1 - ocn, then clean it up ---
   kl = cpl_mct_aVect_indexRA(bun_frac_a%data,"lfrac",perrWith=subName)
   ko = cpl_mct_aVect_indexRA(bun_frac_a%data,"ofrac",perrWith=subName)
   bun_frac_a%data%rAttr(kl,:) = 1.0_r8 - bun_frac_a%data%rAttr(ko,:) 
   where (bun_frac_a%data%rAttr(kl,:) > 1.0_r8) bun_frac_a%data%rAttr(kl,:)=1.0_r8
   where (bun_frac_a%data%rAttr(kl,:) < 0.0_r8) bun_frac_a%data%rAttr(kl,:)=0.0_r8

   !--- NOTE: it's a requirement to elminate land points smaller than .001  ---
   !---       this is to avoid active land cells with tiny land fractions   ---
   !--- NOTE: this may result in some, presumably small, non-conservation   ---
   !--- NOTE: should probably count the number of pnts that get set to zero ---
   where (bun_frac_a%data%rAttr(kl,:) < .001) bun_frac_a%data%rAttr(kl,:) = 0.0_r8

   !----------------------------------------------------------------------------
   ! initialize values on lnd grid: same as atm grid except ice,ocn = 0 
   !----------------------------------------------------------------------------
   ki = cpl_mct_aVect_indexRA(bun_frac_l%data,"ifrac",perrWith=subName)
   ko = cpl_mct_aVect_indexRA(bun_frac_l%data,"ofrac",perrWith=subName)

   bun_frac_l%data%rAttr(: ,:) = bun_frac_a%data%rAttr(:,:)
   bun_frac_l%data%rAttr(ki,:) = 0.0_r8  ! ice frac = zero
   bun_frac_l%data%rAttr(ko,:) = 0.0_r8  ! ocn frac = zero

end subroutine frac_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: frac_set - set/update the surface fraction bundles
!
! !DESCRIPTION:
!    Set/update the fraction bundles based on input (time-variant) ice fraction 
!    and time-invariant land fraction. This set/update routine sets the 
!    time-variant values.  The companion initialization routine must be called 
!    first to set the time-invariant values.
!
! !REVISION HISTORY:
!     2002-aug-21 - B. Kauffman, 1st version.
!
! !INTERFACE: ------------------------------------------------------------------

subroutine frac_set(ifrac_i,map_o2a,domain_a,domain_i,domain_l,domain_o)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   real(R8)        ,intent(in   ) :: ifrac_i(:) ! temporary data array on atm domain
   type(cpl_map   ),intent(inout) :: map_o2a    ! use to map ifrac from ocn -> atm
   type(cpl_domain),intent(in   ) :: domain_a   ! domain of atm fraction bundle
   type(cpl_domain),intent(in   ) :: domain_i   ! domain of ice fraction bundle
   type(cpl_domain),intent(in   ) :: domain_l   ! domain of lnd fraction bundle
   type(cpl_domain),intent(in   ) :: domain_o   ! domain of ocn fraction bundle

!EOP

   !----- local -----
   integer(IN)          :: n           ! generic indicies
   integer(IN)          :: km          ! mask index
   integer(IN)          :: ka,ki,kl,ko ! atm,ice,lnd,ocn frac indicies
   real(R8)             :: minl,maxl   ! local   min/max values
   real(R8)             :: ming,maxg   ! global min/max values
   integer(IN),save     :: t01 = -1    ! shr_timer timer ID
   integer(IN)          :: lSize       ! size of local array

   !----- formats -----
   character(*),parameter :: subName =   '(frac_set) '
   character(*),parameter :: F00     = "('(frac_set) ',8a)"

!-------------------------------------------------------------------------------
! Note: we assume the rAttr indecies are the same for all fraction bundles
!       eg.  cpl_mct_aVect_indexRA(bun_frac_i%data,"ifrac") ==
!            cpl_mct_aVect_indexRA(bun_frac_a%data,"ifrac")
!-------------------------------------------------------------------------------

   if (t01 < 0) call shr_timer_get(t01,'frac_set')
   call shr_timer_start(t01)

   !----------------------------------------------------------------------------
   ! check for erroneous ice fractions 
   !----------------------------------------------------------------------------
   if (dbug > 1) then

      km = cpl_mct_aVect_indexRA(bun_frac_i%dom%lGrid,"mask" ,perrWith=subName)
      ki = cpl_mct_aVect_indexRA(bun_frac_i%data     ,"ifrac",perrWith=subName)

      maxl = 0.0_r8
      minl = 1.0_r8

      lSize = cpl_mct_aVect_lSize(bun_frac_i%data) ! size of local ice grid
      do n = 1,lSize
         if ( nint(bun_frac_i%dom%lGrid%rAttr(km,n)) /= 0 ) then
            maxl = max(maxl,bun_frac_i%data%rAttr(ki,n))
            minl = min(minl,bun_frac_i%data%rAttr(ki,n))
         end if
      end do
      call shr_mpi_min(minl,ming,cpl_comm_comp,subName)
      call shr_mpi_max(maxl,maxg,cpl_comm_comp,subName)
      if (cpl_comm_comp_pid == 0) then
         if (maxg > 1.0_r8) write(6,*) subName,"WARNING: global max ifrac = ",maxg
         if (ming < 0.0_r8) write(6,*) subName,"WARNING: global min ifrac = ",ming
      end if

   end if

   !----------------------------------------------------------------------------
   ! set/update values on ice grid, confine values into [0,1], mask=0 => frac=0
   !----------------------------------------------------------------------------

   km = cpl_mct_aVect_indexRA(bun_frac_i%dom%lGrid,"mask" ,perrWith=subName)
   ki = cpl_mct_aVect_indexRA(bun_frac_i%data     ,"ifrac",perrWith=subName)

   lSize = cpl_mct_aVect_lSize(bun_frac_i%data) ! size of local ice grid
   do n = 1,lSize
      bun_frac_i%data%rAttr(ki,n) = min(1.0_r8,ifrac_i(n))
      bun_frac_i%data%rAttr(ki,n) = max(0.0_r8,ifrac_i(n))
      if ( nint(bun_frac_i%dom%lGrid%rAttr(km,n)) == 0 ) &
         bun_frac_i%data%rAttr(ki,n) = 0.0_r8
   end do

   !----------------------------------------------------------------------------
   ! set/update values on ocn grid (assume ice & ocn have same domain)
   !----------------------------------------------------------------------------
   
   ki = cpl_mct_aVect_indexRA(bun_frac_i%data,"ifrac",perrWith=subName)
   ka = cpl_mct_aVect_indexRA(bun_frac_o%data,"afrac",perrWith=subName)

   lSize = cpl_mct_aVect_lSize(bun_frac_o%data) ! size of local ocn grid
   do n = 1,lSize
      bun_frac_o%data%rAttr(ki,n) =            bun_frac_i%data%rAttr(ki,n)
      bun_frac_o%data%rAttr(ka,n) = 1.0_r8  -  bun_frac_i%data%rAttr(ki,n)
      bun_frac_o%data%rAttr(ka,n) = min(1.0_r8,bun_frac_o%data%rAttr(ka,n))
      bun_frac_o%data%rAttr(ka,n) = max(0.0_r8,bun_frac_o%data%rAttr(ka,n))
   end do

   !----------------------------------------------------------------------------
   ! set/update values on atm grid
   !----------------------------------------------------------------------------

   !--- map ifrac onto atm grid (other mapped fracs are not usefull) ---
   call cpl_map_bun(bun_frac_o,bun_frac_a,map_o2a) 

   ka = cpl_mct_aVect_indexRA(bun_frac_a%data,"afrac",perrWith=subName)
   ki = cpl_mct_aVect_indexRA(bun_frac_a%data,"ifrac",perrWith=subName)
   kl = cpl_mct_aVect_indexRA(bun_frac_a%data,"lfrac",perrWith=subName)
   ko = cpl_mct_aVect_indexRA(bun_frac_a%data,"ofrac",perrWith=subName)

   lSize = cpl_mct_aVect_lSize(bun_frac_a%data) ! size of local atm grid
   do n = 1,lSize
      !--- restore afrac ---
      bun_frac_a%data%rAttr(ka,n) = 1.0_r8 ! must be 1.0
      !--- restore lfrac ---
      bun_frac_a%data%rAttr(kl,n) = bun_frac_l%data%rAttr(kl,n) 
      !--- clean-up ifrac ---
      bun_frac_a%data%rAttr(ki,n) = min( 1.0_r8, bun_frac_a%data%rAttr(ki,n))
      bun_frac_a%data%rAttr(ki,n) = max( 0.0_r8, bun_frac_a%data%rAttr(ki,n))
      !--- compute ofrac = 1.0 - ifrac - lfrac ---
      bun_frac_a%data%rAttr(ko,n) = 1.0_r8 - bun_frac_a%data%rAttr(ki,n) &
                                           - bun_frac_a%data%rAttr(kl,n) 
      !--- clean-up ofrac ---
      bun_frac_a%data%rAttr(ko,n) = min( 1.0_r8, bun_frac_a%data%rAttr(ko,n))
      bun_frac_a%data%rAttr(ko,n) = max( 0.0_r8, bun_frac_a%data%rAttr(ko,n))
   end do

   !----------------------------------------------------------------------------
   ! set accumulation counts to 1
   !----------------------------------------------------------------------------
   bun_frac_a%cnt = 1
   bun_frac_i%cnt = 1
   bun_frac_l%cnt = 1
   bun_frac_o%cnt = 1

   call shr_timer_stop(t01)

end subroutine frac_set

!===============================================================================
!===============================================================================

end module frac_mod
