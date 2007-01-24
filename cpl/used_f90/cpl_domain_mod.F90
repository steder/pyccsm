!===============================================================================
! CVS: $Id: cpl_domain_mod.F90,v 1.1 2005/02/14 20:06:32 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/cpl/used_f90/cpl_domain_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_domain_mod -- fundamental data type definition and methods
!
! !DESCRIPTION:
!     Defines the coupler {\it domain} data types along with associated methods.  
!     The {\it domain} data types is a fundamental coupler data type.
!     A {\it domain} has both a {\it grid} and a {\it decomposition}.  
!     A decomposition is described by a {\it global seg map} (gsMap).
!
! !REVISION HISTORY:
!     2001-aug-15 - B. Kauffman - created module
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_domain_mod

! !USES:

   use shr_sys_mod    ! shared system call wrappers
   use cpl_kind_mod   ! kinds
   use cpl_mct_mod    ! MCT API
   use cpl_comm_mod   ! communicator groups, pids, etc.
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug

   implicit none

   private ! except

! !PUBLIC TYPES:

   public :: cpl_domain

   type cpl_domain
      !--- decomposition-independant data ---
      character(80) :: name   ! = "null" ! name of domain       (eg. "ocean")
      character(80) :: suffix ! = "null" ! netCDF domain suffix (eg. "o")
      integer(IN)   :: n      ! n = ni*nj ~ total number of grid pts (global)
      integer(IN)   :: ni     ! number of 2d array i indicies        (global)
      integer(IN)   :: nj     ! number of 2d array j indicies        (global)

      !--- decomposition-dependant data ---
      type(cpl_mct_aVect) :: lGrid ! grid data      
      type(cpl_mct_gsMap) :: gsMap ! global seg map (defines decomp)
   end type cpl_domain

! !PUBLIC MEMBER FUNCTIONS:

   public cpl_domain_info     ! print some info about a domain
   public cpl_domain_clean    ! clean/dealloc a domain
   public cpl_domain_compare  ! compare two domains for consistency

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP

   character(*),parameter :: modName = "cpl_domain_mod"  ! module name
   integer,parameter      :: pid0 = 0                    ! root PID

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_domain_info -- Write basic domain info to stdout for debugging.
!
! !DESCRIPTION:
!    Write basic, sanity-check domain info to stdout for debugging.
!
! !REVISION HISTORY:
!     2001-Dec-20 - B. Kauffman -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_domain_info(cpl_domain_x)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_domain)   ,target,intent(in) :: cpl_domain_x  ! domain

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName = "(cpl_domain_info) "
   character(*),parameter :: F00 = '("(cpl_domain_info) ",60a)'
   character(*),parameter :: F01 = '("(cpl_domain_info) ", a,i9,2i6)'
   character(*),parameter :: F02 = '("(cpl_domain_info) ", a,i3,2i6)'

!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

   write(6,F00) "<--- domain data type info dump START --->"

   write(6,F00) "name = "    ,trim(cpl_domain_x%name  ),&
                ", suffix = ",trim(cpl_domain_x%suffix)
   write(6,F01) "n,ni,nj ="  ,cpl_domain_x%n,cpl_domain_x%ni,cpl_domain_x%nj
   if (dbug >= 2) then
      write(6,F02) "gsMap comp_id, ngSeg, gSize = ", &
      cpl_domain_x%gsMap%comp_id, &
      cpl_domain_x%gsMap%ngseg,   &
      cpl_domain_x%gsMap%gsize
   end if

   !--- write local grid info ---
   if (dbug == 1) call cpl_mct_aVect_info(1,cpl_domain_x%lGrid,cpl_comm_comp,cpl_comm_comp_pid)
!  if (dbug == 2) call cpl_mct_aVect_info(2,cpl_domain_x%lGrid,cpl_comm_comp,cpl_comm_comp_pid)
   if (dbug >= 2) call cpl_mct_aVect_info(4,cpl_domain_x%lGrid,cpl_comm_comp,cpl_comm_comp_pid)

   write(6,F00) "<--- domain data type info dump END   --->"

   call shr_sys_flush(6)

end subroutine cpl_domain_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_domain_clean -- clear a domain type
!
! !DESCRIPTION:
!    clear a domain type
!
! !REVISION HISTORY:
!     2002-Jan-20 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_domain_clean(dom)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_domain)   ,intent(inout) :: dom  ! domain

!EOP

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------

   dom%name   = "<null>"
   dom%suffix = "<null>"
   dom%n  = 0
   dom%ni = 0
   dom%nj = 0
   call cpl_mct_aVect_clean(dom%lGrid)
   call cpl_mct_gsMap_clean(dom%gsMap)
   dom%gsMap%comp_id = 0
   dom%gsMap%ngseg   = 0
   dom%gsMap%gsize   = 0
   nullify(dom%gsMap%start)
   nullify(dom%gsMap%length)
   nullify(dom%gsMap%pe_loc)

end subroutine cpl_domain_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_domain_compare - compare two domains
!
! !DESCRIPTION:
!    compares two domains, summarizes the differences, aborting is an option.
!
! !REVISION HISTORY:
!    2004-May-21 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_domain_compare(dom1,dom2,enforce_mask,enforce_grid, &
           &                            enforce_area,enforce_aream, enforce_all)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_domain),intent(in) :: dom1          ! domain #1
   type(cpl_domain),intent(in) :: dom2          ! domain #2
   logical,optional,intent(in) :: enforce_mask  ! abort if masks differ wrt zero/nonzero
   logical,optional,intent(in) :: enforce_grid  ! abort if grids differ by eps_grid
   logical,optional,intent(in) :: enforce_area  ! abort if area  differ by eps_area
   logical,optional,intent(in) :: enforce_aream ! abort if aream differ by eps_area
   logical,optional,intent(in) :: enforce_all   ! abort for all of the above

!EOP

   !----- local -----
   character(CL)        :: str      ! generic text string
   integer(IN)          :: n,k      ! generic index
   logical              :: enforce  ! flags abort if inconsistent
   integer(IN)          :: ndiff    ! number of points differing
   integer(IN)          :: npts     ! number of points in a field

   integer(IN)          :: npts1    ! number of points in a dom1 field
   integer(IN)          :: npts2    ! number of points in a dom1 field
   integer(IN)          :: rcode    ! return code
   type(cpl_mct_aVect)  :: gGrid1   ! global/gathered bundle data
   type(cpl_mct_aVect)  :: gGrid2   ! global/gathered bundle data
   real(R8),allocatable :: data1(:) ! temporary real vector
   real(R8),allocatable :: data2(:) ! temporary real vector
   real(R8),allocatable :: mask (:) ! temporary real vector, domain mask

   real(R8)             :: max_diff             ! maximum diff
   real(R8)             :: diff                 ! average diff
   real(R8),parameter   :: eps_mask = 1.0e-6_R8 ! epsilon for masks
   real(R8),parameter   :: eps_grid = 1.0e-2_R8 ! epsilon for grid coords
   real(R8),parameter   :: eps_area = 1.0e-1_R8 ! epsilon for areas
   real(R8)             :: x1,x2                ! temp vars wrt wrap-around lat

   !----- formats -----
   character(*),parameter :: subName = '(cpl_domain_compare) '
   character(*),parameter :: F00 = "('(cpl_domain_compare) ',4a)"
   character(*),parameter :: F01 = "('(cpl_domain_compare) ',a,i6,a)"
   character(*),parameter :: F02 = "('(cpl_domain_compare) ',a,es10.3,a)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   write(6,F00) "domain #1 name: ",dom1%name
   write(6,F00) "domain #2 name: ",dom2%name

   !----------------------------------------------------------------------------
   ! do global gather to create non-distributed aVect (*part* of a bundle)
   !----------------------------------------------------------------------------
   call cpl_mct_aVect_gather(dom1%lGrid,gGrid1,dom1%gsMap,pid0,cpl_comm_comp,rcode)
   call cpl_mct_aVect_gather(dom2%lGrid,gGrid2,dom2%gsMap,pid0,cpl_comm_comp,rcode)
   
   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0) return

   !----------------------------------------------------------------------------
   ! compare size of domain
   !----------------------------------------------------------------------------
   npts1  = cpl_mct_aVect_lsize (gGrid1)
   npts2  = cpl_mct_aVect_lsize (gGrid2)
   npts   = npts1

   if (npts1 == npts2) then
      write(6,F01) "the domain size is = ", npts
   else
      write(6,F01) "domain size #1 = ", npts1
      write(6,F01) "domain size #2 = ", npts2
      write(6,F00) "ERROR: domain size mis-match"
      call shr_sys_abort(subName // "ERROR: domain size mis-match")
   end if

   allocate(data1(npts))
   allocate(data2(npts))
   allocate(mask (npts))

   !----------------------------------------------------------------------------
   ! compare domain mask
   !----------------------------------------------------------------------------

   call cpl_mct_aVect_getRAttr(gGrid1,"mask",data1,rcode)
   call cpl_mct_aVect_getRAttr(gGrid2,"mask",data2,rcode)

   ndiff = 0
   do n=1,npts
      if ( (abs(data1(n)) > eps_mask) .and. (abs(data1(n)) < eps_mask) .or. &
      &    (abs(data1(n)) < eps_mask) .and. (abs(data1(n)) > eps_mask) ) then
          ndiff = ndiff + 1
      end if
   end do

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_mask) .and. enforce_mask) enforce = .true. 
   if (present(enforce_all ) .and. enforce_all ) enforce = .true. 
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain masks"
      call shr_sys_abort(subName // "incompatible domain masks")
   end if

   mask(:) = data1(:) ! save mask to use below

   !----------------------------------------------------------------------------
   ! compare grid points (latitude & longitude)
   !----------------------------------------------------------------------------

   ndiff = 0

   call cpl_mct_aVect_getRAttr(gGrid1,"lat",data1,rcode)
   call cpl_mct_aVect_getRAttr(gGrid2,"lat",data2,rcode)
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         diff = abs(data1(n)-data2(n)) 
         max_diff = max(max_diff,diff)
         if ( diff > eps_grid ) ndiff = ndiff + 1
      end if
   end do
   write(6,F02) "maximum latitude  difference = ",max_diff

   call cpl_mct_aVect_getRAttr(gGrid1,"lon",data1,rcode)
   call cpl_mct_aVect_getRAttr(gGrid2,"lon",data2,rcode)
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         x1 = data1(n)
         x2 = data2(n)
         if (x1 > x2) then ! make sure x1 < x2
            x1 = data2(n)
            x2 = data1(n)
         end if
         do while ( (x1+360.0) < (x2+180.0) ) ! longitude is periodic
            x1 = x1 + 360.0
         end do
         diff = abs(x2 - x1)
         max_diff = max(max_diff,diff)
         if ( diff > eps_grid ) ndiff = ndiff + 1
      end if
   end do
   write(6,F02) "maximum longitude difference = ",max_diff

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_grid) .and. enforce_grid) enforce = .true. 
   if (present(enforce_all ) .and. enforce_all ) enforce = .true. 
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain grid coordinates"
      call shr_sys_abort(subName // "incompatible domain grid coordinates")
   end if

   !----------------------------------------------------------------------------
   ! compare area
   !----------------------------------------------------------------------------

   call cpl_mct_aVect_getRAttr(gGrid1,"area",data1,rcode)
   call cpl_mct_aVect_getRAttr(gGrid2,"area",data2,rcode)

   ndiff = 0
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         if (data2(n)/=0.0) diff = abs((data2(n)-data1(n))/data2(n)) 
         max_diff = max(max_diff,diff)
         if ( diff > eps_area ) ndiff = ndiff + 1
      end if
   end do
   write(6,F02) "maximum relative error of area (model) = ",max_diff

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_area) .and. enforce_area) enforce = .true. 
   if (present(enforce_all ) .and. enforce_all ) enforce = .true. 
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain area (model)"
      call shr_sys_abort(subName // "incompatible domain area (model)")
   end if

   !----------------------------------------------------------------------------
   ! compare aream
   !----------------------------------------------------------------------------

   call cpl_mct_aVect_getRAttr(gGrid1,"aream",data1,rcode)
   call cpl_mct_aVect_getRAttr(gGrid2,"aream",data2,rcode)

   ndiff = 0
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         if (data2(n)/=0.0) diff = abs((data2(n)-data1(n))/data2(n)) 
         max_diff = max(max_diff,diff)
         if ( diff > eps_area ) ndiff = ndiff + 1
      end if
   end do
   write(6,F02) "maximum relative error of area (map)   = ",max_diff

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_area) .and. enforce_area) enforce = .true. 
   if (present(enforce_all ) .and. enforce_all ) enforce = .true. 
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain area (map)"
      call shr_sys_abort(subName // "incompatible domain area (map)")
   end if

   !----------------------------------------------------------------------------
   ! clean-up, deallocate
   !----------------------------------------------------------------------------
   deallocate(data1)
   deallocate(data2)
   deallocate(mask )
   call cpl_mct_aVect_clean(gGrid1)
   call cpl_mct_aVect_clean(gGrid2)
   call shr_sys_flush(6)

end subroutine cpl_domain_compare

!===============================================================================
!===============================================================================

end module cpl_domain_mod
