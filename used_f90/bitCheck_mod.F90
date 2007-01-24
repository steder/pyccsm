!===============================================================================
! CVS: $Id: bitCheck_mod.F90,v 1.1.1.1 2005/02/03 22:29:00 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/bitCheck_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: bitCheck_mod -- do bit-for-bit check
!
! !DESCRIPTION:
!    write out global min, max, and average bundle field data at full precision.  
!    Such data useful evidence that two runs are or are not producing exactly 
!    the same results.  Such evidence can prove two runs are not identical, but
!    at best can only strongly suggest to runs are identical.
!
! !REVISION HISTORY:
!     2003-oct-20 - B. Kauffman - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

module bitCheck_mod

! !USES:

   use cpl_mct_mod       ! mct interface
   use cpl_comm_mod      ! mpi/mph communicator info
   use cpl_fields_mod    ! coupler/model data field indices
   use cpl_bundle_mod    ! defines bundle
   use cpl_domain_mod    ! defines domain
   use cpl_kind_mod      ! defines F90 kinds
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use shr_sys_mod       ! share system routines
   use shr_date_mod      ! defines date data-type
   use shr_mpi_mod       ! layer on MPI

   implicit none

   private ! except

! !PUBLIC TYPES:

   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: bitCheck_write       ! write bitCheck info to stdout
 
! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !--- module variables ---
   integer(IN),parameter  :: pid0 = 0 ! root process pid = zero
   character(*),parameter :: modName = "bitCheck_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: bitCheck_write -- write out global min, max, and avg field data.
!
! !DESCRIPTION:
!    write out global min, max, and average bundle field data at full precision.  
!
! !REVISION HISTORY:
!    2003-Oct-20 - B. Kauffman, initial version (similar to cpl5)
!
! !INTERFACE: ------------------------------------------------------------------

subroutine bitCheck_write(date,caseName,bun,fldName)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date)  ,intent(in) :: date     ! model date
   character(*)    ,intent(in) :: caseName ! case name
   type(cpl_bundle),intent(in) :: bun      ! bundle
   character(*)                :: fldName  ! field name text string

!EOP

   !----- local -----
   type(cpl_mct_aVect) :: gData    ! global/gathered bundle data
   type(cpl_mct_aVect) :: gDom     ! global/gathered bundle domain data
   integer(IN)         :: npts     ! number of points in a field
   integer(IN)         :: rc       ! return code
   integer(IN)         :: i        ! generic index
   integer(IN)         :: n        ! counts number of samples
   integer(IN)         :: kMask    ! aVect index for mask
   integer(IN)         :: kFld     ! aVect index for field
   real(R8)            :: T        ! data value for single cell (temperature?)
   real(R8)            :: mn       ! global min
   real(R8)            :: av       ! global average
   real(R8)            :: mx       ! global max
   integer(IN)         :: y,m,d,s  ! model date: year/month/day/sec

   !----- formats -----
   character(*),parameter :: subName = '(bitCheck_write) '
   character(*),parameter :: F00 = "('(bitCheck_write) ',4a)"
   character(*),parameter :: F01 = "('(bitCheck_write) ',3a,i4)"
   character(*),parameter :: F02 = "('(bitCheck_write) ', &
                             & i4.4,2('-',i2.2),i6,'s ',a,1x,a,1x,3es26.19)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (dbug>2) write(6,F00) 'working on bundle = ',trim(bun%name)

   !----------------------------------------------------------------------------
   ! do global gather to create non-distributed aVect (*part* of a bundle)
   !----------------------------------------------------------------------------
   call cpl_mct_aVect_gather(bun%data     ,gData,bun%dom%gsMap,pid0,cpl_comm_comp,rc)
   call cpl_mct_aVect_gather(bun%dom%lGrid,gDom ,bun%dom%gsMap,pid0,cpl_comm_comp,rc)
   
   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0) return

   !----------------------------------------------------------------------------
   ! find & print global min/max/avg
   !----------------------------------------------------------------------------
   kFld  = cpl_mct_aVect_indexRA(gData,trim(fldName),perrWith=subName)
   kMask = cpl_mct_aVect_indexRA(gDom ,"mask"       ,perrWith=subName)
   npts = cpl_mct_aVect_lsize (gData)

   mx = -1.0e30
   mn =  1.0e30
   av = 0.0
   n  = 0
   do i=1,npts
      if ( nint(gDom%rAttr(kMask,i)) /= 0) then
         T  = gData%rAttr(kFld,i)
         mx = max(mx,T)
         mn = min(mn,T)
         av = av+T
         n  = n+1
      end if
   end do
   av = av/float(n)

   call shr_date_getYMD(date,y,m,d,s) ! model time associate with data

   write(6,F02) y,m,d,s,trim(fldName),trim(caseName),mn,av,mx
   call shr_sys_flush(6)

   !----------------------------------------------------------------------------
   ! clean up global bundle
   !----------------------------------------------------------------------------
   call cpl_mct_aVect_clean(gData)
   call cpl_mct_aVect_clean(gDom )

end subroutine bitCheck_write

!===============================================================================
!===============================================================================

end module bitCheck_mod
