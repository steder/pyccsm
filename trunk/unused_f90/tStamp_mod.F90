!===============================================================================
! CVS: $Id: tStamp_mod.F90,v 1.2 2003/10/24 17:01:25 kauff Exp $
! CVS: $Source: /fs/cgd/csm/models/CVS.REPOS/cpl/cpl6/tStamp_mod.F90,v $
! CVS: $Name: ccsm3_0_rel04 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: tStamp_mod -- routines for logging model date and wall clock time.
!
! !DESCRIPTION:
!     Routines for logging model date, wall clock time, integration time.
!
! !REVISION HISTORY:
!     2003-Aug-28 - B. Kauffman, 1st version.
!
! !INTERFACE: ------------------------------------------------------------------

module tStamp_mod

! !USES:

   use shr_sys_mod         ! share system routines
   use cpl_kind_mod        ! kinds

   implicit none

   private ! except

! !PUBLIC TYPES:

  public :: tStamp_tic

  type tStamp_tic
     integer(IN) :: count ! value of hardware tic counter
     integer(IN) :: accum ! tic counts accumulated since initialization
     integer(IN) :: n     ! number of samples in accumulated count
  end type tStamp_tic

! !PUBLIC MEMBER FUNCTIONS:

   public :: tStamp_write  ! write the time stamp

! !PUBLIC DATA MEMBERS:

  ! no public data

!EOP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: tStamp_write - logs model date and wall clock time to stdout.
!
! !DESCRIPTION:
!     Logs model date and wall clock time to stdout and also the average and
!     instantaneous time difference between calls to this routine.  Generally it 
!     is expected that this routine is called periodically, eg. once per day, in
!     which the average and instantaneous time difference info becomes quite 
!     useful.
!
! !REMARKS:
!     If using externally declared and saved tStamp_tic data, an initial count 
!     value of less than zero implies this is the first call, thus accum and n 
!     are set zero.
!
! !REVISION HISTORY:
!     2002-aug-21 - B. Kauffman, 1st version.
!
! !INTERFACE: ------------------------------------------------------------------

subroutine tStamp_write(str,year,month,day,sec,tic_ext)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)    ,intent(in)             :: str     ! info text string
   integer(IN)     ,intent(in)             :: year    ! model year  (4-digits)
   integer(IN)     ,intent(in)             :: month   ! model month (2-digits)
   integer(IN)     ,intent(in)             :: day     ! model day   (2-digits)
   integer(IN)     ,intent(in)             :: sec     ! model secs  (5-digits)
   type(tStamp_tic),intent(inout),optional :: tic_ext ! external tic count data

!EOP

   !----- local -----
   character( 8) :: dstr             ! date string
   character(10) :: tstr             ! time string
   integer       :: tic_new          ! new system clock tic count
   integer       :: tic_old          ! previous system clock tic count
   integer       :: tic_rate         ! tics per second
   integer       :: tic_max          ! max tic count before roll-over
   integer       :: tic_diff         ! tic_new - tic_old
   integer       :: tic_accum        ! running sum of all tic_diff's
   integer       :: tic_n            ! number of samples in tic_accum
   real(R8)      :: dt               ! tic_diff   in units of seconds
   real(R8)      :: avdt             ! average dt in units of seconds

   type(tStamp_tic),save :: tic_int  ! internal tic count data

   DATA tic_int%count / -999 /
   
!! if some machine doesn't support the DATA statement, do this... (B Kauffman)
!! integer,save :: tic_int_count = -999 ! internal/saved value for tic_count
!! integer,save :: tic_int_accum        ! internal/saved value for tic_accum
!! integer,save :: tic_int_n            ! internal/saved value for tic_n

   !----- formats -----
   character(*), parameter :: F00 = "('(tStamp_write) ',4a)"
   character(*), parameter :: F10 = "('(tStamp_write) ',a, &
   &    '  model date ', i4.4,2('-',i2.2),1x,i5.5,'s',     &
   &    '  wall clock ', a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2, &
   &    '  avg dt',i6,'s  dt',i6,'s'   )"

!-------------------------------------------------------------------------------
! Notes: 
! o tic_old < 0 flags this is the 1st call to this routine, hence dt=0
! Assumptions: 
! o (new count) < (old count) => clock wrapped around once
! o (new count) < (old count) <= clock wrapped around
!-------------------------------------------------------------------------------

   !--- access saved data from internal or external data structure ---
   if (.not. present(tic_ext) ) then
      tic_old   = tic_int%count
      tic_accum = tic_int%accum
      tic_n     = tic_int%n
!!    tic_old   = tic_int_count
!!    tic_accum = tic_int_accum
!!    tic_n     = tic_int_n
   else 
      tic_old   = tic_ext%count
      tic_accum = tic_ext%accum
      tic_n     = tic_ext%n
   end if

   !--- get wall clock date/time and tic-count ---
   call date_and_time(dstr,tstr)
   call system_clock(tic_new,tic_rate,tic_max)

   !--- compute diff and accumulated clock tics ---
   if (tic_old < 0) then
      tic_diff  = 0                      ! 1st call => diff   = 0
      tic_accum = 0                      ! 1st call => accum  = 0
      tic_n     = 0                      ! 1st call => n = 0
   else
      tic_diff = tic_new - tic_old       ! compute diff count
      tic_n    = tic_n   + 1             ! increment counter
   end if
   if (tic_diff < 0) tic_diff = tic_diff + tic_max ! fix clock wrap-around
   tic_accum = tic_accum + tic_diff                ! t = t + dt

   !--- convert clock tics to secs, as appropriate ---
   if (tic_n > 0) then
      dt   = tic_diff / float(tic_rate       )  ! current diff
      avdt = tic_accum/(float(tic_rate)*tic_n)  ! time-avg diff
   else
      dt   = 0.0
      avdt = 0.0
   end if

   !--- write the time stamp ---
   write(6,F10) trim(str),year,month,day,sec, &
   &  dstr(1:4),dstr(5:6),dstr(7:8),  tstr(1:2),tstr(3:4),tstr(5:6), &
   &  nint(avdt),nint(dt)
   call shr_sys_flush(6)

   !--- save new data, either internally or externally ---
   if (.not. present(tic_ext)) then
      tic_int%count = tic_new
      tic_int%accum = tic_accum
      tic_int%n     = tic_n
!!    tic_int_count = tic_new
!!    tic_int_accum = tic_accum
!!    tic_int_n     = tic_n
   else 
      tic_ext%count = tic_new
      tic_ext%accum = tic_accum
      tic_ext%n     = tic_n
   end if

end subroutine tStamp_write

!===============================================================================

end module tStamp_mod

