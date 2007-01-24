!===============================================================================
! CVS $Id: history_mod.F90,v 1.1.1.1 2005/02/03 22:29:01 steder Exp $
! CVS $Source: /home/cvsroot/steder/pyCPL/history_mod.F90,v $
! CVS $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: history_mod -- cpl6 main program history file creation module.
!
! !DESCRIPTION:
!    cpl6 main program history file creation module.
!
! !REVISION HISTORY:
!     2002-Sep-27 - B. Kauffman - created initial version
!
! !REMARKS:
!   This is not a generic low-level module, it is a high-level module 
!   hard-coded to particular bundle declarations and user desires wrt history
!   file content.
!
! !INTERFACE: ------------------------------------------------------------------

module history_mod

! !USES:

   use cpl_control_mod     ! control flags
   use shr_date_mod        ! date/time module
   use cpl_iocdf_mod       ! netCDF file creation
   use frac_mod            ! surface fractions
   use areafact_mod        ! area corrections
   use cpl_kind_mod        ! kinds
   use data_mod            ! lengthly data declarations/inits for main program
   use shr_sys_mod         ! wrappers to system calls
   use shr_timer_mod       ! timing utilities


   implicit none

   private ! except

! !PUBLIC TYPES: 
 
   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: history_write    ! create cpl6 history files
   public :: history_avbundleInit  ! initialize cpl6 tavg history bundles
   public :: history_avwrite  ! create cpl6 tavg history files

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !--- time average history bundles ---
   type(cpl_bundle) :: bun_avXa2c_a  ! everything recv'd from atm
   type(cpl_bundle) :: bun_avXi2c_i  ! everything recv'd from ice
   type(cpl_bundle) :: bun_avXl2c_l  ! everything recv'd from lnd
   type(cpl_bundle) :: bun_avXo2c_o  ! everything recv'd from ocn
   type(cpl_bundle) :: bun_avXr2c_r  ! everything recv'd from runoff
   type(cpl_bundle) :: bun_avXc2a_a  ! everything sent   to   atm
   type(cpl_bundle) :: bun_avXc2i_i  ! everything sent   to   ice
   type(cpl_bundle) :: bun_avXc2l_l  ! everything sent   to   lnd
   type(cpl_bundle) :: bun_avXc2o_o  ! everything sent   to   ocn

   integer(IN),parameter  :: pid0 = 0 ! root process pid = zero

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: history_write -- create desired history file 
!
! !DESCRIPTION:
!    create desired history file from data\_mod data.
!
! !REVISION HISTORY:
!     2002-Sep-27 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine history_write(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date) :: date ! date associated with bundles

!EOP

   !----- local -----
   character(CL),save      :: loc_fn = "unset" ! file name
   integer(IN),save        :: fid              ! file ID
   integer(IN)             :: t01   = -1       ! timer id
   integer(IN)             :: year             ! model date's year
   integer(IN)             :: month            ! model date's month
   integer(IN)             :: day              ! model date's day
   integer(IN)             :: sec              ! model date's seconds
   integer(IN)             :: cDate            ! model date's coded date (yymmss)


   !----- formats -----
   character(*),parameter :: F00 = "('(history_write) ',4a)"
   character(*),parameter :: F01 = "('(history_write) ',a,i8)"
   character(*),parameter :: F02 = "('(history_write) ',3a,i8.8,', ',i5,a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   if (.not. cpl_control_histNow) return

   if (t01 == -1) call shr_timer_get(t01,"history_write")

   call shr_timer_start(t01)

   !--- get date associated with the data ---
   call shr_date_getYMD  (date,year,month,day,sec)
   call shr_date_getCDate(date,cDate,sec)

   !--- construct the local file name and create the file ---
   if (loc_fn == "unset") then
      !oc_fn = "123456789+123456789+123456789.nc"
      loc_fn = "cpl6.hi.yyyy-mm-dd-sssss.nc"
      write(loc_fn( 9:12),'(i4.4)') year
      write(loc_fn(14:15),'(i2.2)') month
      write(loc_fn(17:18),'(i2.2)') day
      write(loc_fn(20:24),'(i5.5)') sec
      loc_fn = trim(cpl_control_caseName) // "." // loc_fn
      write(6,F00) 'creating new file: ',trim(loc_fn)
      call cpl_iocdf_create(loc_fn,desc=cpl_control_caseDesc)
   end if

   !--- append to an existing file ---
   write(6,F02) 'appending to file: ',trim(loc_fn),', date = ',cDate,sec,'s'

   call cpl_iocdf_open(loc_fn,fid)

   call cpl_iocdf_set64bit(cpl_control_hist64bit)

   call cpl_iocdf_append(fid,date,bun_frac_a)
   call cpl_iocdf_append(fid,date,bun_frac_i)
   call cpl_iocdf_append(fid,date,bun_frac_l)
   call cpl_iocdf_append(fid,date,bun_frac_o)

   call cpl_iocdf_append(fid,date,bun_areafact_a)
   call cpl_iocdf_append(fid,date,bun_areafact_i)
   call cpl_iocdf_append(fid,date,bun_areafact_l)
   call cpl_iocdf_append(fid,date,bun_areafact_o)
!  call cpl_iocdf_append(fid,date,bun_areafact_r)

   call cpl_iocdf_append(fid,date,con_Xa2c%bundle)
   call cpl_iocdf_append(fid,date,con_Xi2c%bundle)
   call cpl_iocdf_append(fid,date,con_Xl2c%bundle)
   call cpl_iocdf_append(fid,date,con_Xo2c%bundle)
   call cpl_iocdf_append(fid,date,con_Xr2c%bundle)

   call cpl_iocdf_append(fid,date,con_Xc2a%bundle)
   call cpl_iocdf_append(fid,date,con_Xc2i%bundle)
   call cpl_iocdf_append(fid,date,con_Xc2l%bundle)
   call cpl_iocdf_append(fid,date,con_Xc2o%bundle)
   call cpl_iocdf_append(fid,date,bun_Xc2oSNAP_o )

   !--- close the file ---
   call cpl_iocdf_close(fid)

   loc_fn = "unset"
   call shr_timer_stop(t01)

end subroutine history_write

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: history_avbundleInit -- initialize desired tavg history bundles
!
! !DESCRIPTION:
!    initialize desired time average history bundles
!
! !REVISION HISTORY:
!     2003-Mar-31 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine history_avbundleInit()

   implicit none

! !INPUT/OUTPUT PARAMETERS:

!EOP

   !----- formats -----
   character(*),parameter :: F00 = "('(history_avbundleInit) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   if (cpl_control_avhistType == "none") return

   call cpl_bundle_init(bun_avXa2c_a,"avXa2c_a",cpl_fields_a2c_fields,dom_a)
   call cpl_bundle_init(bun_avXi2c_i,"avXi2c_i",cpl_fields_i2c_fields,dom_i)
   call cpl_bundle_init(bun_avXl2c_l,"avXl2c_l",cpl_fields_l2c_fields,dom_l)
   call cpl_bundle_init(bun_avXo2c_o,"avXo2c_o",cpl_fields_o2c_fields,dom_o)
   call cpl_bundle_init(bun_avXr2c_r,"avXr2c_r",cpl_fields_r2c_fields,dom_r)
   call cpl_bundle_init(bun_avXc2a_a,"avXc2a_a",cpl_fields_c2a_fields,dom_a)
   call cpl_bundle_init(bun_avXc2i_i,"avXc2i_i",cpl_fields_c2i_fields,dom_i)
   call cpl_bundle_init(bun_avXc2l_l,"avXc2l_l",cpl_fields_c2l_fields,dom_l)
   call cpl_bundle_init(bun_avXc2o_o,"avXc2o_o",cpl_fields_c2o_fields,dom_o)

end subroutine history_avbundleInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: history_avwrite -- create desired tavg history file 
!
! !DESCRIPTION:
!    create desired history file from data\_mod data.
!
! !REVISION HISTORY:
!     2003-Mar-30 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine history_avWrite(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date) :: date ! date associated with bundles

!EOP

   !----- local -----
   character(CL),save     :: loc_fn = "unset" ! file name
   integer(IN),save       :: fid              ! file ID
   integer(IN)            :: year             ! model date's year
   integer(IN)            :: month            ! model date's month
   integer(IN)            :: day              ! model date's day
   integer(IN)            :: sec              ! model date's seconds
   integer(IN)            :: cDate            ! model date's coded date (yymmss)
   integer(IN)            :: eDay             ! model elapsed date
   REAL(R8)               :: rEDay            ! model elapsed date + fraction
   REAL(R8)               :: rSec             ! model elapsed seconds on eDay
   type(shr_date),save    :: date_start       ! date of first sample
   type(shr_date),save    :: date_end         ! date of last  sample
   type(shr_date),save    :: date_middle      ! date in middle of all samples
   integer(IN),save       :: nSamples = 0     ! number of sample accumulated
   integer(IN)            :: n                ! generic loop index
   integer(IN)            :: t01   = -1       ! timer id
   integer(IN)            :: t02   = -1       ! timer id
   type(cpl_bundle),save  :: tmp_Xi2c         ! Xi2c with 0 where ifrac = 0
   type(cpl_bundle),save  :: tmp_Xc2i         ! Xc2i with 0 where ifrac = 0
   integer(IN),save       :: k_imask,k_ifrac  ! aVect indecies for mask,frac
   integer(IN),save       :: lSize            ! local size of aVect
   real(R8)               :: iMask,iFrac      ! ice mask, ice fraction
   logical                :: first_call = .true.

   !----- formats -----
   character(*),parameter :: F00 = "('(history_avWrite) ',4a)"
   character(*),parameter :: F01 = "('(history_avWrite) ',a,i8)"
   character(*),parameter :: F02 = "('(history_avWrite) ',3a,i8.8,', ',i5,a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (cpl_control_avhistType == "none") return

   !----------------------------------------------------------------------------
   ! add data to a history file
   !----------------------------------------------------------------------------
   if ( cpl_control_avhistNow .and. nSamples>0 ) then

      if (t01 == -1) call shr_timer_get(t01,"history_avWrite")
      call shr_timer_start(t01)

      !--- form the time average ---
      call cpl_bundle_avg(bun_avXa2c_a)
      call cpl_bundle_avg(bun_avXi2c_i)
      call cpl_bundle_avg(bun_avXl2c_l)
      call cpl_bundle_avg(bun_avXo2c_o)
      call cpl_bundle_avg(bun_avXr2c_r)
      call cpl_bundle_avg(bun_avXc2a_a)
      call cpl_bundle_avg(bun_avXc2i_i)
      call cpl_bundle_avg(bun_avXc2l_l)
      call cpl_bundle_avg(bun_avXc2o_o)

      !--- compute date in middle of all time samples ---
      call shr_date_getEDay(date_start,eDay,sec)
      rEDay =          eDay + sec/86400.0
      call shr_date_getEDay(date_end  ,eDay,sec)
      rEDay = (rEday + eDay + sec/86400.0)/2.0  ! mid date

      !--- approximate this middle-date in a shr_date data-type ---
      eDay =  rEDay                  ! truncated mid date (note F90 truncation)
      rSec = (rEDay - eDay)*86400.0  ! elapes seconds on middle date
      date_middle = shr_date_initEDay(eDay,100)   ! approx middle date (sec=0)
      do n=1,100
         !--- approximated middle date has error < 864s < 15minutes ---
         call shr_date_getEDay(date_middle,eDay,sec)
         if ( sec < rSec ) then
            call shr_date_adv1step(date_middle)
         else
            exit
         end if
      end do

      !--- construct the local file name and create the file ---
      if (loc_fn == "unset") then
         if (cpl_control_avhistType(1:5) == "nyear" .or. &
             cpl_control_avhistType(1:4) == "year") then
            call shr_date_getYMD(date_start,year,month,day,sec)
            loc_fn = "cpl6.ha.yyyy.nc"
            write(loc_fn( 9:12),'(i4.4)') year
         elseif (cpl_control_avhistType(1:6) == "nmonth" .or. &
                 cpl_control_avhistType(1:5) == "month") then
            loc_fn = "cpl6.ha.yyyy-mm.nc"
            call shr_date_getYMD(date_start,year,month,day,sec)
            write(loc_fn( 9:12),'(i4.4)') year
            write(loc_fn(14:15),'(i2.2)') month
         elseif (cpl_control_avhistType(1:4) == "nday" .or. &
                 cpl_control_avhistType(1:4) == "date"  .or. &
                 cpl_control_avhistType(1:5) == "daily") then
            loc_fn = "cpl6.ha.yyyy-mm-dd.nc"
            call shr_date_getYMD(date_start,year,month,day,sec)
            write(loc_fn( 9:12),'(i4.4)') year
            write(loc_fn(14:15),'(i2.2)') month
            write(loc_fn(17:18),'(i2.2)') day
         else
            call shr_date_getYMD(date_end,year,month,day,sec)
            loc_fn = "cpl6.ha.yyyy-mm-dd-sssss.nc"
            write(loc_fn( 9:12),'(i4.4)') year
            write(loc_fn(14:15),'(i2.2)') month
            write(loc_fn(17:18),'(i2.2)') day
            write(loc_fn(20:24),'(i5.5)') sec
          endif
          loc_fn = trim(cpl_control_caseName) // "." // loc_fn
          write(6,F00) 'creating new file: ',trim(loc_fn)
          call cpl_iocdf_create(loc_fn,desc=cpl_control_caseDesc)
      end if

      !--- open, append to, and close an existing file ---
      call shr_date_getCDate(date,cDate,sec)
      write(6,F02) 'appending to file: ',trim(loc_fn),', date = ',cDate,sec,'s'

      call cpl_iocdf_open(loc_fn,fid)

      call cpl_iocdf_append(fid,date_middle,bun_avXa2c_a,date_start,date_end)
      call cpl_iocdf_append(fid,date_middle,bun_avXi2c_i,date_start,date_end)
      call cpl_iocdf_append(fid,date_middle,bun_avXl2c_l,date_start,date_end)
      call cpl_iocdf_append(fid,date_middle,bun_avXo2c_o,date_start,date_end)
      call cpl_iocdf_append(fid,date_middle,bun_avXr2c_r,date_start,date_end)

      call cpl_iocdf_append(fid,date_middle,bun_avXc2a_a,date_start,date_end)
      call cpl_iocdf_append(fid,date_middle,bun_avXc2i_i,date_start,date_end)
      call cpl_iocdf_append(fid,date_middle,bun_avXc2l_l,date_start,date_end)
      call cpl_iocdf_append(fid,date_middle,bun_avXc2o_o,date_start,date_end)
   
      call cpl_iocdf_close(fid)

      !--- zero-out the accumation bundle ---
      call cpl_bundle_zero(bun_avXa2c_a)
      call cpl_bundle_zero(bun_avXi2c_i)
      call cpl_bundle_zero(bun_avXl2c_l)
      call cpl_bundle_zero(bun_avXo2c_o)
      call cpl_bundle_zero(bun_avXr2c_r)
      call cpl_bundle_zero(bun_avXc2a_a)
      call cpl_bundle_zero(bun_avXc2i_i)
      call cpl_bundle_zero(bun_avXc2l_l)
      call cpl_bundle_zero(bun_avXc2o_o)
      nSamples = 0
      loc_fn = "unset"

      call shr_timer_stop(t01)

   end if

   !----------------------------------------------------------------------------
   ! add a sample to the accumulation 
   !----------------------------------------------------------------------------
   if (t02 == -1) call shr_timer_get(t02,"history_avWrite_accum")
   call shr_timer_start(t02)

   !--- for ice-domain bundles: where (ifrac = 0) set (data-sample = 0) ---
   if (first_call) then
      call cpl_bundle_initv(tmp_Xi2c,"tmp_Xi2c",con_Xi2c%bundle,con_Xi2c%bundle%dom)
      call cpl_bundle_initv(tmp_Xc2i,"tmp_Xc2i",con_Xc2i%bundle,con_Xi2c%bundle%dom)

      k_iMask = cpl_mct_aVect_indexRA(tmp_Xi2c%dom%lGrid,'mask'    )
      k_ifrac = cpl_mct_aVect_indexRA(tmp_Xi2c%data     ,'Si_ifrac')
      lSize   = cpl_mct_aVect_lSize  (tmp_Xi2c%data)
      first_call = .false.
   end if
   call cpl_bundle_fcopy(con_Xi2c%bundle,tmp_xi2c)
   call cpl_bundle_fcopy(con_Xc2i%bundle,tmp_xc2i)
   tmp_xi2c%cnt = 1  ! set number of samples to 1 (should bundle_fcopy do this?)
   tmp_xc2i%cnt = 1  ! set number of samples to 1 (should bundle_fcopy do this?)
   do n = 1,lSize
      iMask = tmp_Xi2c%dom%lGrid%rAttr(k_iMask,n) 
      iFrac = tmp_Xi2c%data     %rAttr(k_iFrac,n) 
      if (iMask /=0 .and. iFrac < 0.001) then
         tmp_Xi2c%data%rAttr(:,n) = 0.0
         tmp_Xc2i%data%rAttr(:,n) = 0.0
      end if
   enddo

   !--- accumulate samples ---
   nSamples = nSamples + 1
   date_end = date
   if (nSamples == 1) date_start = date
   call cpl_bundle_accum(con_Xa2c%bundle,outbun=bun_avXa2c_a)
   call cpl_bundle_accum(tmp_Xi2c       ,outbun=bun_avXi2c_i)
   call cpl_bundle_accum(con_Xl2c%bundle,outbun=bun_avXl2c_l)
   call cpl_bundle_accum(con_Xo2c%bundle,outbun=bun_avXo2c_o)
   call cpl_bundle_accum(con_Xr2c%bundle,outbun=bun_avXr2c_r)
   call cpl_bundle_accum(con_Xc2a%bundle,outbun=bun_avXc2a_a)
   call cpl_bundle_accum(tmp_Xc2i       ,outbun=bun_avXc2i_i)
   call cpl_bundle_accum(con_Xc2l%bundle,outbun=bun_avXc2l_l)
!  call cpl_bundle_accum(con_Xc2o%bundle,outbun=bun_avXc2o_o)
   call cpl_bundle_accum(bun_Xc2oSNAP_o ,outbun=bun_avXc2o_o)
   
   call shr_timer_stop(t02)

end subroutine history_avWrite

!===============================================================================

end module history_mod
