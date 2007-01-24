!===============================================================================
! CVS $Id: timeCheck.F90,v 1.6 2003/05/30 17:26:49 kauff Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/cpl/cpl6/timeCheck.F90,v $
! CVS $Name: ccsm3_0_rel04 $
!===============================================================================
!BOP ===========================================================================
!
! !ROUTINE: timeCheck -- verify/enforce component time coordination.
!
! !DESCRIPTION:
!    Verify/enforce component time coordination.
!
! !REMARKS: HISTORY:
! date/time for component models is "global" data accessed from module data_mod.
!
! !REVISION HISTORY:
!     2002-Nov-20 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine timeCheck(date_c,print,enforce)

! !USES:

   use cpl_kind_mod        ! access to F90 kind declarations
   use shr_cal_mod         ! access to calendar routines
   use shr_date_mod        ! access to date     routines
   use data_mod            ! lengthly data declarations/inits for main program
   use shr_sys_mod         ! wrappers to system calls

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in)  :: date_c  ! official coupler date
   logical       ,intent(in)  :: print   ! true => print out the cpl date
   logical       ,intent(in)  :: enforce ! true => abort if time uncoordinated

!EOP

   !----- local -----
   integer(IN) :: sec           ! model date's seconds
   integer(IN) :: cDate         ! model date's coded date (yymmss)
   integer(IN) :: sec_a         ! atm infobuf date: seconds
   integer(IN) :: sec_i         ! ice infobuf date: seconds
   integer(IN) :: sec_l         ! lnd infobuf date: seconds
   integer(IN) :: sec_o         ! ocn infobuf date: seconds
   integer(IN) :: cDate_a       ! atm infobuf date: coded date (yymmss)
   integer(IN) :: cDate_i       ! ice infobuf date: coded date (yymmss)
   integer(IN) :: cDate_l       ! lnd infobuf date: coded date (yymmss)
   integer(IN) :: cDate_o       ! ocn infobuf date: coded date (yymmss)
   integer(IN) :: y,m,d         ! year,month,day
   logical :: uncoordinated ! flags components uncoordinated in time

   !----- formats -----
   character(*),parameter :: F00 = "('(timeCheck) ',4a)"
   character(*),parameter :: F01 = "('(timeCheck) ',a,' date ',&
   &                         i4.4,2('-',i2.2),', ',i5,' sec')"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !--- get the date & time for all components ---
   call shr_date_getCDate(date_c,cDate,sec)
   cDate_a = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_cDate)
   cDate_i = con_Xi2c%infobuf%ibuf(cpl_fields_ibuf_cDate)
   cDate_l = con_Xl2c%infobuf%ibuf(cpl_fields_ibuf_cDate)
   cDate_o = con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_cDate)
   sec_a   = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_sec)
   sec_i   = con_Xi2c%infobuf%ibuf(cpl_fields_ibuf_sec)
   sec_l   = con_Xl2c%infobuf%ibuf(cpl_fields_ibuf_sec)
   sec_o   = con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_sec)

   !--- check for time coordination ---
   uncoordinated = .false.
   if (cdate /= cDate_a) uncoordinated = .true.
   if (sec   /=   sec_a) uncoordinated = .true.
   if (cdate /= cDate_i) uncoordinated = .true.
   if (sec   /=   sec_i) uncoordinated = .true.
   if (cdate /= cDate_l) uncoordinated = .true.
   if (sec   /=   sec_l) uncoordinated = .true.
   if (cdate /= cDate_o) uncoordinated = .true.
   if (sec   /=   sec_o) uncoordinated = .true.

   !--- print out information ?? ---
   if ( uncoordinated ) then
      call shr_date_getYMD(date_c,y,m,d,sec)
      write(6,F01) "cpl:",y,m,d,sec
      call shr_cal_date2ymd(cdate_a,y,m,d)
      write(6,F01) "atm:",y,m,d,sec_a
      call shr_cal_date2ymd(cdate_i,y,m,d)
      write(6,F01) "ice:",y,m,d,sec_i
      call shr_cal_date2ymd(cdate_l,y,m,d)
      write(6,F01) "lnd:",y,m,d,sec_l
      call shr_cal_date2ymd(cdate_o,y,m,d)
      write(6,F01) "ocn:",y,m,d,sec_o
      write(6,F00) "WARNING: component model dates are uncoordinated"
   else if ( print) then
      call shr_date_getYMD(date_c,y,m,d,sec)
      write(6,F01) "cpl:",y,m,d,sec
   end if
   call shr_sys_flush(6)

   !--- abort if dates are uncoordinated? ---
   if ( enforce .and. uncoordinated ) then
      write(6,F00)'ERROR: models uncoordinated in time'
      call shr_sys_flush(6)
      call shr_sys_abort('(timeCheck) models uncoordinated in time')
   end if

end subroutine timeCheck

!===============================================================================
!===============================================================================
