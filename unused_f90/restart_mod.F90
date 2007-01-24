!===============================================================================
! CVS $Id: restart_mod.F90,v 1.13 2004/04/02 17:10:43 kauff Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/cpl/cpl6/restart_mod.F90,v $
! CVS $Name: ccsm3_0_rel04 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: restart_mod -- cpl6 main program read/write restart file module.
!
! !DESCRIPTION:
!    cpl6 program read/write restart file module.
!
! !REVISION HISTORY:
!     2002-Nov-06 - B. Kauffman - created initial version
!
! !REMARKS:
!   This is not a generic low-level module, it is a high-level module 
!   hard-coded to particular bundle declarations and restart data needs for a
!   particular version of the coupler.
!
! !INTERFACE: ------------------------------------------------------------------

module restart_mod

! !USES:

   use cpl_kind_mod        ! access to F90 kind declarations
   use cpl_control_mod, dbug=>cpl_control_infoDBug
   use cpl_comm_mod        ! mpi communicator groups & related
   use diag_mod            ! runtime diagnostic subsystem module
   use shr_date_mod        ! date/time module
   use cpl_iobin_mod       ! binary data file creation
   use frac_mod            ! surface fractions
   use data_mod            ! lengthly data declarations/inits for main program
   use shr_sys_mod         ! wrappers to system calls
   use shr_timer_mod       ! timing utilities
   use shr_file_mod        ! file get/put
   use shr_mpi_mod         ! mpi layer

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: restart_write    ! write a restart file
   public :: restart_read     ! read  a restart file
   public :: restart_readDate ! read  a restart file, read date only

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !--- module variables ---
   integer(IN),parameter  :: pid0 = 0 ! root process pid = zero

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: restart_write -- create desired restart file.
!
! !DESCRIPTION:
!    create desired restart file, update restart pointer file.
!
! !REMARKS:
!    Acceses data file from module data_mod .
!
! !REVISION HISTORY:
!     2002-Nov-06 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine restart_write(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date) :: date ! date associated with bundles

!EOP

   !----- local -----
   character(CL),save :: loc_fn = "unset" ! file name
   character(CL)      :: str              ! generic text string
   character( 8)      :: dstr             ! F90 wall clock date str yyyymmdd
   character(10)      :: tstr             ! F90 wall clock time str hhmmss.sss
   integer(IN),save   :: fid              ! file ID (a fortran unit?)
   integer(IN)        :: t01   = -1       ! timer id
   integer(IN)        :: year             ! model date's year
   integer(IN)        :: month            ! model date's month
   integer(IN)        :: day              ! model date's day
   integer(IN)        :: sec              ! model date's seconds
   integer(IN)        :: cDate            ! model date's coded date (yymmss)
   integer(IN)        :: rCode            ! return code
   logical            :: open             ! true if file unit is open
   integer(IN)        :: nUnit            ! a file unit number
   real(R8)           :: data(8*6*3)      ! diagnostic data ~ must coordinate
   character(CL)      :: dataName         ! diagnostic data, data ID/name
   integer(IN)        :: i,j,k,n          ! generic indicies

   !----- formats -----
   character(*),parameter :: subname = "(restart_write)"
   character(*),parameter :: F00 = "('(restart_write) ',4a)"
   character(*),parameter :: F01 = "('(restart_write) ',a,i8)"
   character(*),parameter :: F02 = "('(restart_write) ',3a,i8.8,', ',i5,a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   if (.not. cpl_control_restNow) return

   if (dbug > 1) then
      write(6,F00) "cpl_control_caseName  = ", cpl_control_caseName 
      write(6,F00) "cpl_control_restType  = ", cpl_control_restType
      write(6,F01) "cpl_control_restCDate = ", cpl_control_restCDate
      write(6,F00) "cpl_control_restPFn   = ", cpl_control_restPFn  
      write(6,F00) "cpl_control_restBFn   = ", cpl_control_restBFn  
      call shr_sys_flush(6)
   end if

   call shr_mpi_barrier(cpl_comm_comp,subName)

   if (t01 == -1) call shr_timer_get(t01,"restart_write")
   call shr_timer_start(t01)

   !--- get date associated with the data ---
   call shr_date_getYMD  (date,year,month,day,sec)
   call shr_date_getCDate(date,cDate,sec)

   !----------------------------------------------------------------------------
   ! construct the local file name and create the file 
   !----------------------------------------------------------------------------
   loc_fn     =  "unset"  ! => create new file for each write
   if (loc_fn == "unset") then
      !oc_fn = "123456789+123456789+123456789.nc"
      loc_fn = "cpl6.r.yyyy-mm-dd-sssss "
      write(loc_fn( 8:11),'(i4.4)') year
      write(loc_fn(13:14),'(i2.2)') month
      write(loc_fn(16:17),'(i2.2)') day
      write(loc_fn(19:23),'(i5.5)') sec
      loc_fn = trim(cpl_control_caseName) // "." // loc_fn

      write(6,F00) 'creating new file: ',trim(loc_fn)
      call cpl_iobin_create(loc_fn,trim(cpl_control_caseDesc))
      call shr_sys_flush(6)
   end if

   !----------------------------------------------------------------------------
   ! open an existing file & append data
   !----------------------------------------------------------------------------
   write(6,F02) 'appending to file: ',trim(loc_fn),', date = ',cDate,sec,'s'
   call shr_sys_flush(6)

   fid = -999
   if ( cpl_comm_comp_pid == pid0 ) call cpl_iobin_open(loc_fn,fid)

   !--- append bundle data ---
   call cpl_iobin_appendBun(fid,date,con_Xa2c%bundle) ! everything from atm
   call cpl_iobin_appendBun(fid,date,con_Xi2c%bundle) ! everything from ice
   call cpl_iobin_appendBun(fid,date,con_Xl2c%bundle) ! everything from lnd
   call cpl_iobin_appendBun(fid,date,con_Xo2c%bundle) ! everything from ocn
   call cpl_iobin_appendBun(fid,date,con_Xr2c%bundle) ! everything from run-off
   call cpl_iobin_appendBun(fid,date,con_Xc2o%bundle) ! everything to   ocn
   call cpl_iobin_appendBun(fid,date,bun_aoflux_o   ) ! atm/ocn fluxes calc by cpl
   call cpl_iobin_appendBun(fid,date,bun_oalbedo_o  ) ! ocn albedos    calc by cpl

   !--- append t-avg diagnostic data, master process only ---
   if ( cpl_comm_comp_pid == pid0 ) then
      data(1)  = diag_ns
      data(2)  = diag_eday0
      data(3)  = diag_eday1
      dataName = "diag_ns "
      call cpl_iobin_appendReal(fid,date,dataName,data(1:3),rcode)

      !--- somewhat kludgy re-shaping of diag_datas wrt appendReal API ---
      n=0
      do i=1,8
      do j=1,6
      do k=1,3
         n=n+1
         data(n) = diag_datas(i,j,k)
      end do
      end do
      end do

      dataName = "diag_datas "
      call cpl_iobin_appendReal(fid,date,dataName,data,rcode)

   end if

   !--- close the file ---
   if ( cpl_comm_comp_pid == pid0 ) call cpl_iobin_close(fid)

   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe only
   !----------------------------------------------------------------------------
   call shr_timer_stop(t01)
   if (cpl_comm_comp_pid /= pid0 ) return
   call shr_timer_start(t01)

   !----------------------------------------------------------------------------
   ! update the restart pointer file
   !----------------------------------------------------------------------------
   write(6,F00) 'updating restart pointer file: ',trim(cpl_control_restPFn)
   call shr_sys_flush(6)

   !--- find an unused unit number ---
   do nUnit=10,99
      inquire(fid,opened=open)
      if (.not. open) exit
   end do
   if (open) then
      write(6,F00) 'ERROR: couldn''t find unused unit number'
      call shr_sys_abort(subName//": all units open?")
   end if

   !--- create a local pointer file ---
   call date_and_time(dstr,tstr)
   open (nUnit,file='rpointer',form="FORMATTED",status="REPLACE")
   write(nUnit,'(a)') trim(loc_fn)
   write(nUnit,'(a)') trim(cpl_control_caseName)
   write(nUnit,'(a)') trim(cpl_control_caseDesc)
   write(nUnit,*) 'Pointer file created: '                    &
      &     //dstr(1:4)//'-'//dstr(5:6)//'-'//dstr(7:8)//' '  &
      &     //tstr(1:2)//':'//tstr(3:4)//':'//tstr(5:6)
   close(nUnit)

   !--- copy pointer file elsewhere, if necessary ---
   call shr_file_put(rcode,"rpointer",trim(cpl_control_restPFn))

   call shr_sys_flush(6)
   call shr_timer_stop(t01)

end subroutine restart_write

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: restart_read -- read desired date from restart file 
!
! !DESCRIPTION:
!    Read desired date from restart file.
!
! !REVISION HISTORY:
!     2002-Nov-06 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine restart_read(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(inout) :: date   ! date associated with restart data

!EOP

   !----- local -----
   character(CL),save :: loc_fn = "unset" ! file name
   character(CL)      :: str              ! generic text string
   integer(IN),save   :: fid              ! file ID (a fortran unit?)
   integer(IN)        :: t01   = -1       ! timer id
   integer(IN)        :: rCode            ! return code
   integer(IN)        :: i,j,k,n          ! generic indicies
   logical            :: open             ! true iff file/unit is open
   real(R8)           :: data(8*6*3)      ! diagnostic data ~ must coordinate
   character(CL)      :: dataName         ! diagnostic data, data ID/name

   !----- formats -----
   character(*),parameter :: subName = '(restart_read) '
   character(*),parameter :: F00 = "('(restart_read) ',4a)"
   character(*),parameter :: F01 = "('(restart_read) ',a,i3)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   if (t01 == -1) call shr_timer_get(t01,"restart_read")
   call shr_timer_start(t01)

   if ( trim(cpl_control_restType) == "initial" ) then
      !-------------------------------------------------------------------------
      ! no data file to read 
      !-------------------------------------------------------------------------

      write(6,F00) 'restart type = ',trim(cpl_control_restType), &
      &            ' => no restart file to read'

   else if ( trim(cpl_control_restType) == "continue" .or. &
   &         trim(cpl_control_restType) == "branch"        ) then 
      !-------------------------------------------------------------------------
      ! read data found in a restart file
      !-------------------------------------------------------------------------

      !--- determine restart data file name ---
      if (  trim(cpl_control_restType) == "branch" ) then
         write(6,F00) 'restart type = ',trim(cpl_control_restType), &
         &            ' => restart file specified by namelist variable'
         loc_fn = cpl_control_restBFn
      else
         write(6,F00) 'restart type = ',trim(cpl_control_restType), &
         &            ' => restart file specified by pointer file'

         if (cpl_comm_comp_pid == pid0 ) then
            call shr_file_get(rcode,"rpointer",trim(cpl_control_restPFn))
         end if
         call shr_sys_flush(6)
         call shr_mpi_barrier(cpl_comm_comp,subName)

         !--- find an unused unit number ---
         do fid=10,99
            inquire(fid,opened=open)
            if (.not. open) exit
         end do
         if (open) then
            write(6,F00) 'ERROR: couldn''t find unused unit number'
            call shr_sys_abort(subName//": all file unit numbers in use?")
         end if
         if (dbug>2) write(6,F01) 'using unit number ',fid

         !--- get restart file name from pointer file ---
         open (fid,file="rpointer",STATUS='OLD',IOSTAT=rcode)
         if (rcode /= 0) then
            write(6,F00) 'ERROR: opening file = rpointer'
            call shr_sys_abort(subName//": ERROR openint rpointer file")
         end if
         read (fid,*) loc_fn
         close(fid)
      end if

      !--- open restart file ---
      write(6,F00) 'open restart file: ',trim(loc_fn)
      call shr_sys_flush(6)
      call cpl_iobin_open(loc_fn,fid)

      !--- read the file, skipping over some bundles, as necessary ---
      if (cpl_control_icData_a) then
         write(6,F00) '* NOT reading atm IC data from cpl restart file'
      else
         call cpl_iobin_readBun(fid,date,con_Xa2c%bundle) ! everything from atm
      end if
      if (cpl_control_icData_i) then
         write(6,F00) '* NOT reading ice IC data from cpl restart file'
      else
         call cpl_iobin_readBun(fid,date,con_Xi2c%bundle) ! everything from ice
      end if
      if (cpl_control_icData_l) then
         write(6,F00) '* NOT reading lnd IC data from cpl restart file'
      else
         call cpl_iobin_readBun(fid,date,con_Xl2c%bundle) ! everything from lnd
      end if
      if (cpl_control_icData_o) then
         write(6,F00) '* NOT reading ocn IC data from cpl restart file'
      else
         call cpl_iobin_readBun(fid,date,con_Xo2c%bundle) ! everything from ocn
      end if
      if (cpl_control_icData_r) then
         write(6,F00) '* NOT reading roff IC data from cpl restart file'
      else
         call cpl_iobin_readBun(fid,date,con_Xr2c%bundle) ! everything from roff
      end if

      call cpl_iobin_readBun(fid,date,con_Xc2o%bundle) ! everything to ocn
      call cpl_iobin_readBun(fid,date,bun_aoflux_o   ) ! atm/ocn fluxes
      call cpl_iobin_readBun(fid,date,bun_oalbedo_o  ) ! ocn albedos

      !--- read t-avg diagnostic data, master process only ---
      if ( cpl_comm_comp_pid == pid0 ) then

         dataName = "diag_ns "
         call cpl_iobin_readReal(fid,date,dataName,data(1:3),rcode)
         diag_ns    = data(1)
         diag_eday0 = data(2)
         diag_eday1 = data(3)

         dataName = "diag_datas "
         call cpl_iobin_readReal(fid,date,dataName,data(:)  ,rcode)
         !--- somewhat kludgy re-shaping of diag_datas wrt readReal API ---
         n=0
         do i=1,8
         do j=1,6
         do k=1,3
            n=n+1
            diag_datas(i,j,k) = data(n)
         end do
         end do
         end do

      end if

      !--- close the file ---
      call cpl_iobin_close(fid)

   else 
      !-------------------------------------------------------------------------
      ! unrecognized option
      !-------------------------------------------------------------------------
      write(6,F00) 'ERROR: unrecognized restart type: ',trim(cpl_control_restType)
      call shr_sys_abort(subName//": ERROR restart type")
   end if

   call shr_timer_stop(t01)
   call shr_sys_flush(6)

end subroutine restart_read

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: restart_readDate -- read model date from restart file 
!
! !DESCRIPTION:
!    Read model date from restart file.
!
! !REMARKS:
!    All processors need to read the date from the restart file.
!
! !REVISION HISTORY:
!     2002-Dec-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine restart_readDate(cDate)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(out) :: cDate  ! start date from restart file

!EOP

   !----- local -----
   character(CL),save :: loc_fn = "unset" ! file name
   character(CL)      :: str              ! generic text string
   integer(IN),save   :: fid              ! file ID (a fortran unit?)
   integer(IN)        :: t01   = -1       ! timer id
   integer(IN)        :: rCode            ! return code
   logical            :: open             ! true iff file/unit is open
   integer(IN)        :: sec              ! model date's seconds

   !----- formats -----
   character(*),parameter :: subName = '(restart_readDate) '
   character(*),parameter :: F00 = "('(restart_readDate) ',4a)"
   character(*),parameter :: F01 = "('(restart_readDate) ',a,i3)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   if (t01 == -1) call shr_timer_get(t01,"restart_readDate")
   call shr_timer_start(t01)

   if ( trim(cpl_control_restType) == "initial" ) then
      !-------------------------------------------------------------------------
      ! no file to read 
      !-------------------------------------------------------------------------

      write(6,F00) 'restart type = ',trim(cpl_control_restType), &
      & ' => start date specified by input namelist'

      cDate = cpl_control_restCDate

   else if ( trim(cpl_control_restType) == "continue" &
   &   .or.  trim(cpl_control_restType) == "branch"   ) then
      !-------------------------------------------------------------------------
      ! read data found in a restart file
      !-------------------------------------------------------------------------

      !--- determine restart data file name ---
      if (  trim(cpl_control_restType) == "branch" ) then
         write(6,F00) 'restart type = ',trim(cpl_control_restType), &
         &            ' => restart file specified by namelist variable'
         loc_fn = cpl_control_restBFn
      else
         write(6,F00) 'restart type = ',trim(cpl_control_restType), &
         &            ' => restart file specified by pointer file'

         if (cpl_comm_comp_pid == pid0 ) then
            call shr_file_get(rcode,"rpointer",trim(cpl_control_restPFn))
         end if
         call shr_sys_flush(6)
         call shr_mpi_barrier(cpl_comm_comp,subName)

         !--- find an unused unit number ---
         do fid=10,99
            inquire(fid,opened=open)
            if (.not. open) exit
         end do
         if (open) then
            write(6,F00) 'ERROR: couldn''t find unused unit number'
            call shr_sys_abort(subName//": all file unit numbers in use?")
         end if
         if (dbug>2) write(6,F01) 'using unit number ',fid

         !--- get restart file name from pointer file ---
         open (fid,file="rpointer",STATUS='OLD',IOSTAT=rcode)
         if (rcode /= 0) then
            write(6,F00) 'ERROR: opening file = rpointer'
            call shr_sys_abort(subName//": ERROR openint rpointer file")
         end if
         read (fid,*) loc_fn
         close(fid)
      end if

      !--- open, read, & close the restart file ---
      write(6,F00) 'reading start date from restart file: ',trim(loc_fn)
      call shr_sys_flush(6)
      call cpl_iobin_open(loc_fn,fid)
      call cpl_iobin_readDate(fid,cDate,sec,rcode)
      call cpl_iobin_close(fid)

   else 
      !-------------------------------------------------------------------------
      ! unrecognized option
      !-------------------------------------------------------------------------
      write(6,F00) 'ERROR: unrecognized restart type: ',trim(cpl_control_restType)
      call shr_sys_abort(subName//": ERROR restart type")
   end if

   call shr_timer_stop(t01)

end subroutine restart_readDate

!===============================================================================
!===============================================================================

end module restart_mod
