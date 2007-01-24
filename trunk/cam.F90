#include <misc.h>
#include <params.h>

#if ( !defined SCAM )
program cam
!-----------------------------------------------------------------------
!
! Purpose: Entry point for NCAR CAM 
!
!-----------------------------NOTICE------------------------------------
!
!            Community Atmospheric Model, $Name: ccsm3_0_rel04 $ (CAM)
!
! 
! Method: Call appropriate initialization, time-stepping, and finalization routines.
! 
! Author: CAM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => SHR_KIND_R8
   use pmgrid
   use dycore,  only: dycore_is
   use history, only: bldfld, intht
   use units
   use restart,     only: read_restart
   use time_manager, only: get_nstep
   use abortutils, only: endrun
#if ( defined SPMD || defined COUP_CSM )
   use mpishorthand, only: mpicom, nsend, nrecv, nwsend, nwrecv
#endif
   use shr_sys_mod, only: shr_sys_flush
   use runtime_opts, only: runtime_options

#if (defined COUP_CSM)
   use shr_msg_mod
   use cpl_interface_mod
   use cpl_comm_mod
   use cpl_fields_mod
#endif
#if (defined BFB_CAM_SCAM_IOP)
   use history, only: initialize_iop_history
#endif
!-----------------------------------------------------------------------
   implicit none

#include <comctl.h>
#include <comlun.h>
#include <gpt.inc>

#ifdef SUNOS
!#include <floatingpoint.h>
#endif
!
! Local workspace
!
#ifdef OSF1
#include <for_fpe_flags.f>
   integer(4) old_fpe_flags   ! old settings of floating point exception flags
   integer(4) new_fpe_flags   ! new settings of floating point exception flags
   integer(4) for_set_fpe     ! function to set the floating point exceptions
#endif
   character*8 cdate          ! System date
   character*8 ctime          ! System time
   character*13 filenam
   integer iu
   integer :: nstep           ! Current timestep number.
!------------------------------Externals--------------------------------
#if ( defined SUNOS )
!      integer iexcept, ieee_handler
#endif
!-----------------------------------------------------------------------

#ifdef OSF1
!
! Compaq floating point exception handler 
! Terminate if hit invalid, divide by zero, or overflow.
!
  new_fpe_flags = FPE_M_TRAP_INV + FPE_M_TRAP_DIV0 + FPE_M_TRAP_OVF
  old_fpe_flags = for_set_fpe(new_fpe_flags)
#endif

#if ( defined SUNOS )
!
! SUN: Trap ieee exceptions for debugging purposes
!      iexcept = ieee_handler( 'set', 'common', SIGFPE_ABORT )
!      if ( iexcept /= 0 ) write(6,*)'ieee trapping not supported here'
!
#endif
!
! Initialize timing library.  2nd arg 0 means disable, 1 means enable
!
   call t_setoptionf (usrsys, 0)
   call t_initializef ()

   call t_startf ('total')
   call t_startf ('initialization')

!
! Initialize internal/external MPI if appropriate
!
#if ( defined COUP_CSM )
   call shr_msg_stdio('atm')
   !call cpl_interface_init(cpl_fields_atmname,mpicom)
   call cpl_comm_init_wo_mph(cpl_fields_atmname,mpicom)
#endif
!
! Set up spectral arrays
!
   call trunc
!
! Initialize SPMD environment if applicable
!
#if ( defined SPMD )
   call spmdinit ()
#endif
!
! Print Model heading and copyright message
!
   if (masterproc) then
      write(6,*)'------------------------------------------------------------'
      write(6,*)'NCAR Community Atmospheric Model (CAM)'
      write(6,*)'$Name: ccsm3_0_rel04 $ '
      write(6,*)'$Date: 2004/05/20 18:36:01 $'
      write(6,*)'------------------------------------------------------------'
      write(6,*)'(Online documentation is available on the CAM'
      write(6,*)' home page: http://www.ccsm.ucar.edu/models/atm-cam/'
      write(6,*)' License information is available as a link from above or from:'
      write(6,*)' home page: http://www.ccsm.ucar.edu/models/atm-cam/license.html)'
      write(6,*)'------------------------------------------------------------'
   end if
!
! Fetch and print current date and time
!
   call datetime(cdate,ctime)
   if (masterproc) then 
      write(6,*) 'DATE ',cdate, ' TIME ', ctime
      write(6,*)'------------------------------------------------------------'
      if (dycore_is ('EUL')) then
         write(6,*)'DYCORE is EUL'
      else if (dycore_is ('SLD')) then
         write(6,*)'DYCORE is SLD'
      else if (dycore_is ('LR')) then
         write(6,*)'DYCORE is LR'
      end if
   end if
!
! Set defaults then override with user-specified input
!
   call runtime_options ()   ! used to be called preset() and parse_namelist()

!
! Initialize SPMD decompositions if applicable
!
#if ( defined SPMD )
   call decompinit ()
#endif
!
! Define fortran unit numbers
!
   nsds    = getunit ()
   nrg     = getunit ()
   nrg2    = getunit ()
   luhrest = getunit ()

   if (masterproc) then
      write(6,*) '**** Summary of Logical Unit assignments ****'
      write(6,*)
      write(6,*) '   Restart pointer unit (nsds)     = ', nsds
      write(6,*) '   Master restart unit (nrg)       = ', nrg
      write(6,*) '   Abs/ems unit for restart (nrg2) = ', nrg2
      write(6,*) '   History restart unit (luhrest)  = ', luhrest
   end if
!
! Initialize index values for advected and non-advected tracers
!
   call initindx ()
!
! Do appropriate dynamics and history initialization depending on whether initial, restart, or 
! branch.  On restart run intht need not be called because all the info is on restart dataset.
!
   select case (nsrest)
   case (0)                ! initial run
      call inital ()       ! dynamics (mostly) init
      call inti ()         ! physics init
      call bldfld ()       ! master field list
      call intht ()        ! set up history tape contents for this run
   case (1)                ! restart
      call read_restart () ! read restart file(s)
#if (defined BFB_CAM_SCAM_IOP)
      call initialize_iop_history()
#endif
      call inti ()         ! physics init
      call bldfld ()       ! master field list
   case (3)                ! branch
      call read_restart () ! read restart file(s), minus history info
#if (defined BFB_CAM_SCAM_IOP)
      call initialize_iop_history()
#endif
      call inti ()         ! physics init
      call bldfld ()       ! master field list
      call intht ()        ! set up history tape contents for this run
   case default
      write(6,*)'CAM: nsrest=', nsrest,' must be 0, 1, or 3'
      call endrun ()
   end select
!
! Initialize external models or datasets depending upon whether coupled
!
   call initext ()
   call t_stopf ('initialization')
!
! Invoke driving routine for time integration
!
   call t_startf('stepon')
   call stepon ()
   call t_stopf('stepon')
!
! End the run cleanly
!
   call t_stopf('total')
   call t_prf(iam)

#if ( defined SPMD )
   if (.false.) then
      write(0,*)'The following stats are exclusive of initialization/boundary datasets'
      write(0,*)'Number of messages sent by proc ',iam,' is ',nsend
      write(0,*)'Number of messages recv by proc ',iam,' is ',nrecv
   end if
#endif

! This flush attempts to ensure that asynchronous diagnostic prints from all 
! processes do not get mixed up with the "END OF MODEL RUN" message printed 
! by masterproc below.  The test-model script searches for this message in the 
! output log to figure out if CAM completed successfully.  This problem has 
! only been observed with the Linux Lahey compiler (lf95) which does not 
! support line-buffered output.  
   call shr_sys_flush( 6 )   ! Flush all output to standard output

   if (masterproc) then
      nstep = get_nstep()
      write (6,9300) nstep-1,nstep
9300  format (//'Number of completed timesteps:',i6,/,'Time step ',i6, &
                ' partially done to provide convectively adjusted and ', &
                'time filtered values for history tape.')
      write(6,*)'------------------------------------------------------------'
      write(6,*)'******* END OF MODEL RUN *******'
   end if

#if ( defined COUP_CSM )
   call cpl_interface_finalize(cpl_fields_atmname)
#elif ( defined SPMD )
   call mpibarrier (mpicom)
   call mpifinalize
#endif

#if ( defined SPMD )
   iu = getunit ()
   write(filenam,'(a10,i3.3)') 'spmdstats.', iam
   open (unit=iu, file=filenam, form='formatted', status='replace')
   write (iu,*)'iam ',iam,' msgs  sent =',nsend
   write (iu,*)'iam ',iam,' msgs  recvd=',nrecv
   write (iu,*)'iam ',iam,' words sent =',nwsend
   write (iu,*)'iam ',iam,' words recvd=',nwrecv
#endif
!   call print_memusage
   stop 0
end program cam
#else
subroutine cam
end subroutine cam
#endif    !  defined SCAM

