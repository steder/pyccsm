!===============================================================================
! CVS: $Id: cpl_infobuf_mod.F90,v 1.1.1.1 2005/02/03 22:29:02 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/dead6build/cpl_infobuf_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_infobuf_mod -- information buffer module
!
! !DESCRIPTION:
!     This provides methods for the generic CCSM infobuf type.
!     The infobuf, or "information buffer", is used to exchange control flags
!     and other miscellaneous control information that is typically sent along
!     with multi-dimensional field data.
!
! !REMARKS:
!     The coupler and components use this module to carry out fundamental
!     communication.  The coupler calls this mainly from cpl_msg and the
!     components call this from cpl_interface.
!
! !REVISION HISTORY:
!     2002-Dec-5  - T. Craig - Moved cpl_coupling_ibuf methods here.
!     2003-Jan-10 - R. Jacob - change this module to work with an infobuf type.
!     2003-Jan-15 - T. Craig - Renamed this to infobuf from ibuf
!
! !INTERFACE:  -----------------------------------------------------------------

module cpl_infobuf_mod

! !USES:

   use cpl_kind_mod
   use cpl_fields_mod
   use shr_timer_mod
   use shr_sys_mod
   use shr_mpi_mod

   implicit none

   private ! except

! !PUBLIC TYPES:

   integer(IN),parameter,public :: cpl_infobuf_ibufSize = cpl_fields_ibuf_total
   integer(IN),parameter,public :: cpl_infobuf_rbufSize = cpl_fields_rbuf_total

   public :: cpl_infobuf

   type cpl_infobuf
      integer(IN) :: ibuf(cpl_infobuf_ibufSize) ! integer data
      real(R8)    :: rbuf(cpl_infobuf_rbufSize) ! real    data
   end type cpl_infobuf

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_infobuf_init   ! initialize infobuf to default values
   public :: cpl_infobuf_send   ! send an infobuf
   public :: cpl_infobuf_recv   ! recv an infobuf
   public :: cpl_infobuf_bcast  ! broadcast an infobuf

! !PUBLIC DATA MEMBERS:

   integer(IN),parameter,public :: cpl_infobuf_iDefault = 0
   integer(IN),parameter,public :: cpl_infobuf_rDefault = 0.0
  !integer(IN),parameter,public :: cpl_infobuf_ibufSize = ! must define above
  !integer(IN),parameter,public :: cpl_infobuf_rbufSize = ! must define above

!EOP

   character(*),parameter :: modName = "cpl_infobuf_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_infobuf_init -- initialize to default values
!
! !DESCRIPTION:
!    Initialize to default values, namely zero.
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2003-Jan-15 - B. Kauffman -- initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_infobuf_init(infobuf)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_infobuf), intent(out):: infobuf   ! info buffer

!EOP

   !----- formats -----
   character(*),parameter :: subName = "(cpl_infobuf_init) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   infobuf%ibuf(:) = cpl_infobuf_iDefault
   infobuf%rbuf(:) = cpl_infobuf_rDefault

end subroutine cpl_infobuf_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_infobuf_execute -- generic send/receive infobuf
!
! !DESCRIPTION:
!     Do work on infobuf array.
!     Private.
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2002-Sep-05 - T. Craig -- Combined infobuf send, recv into 
!          one subroutine
!     2002-Dec-05 - T. Craig -- added bcast to execute
!                            
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_infobuf_execute(srtype,infobuf,pid,tag,comm)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: srtype        ! 'send', 'recv', 'bcast'
   type(cpl_infobuf), intent(inout) :: infobuf          ! info buffer
   integer(IN), intent(in)      :: pid           ! proc id
   integer(IN), intent(in)      :: tag           ! tag
   integer(IN), intent(in)      :: comm          ! mpi communicator
!EOP

   !----- local -----

   integer(IN)         :: ierr                     ! mpi return code
   integer(IN),save    :: t0,t1,t2,t3,t4,t5        ! timers
   logical,save        :: first_call=.true.    

   !----- formats -----
   character(*),parameter :: subName = "(cpl_infobuf_execute) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (first_call) then
     first_call = .false.
     call shr_timer_get(t0,'cpl_infobuf_mod t0 ')
     call shr_timer_get(t1,'cpl_infobuf_mod t1 ')
     call shr_timer_get(t2,'cpl_infobuf_mod t2 ')
     call shr_timer_get(t3,'cpl_infobuf_mod t3 ')
     call shr_timer_get(t4,'cpl_infobuf_mod t4 ')
     call shr_timer_get(t5,'cpl_infobuf_mod t5 ')
   endif

   if (srtype == 'send') then
     call shr_timer_start(t0)
     call shr_mpi_send(infobuf%ibuf,pid,tag,comm,subName//" MPI in i_infobuf send")
     call shr_timer_stop(t0)
     call shr_timer_start(t1)
     call shr_mpi_send(infobuf%rbuf,pid,tag,comm,subName//" MPI in r_infobuf send")
     call shr_timer_stop(t1)
   else if (srtype == 'recv') then
     call shr_timer_start(t2)
     call shr_mpi_recv(infobuf%ibuf,pid,tag,comm,subName//" MPI in i_infobuf recv")
     call shr_timer_stop(t2)
     call shr_timer_start(t3)
     call shr_mpi_recv(infobuf%rbuf,pid,tag,comm,subName//" MPI in r_infobuf recv")
     call shr_timer_stop(t3)
   else if (srtype == 'bcast') then
     call shr_timer_start(t4)
     call shr_mpi_bcast(infobuf%ibuf,comm,subName//" MPI in ibuf bcast")
     call shr_timer_stop(t4)
     call shr_timer_start(t5)
     call shr_mpi_bcast(infobuf%rbuf,comm,subName//" MPI in ibuf bcast")
     call shr_timer_stop(t5)
   else
     write(6,*) "ERROR: ",subName," srtype unknown ",trim(srtype)
   endif

end subroutine cpl_infobuf_execute

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_infobuf_send -- generic receive infobuf
!
! !DESCRIPTION:
!     Send infobuf array
!     Calls cpl\_coupling\_infobuf
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2002-Aug-05 - T. Craig -- abstracted mpi_send call into subroutine
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_infobuf_send(infobuf,pid,tag,comm)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_infobuf), intent(inout):: infobuf                     ! info buffer
   integer(IN), intent(in)     :: pid                      ! proc id
   integer(IN), intent(in)     :: tag                      ! tag
   integer(IN), intent(in)     :: comm                     ! mpi communicator
!EOP

   !----- formats -----
   character(*),parameter :: subName = "(cpl_infobuf_send) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_infobuf_execute('send',infobuf,pid,tag,comm)

end subroutine cpl_infobuf_send

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_infobuf_recv -- generic receive infobuf
!
! !DESCRIPTION:
!     Receive infobuf array
!     Calls cpl\_coupling\_infobuf
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2002-Aug-05 - T. Craig -- abstracted mpi_recv call into subroutine
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_infobuf_recv(infobuf,pid,tag,comm)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_infobuf), intent(out):: infobuf                     ! info buffer
   integer(IN), intent(in)   :: pid                      ! proc id
   integer(IN), intent(in)   :: tag                      ! tag
   integer(IN), intent(in)   :: comm                     ! mpi communicator
!EOP

   !----- formats -----
   character(*),parameter :: subName = "(cpl_infobuf_recv) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_infobuf_execute('recv',infobuf,pid,tag,comm)

end subroutine cpl_infobuf_recv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_infobuf_bcast -- generic bcast of infobuf
!
! !DESCRIPTION:
!     Broadcast infobuf array
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2002-Aug-05 - T. Craig -- abstracted mpi_bcast call into subroutine
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_infobuf_bcast(infobuf,pid,comm)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_infobuf), intent(inout) :: infobuf                 ! integer buffer
   integer(IN), intent(in)      :: pid                  ! proc id
   integer(IN)                  :: tag                  ! tag
   integer(IN), intent(in)      :: comm                 ! mpi communicator
!EOP

   !----- formats -----
   character(*),parameter :: subName = "(cpl_infobuf_bcast) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_infobuf_execute('bcast',infobuf,pid,tag,comm)

end subroutine cpl_infobuf_bcast

!===============================================================================
!===============================================================================

end module cpl_infobuf_mod










