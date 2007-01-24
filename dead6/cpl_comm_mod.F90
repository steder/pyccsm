!===============================================================================
! CVS: $Id: cpl_comm_mod.F90,v 1.1.1.1 2005/02/03 22:29:01 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/dead6/cpl_comm_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_comm_mod -- "communicator module" defines comm groups and model ID's
!
! !DESCRIPTION:
!     Sets up communicator groups and MPH component ID's (cid).
!     Also declares and defines handy data wrt number of pe's, PID's, CID's
!     relative to world and component communicator groups.
!
! !REVISION HISTORY:
!     2001-Aug-20 - B. Kauffman - new naming convention
!     2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_comm_mod

! !USES:

  use cpl_kind_mod   ! kinds
  use shr_sys_mod    ! system calls
  use shr_mpi_mod    ! mpi layer

  implicit none

  private ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_comm_init_mike

! !PUBLIC DATA MEMBERS:

   integer(IN),public :: cpl_comm_wrld         ! = MPI_COMM_WORLD, global comm grp
   integer(IN),public :: cpl_comm_wrld_npe     ! number of pe's in MPI_COMM_WORLD
   integer(IN),public :: cpl_comm_wrld_pid     ! this comp pid in MPI_COMM_WORLD

   integer(IN),public :: cpl_comm_comp         ! this comp communicator group
   integer(IN),public :: cpl_comm_comp_npe     ! number of pe's in comp comm group
   integer(IN),public :: cpl_comm_comp_pid     ! this comp's pid in comp comm group

   integer(IN),public :: cpl_comm_mph_cid      ! MPH component ID, this component
   integer(IN),public :: cpl_comm_mph_cid_atm  ! MPH component ID, atm
   integer(IN),public :: cpl_comm_mph_cid_ice  ! MPH component ID, ice
   integer(IN),public :: cpl_comm_mph_cid_lnd  ! MPH component ID, lnd
   integer(IN),public :: cpl_comm_mph_cid_ocn  ! MPH component ID, ocn
   integer(IN),public :: cpl_comm_mph_cid_cpl  ! MPH component ID, cpl

   integer(IN),public :: cpl_comm_wrld_pe0     ! comm world pe0, this component
   integer(IN),public :: cpl_comm_wrld_pe0_atm ! comm world pe0, atm
   integer(IN),public :: cpl_comm_wrld_pe0_ice ! comm world pe0, ice
   integer(IN),public :: cpl_comm_wrld_pe0_lnd ! comm world pe0, lnd
   integer(IN),public :: cpl_comm_wrld_pe0_ocn ! comm world pe0, ocn
   integer(IN),public :: cpl_comm_wrld_pe0_cpl ! comm world pe0, cpl

!EOP

   save

   character(*),parameter :: modName = 'cpl_comm_mod'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_comm_init -- initialize the coupling/mpi environment.
!
! !DESCRIPTION:
!    This routine calls mpi\_init, establishes coupled component communicator
!    groups, coupled component id's, and pid's relative to world and component
!    communicator groups.
! 
! !REMARKS:
!
! !REVISION HISTORY:
!     2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob -- first prototype
!     2001-Dec-10 - R. Jacob -- switch arguments in cpl_mct_world_init to
!                   to match new version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_comm_init_mike(name,comm)

! !USES:
   use cpl_fields_mod      ! contains valid component name strings
   use mph_module,only : mph_components
   use mph_module,only : mph_global_proc_id
   use mph_module,only : mph_local_proc_id
   use mph_module,only : mph_total_components
   use mph_module,only : mph_comp_id
   use mph_module,only : mph_local_totprocs
   use mph_module,only : mph_global_totprocs
   use mph_module,only : mph_global_id
   use mph_module,only : mph_comp_name
   use mph_module,only : mph_global_world
   use m_MCTWorld    ,only : cpl_mct_world_init   => init

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: name ! name of component name
   integer(IN) ,intent(out) :: comm ! communicator group for component

!EOP

   integer(IN)      :: i, n    ! generic loop index
   integer,parameter :: MAXNAMELEN = 256
   integer(IN) :: color
!  List of Names
   character, allocatable, dimension(:,:) :: names, unique_names
!  MPI Variables
   integer(IN)      :: myid, numprocs, status(MPI_STATUS_SIZE), ierr
   character(MAXNAMELEN)    :: temp_buff
   !----- formats -----
   character(*),parameter :: subname = "(cpl_comm_init) "
   character(*),parameter :: F00 = "('(cpl_comm_init) ',4a)"
   character(*),parameter :: F02 = "('(cpl_comm_init) ',a,6i4)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "setting up communicators, name = ",trim(name)

! Setup MPI
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

! Allocate Memory now that we know how many processors (&names) to expect
ALLOCATE( names(numprocs, MAXNAMELEN) )

! Node 0 gathers names, the other nodes send their's in.
if myid .eq. then
  names(0) = name ! Node 0's name is first
  recieve: do n=0,numprocs
    MPI_Recv( temp_buff, MAXNAMELEN, MPI_CHAR, n+1, 7, MPI_COMM_WORLD, status, ierr )
    ! Just add them?
    ! Names correspond to MPI_COMM_WORLD.rank
    names(n) = temp_buff
    ! Check for unique names?
    !check: do i=0,numprocs
    !  if names(i) .eq. temp_buff then
    !    cycle
    !  else
    !    names(n) = temp_buff
    !    exit
    !  end if
    !end do check
  end do recieve
  ! Broadcast list of names to all processors:
  MPI_Bcast( names, numprocs * MAXNAMELEN, MPI_CHAR, 0, MPI_COMM_WORLD, ierr )
! Non 0 node:
else
  MPI_Send( name, MAXNAMELEN, MPI_CHAR, 0, 7, MPI_COMM_WORLD, status, ierr )
  MPI_Bcast( names, numprocs * MAXNAMELEN, MPI_CHAR, 0, MPI_COMM_WORLD, ierr )
end if

!  Replacing the Call to MPH
!  1). Do a gather operation to get all the names
!  2). Create List of Unique Names
!  3). Broadcast names to every node
!  4). Loop over names and create communicators

!  Loop:
color = 0
commcreate: do n=0,numprocs
    if names[i] .eq. myname then
        MPI_COMM_SPLIT( MPI_COMM_WORLD, color, myid, comm, ierr )
        exit
    else
        color += 1
    end if
end do commcreate

! CLEAN UP
deallocate( names, unique_names )

end subroutine cpl_comm_init_mike(name,comm)

end module cpl_comm_mod

