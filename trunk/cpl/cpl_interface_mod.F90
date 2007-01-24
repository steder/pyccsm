!===============================================================================
! CVS: $Id: cpl_interface_mod.F90,v 1.1.1.1 2005/02/03 22:29:02 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/cpl/cpl_interface_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_interface_mod -- cpl6 API module for use in component models.
!
! !DESCRIPTION:
!     This module represents a major subsystem of cpl6.
!     This module contains cpl6 wrappers to insert in component model code.  
!     These wrappers present component models with an API using only standard
!     Fortran 90 datatypes, thus hiding the use of derived data types and lower
!     level libraries.  These wrappers are contain/wrap all the code necessary
!     for a component model to connect to, and exchange data with, the CCSM 
!     Coupler version 6.
!
! !REMARKS:
!     Component models communicate to the coupler via this module.  The coupler
!     communicates to the component models via it's cpl_msg_mod module.
!     Each routine in this module as a corresponding routine in cpl_msg_mod
!     -- these modules, and the order in which they are invoked, 
!     must be carefully coordinated. 

! !REVISION HISTORY:
!     2003-Jan-15 - T. Craig - change ibuf to infobuf datatype module
!     2002-Dec-05 - T. Craig - changed call from cpl_coupling to cpl_contract
!     2002-Sep-10 - T. Craig - abstracted functionality into cpl_coupling_mod.F90
!     2001-Aug-16 - B. Kauffman - reorganized code according to arch document.
!     2001-Mar-20 - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_interface_mod

! !USES:

   use cpl_mct_mod         ! mct interface
   use cpl_comm_mod        ! mpi/mph communicator info
   use cpl_fields_mod      ! coupler/model data field indicies
   use cpl_bundle_mod      ! defines bundle
   use cpl_domain_mod      ! defines domain
   use cpl_infobuf_mod     ! defines infobuf
   use cpl_contract_mod    ! defines contract
   use cpl_kind_mod        ! defines cpl kinds
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use shr_sys_mod         ! share system routines
   use shr_timer_mod       ! share timer routines
   use shr_mpi_mod         ! mpi layer

   implicit none

   private   ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_interface_init
   public :: cpl_interface_finalize
   public :: cpl_interface_ibufSend
   public :: cpl_interface_ibufRecv
   public :: cpl_interface_infobufSend
   public :: cpl_interface_infobufRecv
   public :: cpl_interface_contractSend
   public :: cpl_interface_contractRecv
   public :: cpl_interface_contractInit
   public :: cpl_interface_dbugSet   ! set this module's internal dbug level

   interface cpl_interface_ibufSend; module procedure cpl_interface_infobufSend; end interface
   interface cpl_interface_ibufRecv; module procedure cpl_interface_infobufRecv; end interface

! !PUBLIC DATA MEMBERS:

  ! none

!EOP


   !--- module variables ---
   character(*),parameter :: modName = "cpl_interface_mod"

   integer(IN),save  :: timer01,timer02,timer03,timer04    ! timers
   integer(IN),save  :: timer11,timer12,timer13,timer14    ! timers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_init -- initialize the coupling/mpi environment.
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
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_init(name,comm)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: name ! name of component name
   integer(IN) ,intent(out) :: comm ! communicator group for component

!EOP

   integer(IN)          :: n    ! generic loop index

   character(*),parameter :: subName = '(cpl_interface_init) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_comm_init(name, comm)

end subroutine cpl_interface_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_finalize -- terminate the coupling/mpi environment.
!
! !DESCRIPTION:
!    Calls mpi\_finalize().
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2001-mmm-dd -
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_finalize(cname)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in) :: cname  ! component name
   integer(IN)             :: rcode  ! return code

!EOP

   character(*),parameter :: subName = '(cpl_interface_finalize) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   call shr_mpi_finalize(subName//" MPI finalize")

end subroutine cpl_interface_finalize

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractInit -- initialize contract
!
! !DESCRIPTION:
!     Initialize a contract, but NOT the router.
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-Jul-30 - T.Craig -- prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_contractInit(contract,my_name,other_name,fields,ibufi,buf,ibufr,bunname,decomp)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_contract)   ,intent(out)   :: contract   ! contract
   character(*)         ,intent(in)    :: my_name    ! my component name
   character(*)         ,intent(in)    :: other_name ! other component name
   character(*)         ,intent(in)    :: fields     ! fields char string for bun
   integer(IN) ,optional,intent(inout) :: ibufi(:)   ! info buffer ints
   real(R8)    ,optional,intent(in)    :: buf(:,:)   ! data buffer
   real(R8)    ,optional,intent(inout) :: ibufr(:)   ! info buffer reals
   character(*),optional,intent(in)    :: bunname
   integer(IN) ,optional,intent(in)    :: decomp     ! decomposition type

!EOP

   !--- local ---
   character(CL)          :: bn  ! bundle name
   character(*),parameter :: subName = '(cpl_interface_contractInit) '
   integer(IN)            :: decomp_type

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   decomp_type = 1
   if (present(decomp)) then
     decomp_type = decomp
   endif

   if (present(ibufi)) then
     contract%infobuf%ibuf = ibufi
   endif

   if (present(ibufr)) then
     contract%infobuf%rbuf = ibufr
   endif

   if (present(buf)) then
      call cpl_contract_Init('send',contract,my_name,other_name,buf,decomp=decomp_type)
   else
      call cpl_contract_Init('recv',contract,my_name,other_name,decomp=decomp_type)
   endif

   if (present(bunname)) then
     bn = bunname
   else
     bn = "undef"
   endif
   call cpl_bundle_init(contract%bundle,bn,fields,contract%domain)

end subroutine cpl_interface_contractInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_infobufSend -- send initial data/msg to coupler.
!
! !DESCRIPTION:
!     Send time-invariant data such as a domain description.
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2001-Aug-16 -
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_infobufSend(cname,ibufi,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)       :: cname      ! component name
   integer(IN),optional,intent(inout)    :: ibufi(:)   ! info buffer ints
   real(R8)   ,optional,intent(inout)    :: ibufr(:)   ! info buffer reals

!EOP

   integer(IN)             :: pid          ! mpi process id
   integer(IN),parameter   :: tag=1002     ! mpi msg tag
   type(cpl_infobuf)       :: infobuf      ! local info buffer

   character(*),parameter :: subName = '(cpl_interface_infobufSend) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   !--- identify cpl_comm_wrld pid for other component's pe0 ---
   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   if (present(ibufi)) then
     infobuf%ibuf = ibufi
   endif

   if (present(ibufr)) then
     infobuf%rbuf = ibufr
   endif

   !--- send ibuf (tell cpl how big local & global component data is) ---
   if (cpl_comm_comp_pid == 0) then
     call cpl_infobuf_send(infobuf,pid,tag,cpl_comm_wrld)
   endif


end subroutine cpl_interface_infobufSend

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_infobufRecv -- receive initial data/msg from cpl.
!
! !DESCRIPTION:
!     Receive time-invariant data from coupler.
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2001-mmm-dd - 
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_infobufRecv(cname,ibufi,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)  :: cname  ! component name
   integer(IN),optional,intent(out) :: ibufi(cpl_infobuf_ibufSize) ! info-buffer
   real(R8)   ,optional,intent(out) :: ibufr(cpl_infobuf_rbufSize) ! info-buffer

!EOP

   integer(IN) :: pid                     ! mph process ID
   integer(IN),parameter :: pid0 = 0      ! component's root pid
   integer(IN),parameter :: tag = 1002    ! mpi msg tag
   type(cpl_infobuf)     :: infobuf       ! local infobuf

   character(*),parameter :: subName = '(cpl_interface_infobufRecv) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   !--- recv info-buffer data ---
   if (cpl_comm_comp_pid == pid0) then
     call cpl_infobuf_recv(infobuf,pid,tag,cpl_comm_wrld)
   endif
   call cpl_infobuf_bcast(infobuf,pid0,cpl_comm_comp)

   if (present(ibufi)) then
     ibufi = infobuf%ibuf
   endif

   if (present(ibufr)) then
     ibufr = infobuf%rbuf
   endif

end subroutine cpl_interface_infobufRecv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractSend -- send data/msg to coupler.
!
! !DESCRIPTION:
!     Send time-variant data such as state variables and forcing fields.
! 
! !REVISION HISTORY:
!    2001-mmm-dd - 
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_contractSend(cname,contract,ibufi,buf,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)    :: cname      ! component name
   type(cpl_contract)  ,intent(inout) :: contract   ! data buffer and domain
   integer(IN),optional,intent(inout) :: ibufi(:)   ! info buffer ints
   real(R8)   ,optional,intent(inout) :: ibufr(:)   ! info buffer reals
   real(R8)   ,optional,intent(inout) ::  buf(:,:)  ! data buffer

!EOP

   logical,save :: first_call = .true.    ! first time in subroutine
   integer(IN) :: AVsiz                   ! size of a field in the AttrVect
   integer(IN) :: AVnum                   ! number of fields in AttrVect
   integer(IN) :: i,j,n                   ! dummy variables
   integer(IN) :: pid                     ! mpi process ID
   integer(IN),parameter :: tag = 1003    ! msg tag

   character(*),parameter :: subName = '(cpl_interface_contractSend) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   !--- setup timers ---
   if (first_call) then
     first_call = .false.
     call shr_timer_get(timer01,'comp_send total')
     call shr_timer_get(timer02,'comp_send reorder')
     call shr_timer_get(timer03,'comp_send cpl_contract_Send')
     call shr_timer_get(timer04,'comp_send diagnostics')
   endif

   call shr_timer_start(timer01)

   !--- copy data into contract ---
   if (present(ibufi)) then
     contract%infobuf%ibuf = ibufi
   endif

   if (present(ibufr)) then
     contract%infobuf%rbuf = ibufr
   endif

   if (present(buf)) then
     AVnum = cpl_mct_aVect_nRAttr(contract%bundle%data)
     AVsiz = cpl_mct_aVect_lsize(contract%bundle%data)

     if (AVsiz /= size(buf,1) .or. AVnum /= size(buf,2)) then
       write(6,*) subName,' ERROR in buffer/contract buffer size:',  &
                  trim(contract%bundle%name)
       write(6,*) subName,' sending buffer, size:',size(buf,1),size(buf,2)
       write(6,*) subName,'contract buffer, size:',AVsiz,AVnum
       call shr_sys_flush(6)
     endif

     call shr_timer_start(timer02)
     !--- reorder bundle data as per mct data structure ---
#ifdef CPP_VECTOR
     do j=1,AVnum
     do i=1,AVsiz
#else
     do i=1,AVsiz
     do j=1,AVnum
#endif
       contract%bundle%data%rattr(j,i)=buf(i,j)
     enddo
     enddo
     call shr_timer_stop(timer02)
   endif

   !--- diagnostics ---
   if ( dbug >= 2) then
     call shr_timer_start(timer04)
     call cpl_bundle_gsum(contract%bundle,contract%bundle%dom%lGrid,'aream', &
                                          contract%bundle%dom%lGrid,'mask', &
                          scalar=cpl_const_rearth2,istr='send '//trim(cname))
     call shr_timer_stop(timer04)
   endif

   call shr_timer_start(timer03)
   
   !--- identify pid for component pe0 ---
   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   !--- send data ---
   call cpl_contract_Send(contract,cpl_comm_comp_pid,cpl_comm_comp,pid)

   call shr_timer_stop(timer03)
   call shr_timer_stop(timer01)

end subroutine cpl_interface_contractSend

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractRecv -- receive data/msg from coupler.
!
! !DESCRIPTION:
!     Receive time-variant data such as state variables and forcing fields.
! 
! !REVISION HISTORY:
!    2001-mmm-dd - 
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_contractRecv(cname,contract,ibufi,buf,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)    :: cname      ! component name
   type(cpl_contract)  ,intent(inout) :: contract   ! data buffer and domain
   integer(IN),optional,intent(out)   :: ibufi(:)   ! info buffer ints
   real(R8)   ,optional,intent(out)   :: ibufr(:)   ! info buffer reals
   real(R8)   ,optional,intent(out)   ::  buf(:,:)  ! data buffer

!EOP

   logical,save:: first_call = .true.     ! first time in subroutine
   integer(IN) :: AVsiz                   ! size of a field in the AttrVect
   integer(IN) :: AVnum                   ! number of fields in AttrVect
   integer(IN) :: i,j,n                   ! dummy variables
   integer(IN) :: pid                     ! generic pid
   integer(IN),parameter :: pid0 = 0      ! component's root pid
   integer(IN),parameter :: tag = 1003    ! mpi msg tag

   character(*),parameter :: subName = '(cpl_interface_contractRecv) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   if (first_call) then
     call shr_timer_get(timer11,'comp_recv total')
     call shr_timer_get(timer12,'comp_recv reorder')
     call shr_timer_get(timer13,'comp_recv cpl_contract_Recv')
     call shr_timer_get(timer14,'comp_recv diagnostics')
     first_call = .false.
   endif

   if (present(buf)) then
     AVnum = cpl_mct_aVect_nRAttr(contract%bundle%data)
     AVsiz = cpl_mct_aVect_lsize(contract%bundle%data)
     if (AVsiz /= size(buf,1) .or. AVnum /= size(buf,2)) then
       write(6,*) subName,' ERROR in buffer/contract buffer size'
       write(6,*) subName,' sending buffer, size:',size(buf,1),size(buf,2)
       write(6,*) subName,'contract buffer, size:',AVsiz,AVnum
       call shr_sys_flush(6)
     endif
   endif

   call shr_timer_start(timer11)
   call shr_timer_start(timer13)

   !--- identify pid for component pe0 ---
   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   !--- receive data ---
   call cpl_contract_Recv(contract,cpl_comm_comp_pid,cpl_comm_comp,pid)
   call shr_timer_stop(timer13)

   !--- diagnostics ---
   if ( dbug >= 2) then
     call shr_timer_start(timer14)
     call cpl_bundle_gsum(contract%bundle,contract%bundle%dom%lGrid,'aream', &
                          scalar=cpl_const_rearth2,istr='recv '//trim(cname))
     call shr_timer_stop(timer14)
   endif

   !--- copy data out of contract ---
   if (present(ibufi)) then
     ibufi = contract%infobuf%ibuf
   endif

   if (present(ibufr)) then
     ibufr = contract%infobuf%rbuf
   endif

   if (present(buf)) then
     call shr_timer_start(timer12)
     !--- reorder bundle data as per mct data structure ---
#ifdef CPP_VECTOR
     do j=1,AVnum
     do i=1,AVsiz
#else
     do i=1,AVsiz
     do j=1,AVnum
#endif
       buf(i,j)=contract%bundle%data%rattr(j,i)
     enddo
     enddo
     call shr_timer_stop(timer12)
   endif

   call shr_timer_stop(timer11)

end subroutine cpl_interface_contractRecv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_dbugSet  -- set this module's internal debug level.
!
! !DESCRIPTION:
!    Set this module's internal debug level: 0,1,2,3 (lowest to highest). 
!
! !REVISION HISTORY:
!    2003-Jan-21 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_dbugSet(level)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in) :: level  ! requested debug level

!EOP

   !----- local -----
   integer(IN):: newLevel  ! new debug level

   !----- formats -----
   character(*),parameter :: F00 = "('(cpl_interface_dbugSet) ',a,i1,a)"

!-------------------------------------------------------------------------------
!  Set module's internal debug level: 0,1,2, or 3 (lowest to highest)
!-------------------------------------------------------------------------------

   newLevel = max(0,min(3,level))

   !--- correct invalid level values ---
   if (newLevel /= level) then
      write(6,F00) 'WARNING: level ',level,' not in {0,1,2,3} '
      write(6,F00) 'WARNING: resetting level to ',newLevel
   end if

   if (dbug>0 .OR. dbug/=newLevel) write(6,F00) 'set debug level to ',newLevel

   dbug = newLevel

end subroutine cpl_interface_dbugSet

!===============================================================================
!===============================================================================

end module cpl_interface_mod
