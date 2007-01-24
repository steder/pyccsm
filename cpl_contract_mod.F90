!===============================================================================
! CVS: $Id$
! CVS: $Source$
! CVS: $Name$
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_contract_mod -- coupler/component contract type
!
! !DESCRIPTION:
!     Contains the contract type.
!     This module represents a major subsystem of cpl6.
!     This is the module/package that contains fundamental, low-level 
!     coupler/component data exchange routines.  It's functionality
!     relies heavily on the lower level MCT, MPH, and MPI libraries.
!
! !REMARKS:
!     The coupler and components use this module to carry out fundamental
!     communication.  The coupler calls this mainly from cpl_msg and the
!     components call this from cpl_interface.
!
! !REVISION HISTORY:
!     2002-Jul-16 - T. Craig - abstracted basic functionality from
!                   cpl_msg and cpl_interface to this layer.
!     2002 Aug 01 - T. Craig - prototype for contract datatype
!     2002 Dec 05 - T. Craig - combined cpl-coupling module and cpl_contract
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_contract_mod

! !USES:

   use shr_timer_mod       ! timers
   use cpl_kind_mod        ! kinds
   use cpl_mct_mod         ! mct interface
   use cpl_comm_mod        ! mpi/mph communicator info
   use cpl_fields_mod      ! fields module
   use cpl_bundle_mod      ! defines bundle
   use cpl_domain_mod      ! defines domain
   use cpl_infobuf_mod     ! defines infobuf 
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug

   implicit none

   private ! except

! !PUBLIC TYPES:

   public :: cpl_contract

   type cpl_contract
      type(cpl_infobuf)     :: infobuf       ! infobuf that goes with contract
      type(cpl_bundle)      :: bundle        ! bundle
      type(cpl_domain)      :: domain        ! domain info (grid with decomp)
      type(cpl_mct_Router)  :: rtr           ! memory for new router if needed
   end type cpl_contract

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_contract_execute
   public :: cpl_contract_send
   public :: cpl_contract_recv
   public :: cpl_contract_init
   public :: cpl_contract_initSend
   public :: cpl_contract_initRecv

! !PUBLIC DATA MEMBERS:

  ! none

!EOP

   character(*),parameter :: modName = "cpl_contract_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_contract_execute -- send/recv data/msg to component.
!
! !DESCRIPTION:
!    Send/Recv msg containing integer parameters
!    Send/Recv msg containing time-variant data such as state and forcing fields.
!
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-09-10 - T.Craig - merged cpl_contract_send, cpl_contract_recv
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_contract_execute(srtype,contract,mypid,comm,otherpid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)      ,  intent(in)   :: srtype   ! 'send' or 'recv'
   type(cpl_contract),  intent(inout):: contract ! contract
   integer(IN),         intent(in)   :: mypid    ! my mpi process ID
   integer(IN),         intent(in)   :: comm     ! local communicator group
   integer(IN),         intent(in)   :: otherpid ! mpi process ID to send to

!EOP

   !--- local ---
   integer(IN),parameter :: tag=1003      ! generic mpi tag
   integer(IN),parameter :: pid0=0        ! root pid

   logical,save :: first_call = .true.    ! first time in subroutine
   integer(IN),save  :: timer00                            ! timers
   integer(IN),save  :: timer01,timer02,timer03,timer04    ! timers
   integer(IN),save  :: timer11,timer12,timer13,timer14    ! timers

   !--- formats ---
   character(*),parameter :: subName = "(cpl_contract_execute) "
   character(*),parameter :: F00 = "('(cpl_contract_execute) ',8a)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   !--- setup timers ---
   if (first_call) then
     first_call = .false.
     call shr_timer_get(timer00,'cpl_contract_execute timer00, total')
     call shr_timer_get(timer01,'cpl_contract_execute timer01, ibuf send')
     call shr_timer_get(timer02,'cpl_contract_execute timer02, mct send')
     call shr_timer_get(timer03,'cpl_contract_execute timer03')
     call shr_timer_get(timer04,'cpl_contract_execute timer04')
     call shr_timer_get(timer11,'cpl_contract_execute timer11, zero')
     call shr_timer_get(timer12,'cpl_contract_execute timer12, ibuf recv')
     call shr_timer_get(timer13,'cpl_contract_execute timer13, ibuf bcast')
     call shr_timer_get(timer14,'cpl_contract_execute timer14, mct recv')
   endif

   call shr_timer_start(timer00)
   if (srtype == 'send') then
     !--- send info-buffer data ---
     call shr_timer_start(timer01)
     if (mypid == pid0) then
       call cpl_infobuf_send(contract%infobuf,otherpid,tag,cpl_comm_wrld)
     call shr_timer_stop(timer01)
     endif

     !--- send bundle data ---
     call shr_timer_start(timer02)
     call cpl_mct_send(contract%bundle%data,contract%rtr)
     call shr_timer_stop(timer02)
   endif
   if (srtype == 'recv') then
     call shr_timer_start(timer11)
     call cpl_bundle_zero(contract%bundle)
     call shr_timer_stop(timer11)
     !--- recv info-buffer data ---
     call shr_timer_start(timer12)
     if (mypid == pid0) then
       call cpl_infobuf_recv(contract%infobuf,otherpid,tag,cpl_comm_wrld)
     endif
     call shr_timer_stop(timer12)
     call shr_timer_start(timer13)
     call cpl_infobuf_bcast(contract%infobuf,pid0,comm)
     call shr_timer_stop(timer13)

     !--- recv bundle data ---
     call shr_timer_start(timer14)
     call cpl_mct_recv(contract%bundle%data,contract%rtr)
     contract%bundle%cnt = 1
     call shr_timer_stop(timer14)
   endif
   call shr_timer_stop(timer00)

end subroutine cpl_contract_execute

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_contract_send -- send data/msg to component.
!
! !DESCRIPTION:
!    Send msg containing integer parameters
!    Send msg containing time-variant data such as state and forcing fields.
!    Calls cpl\_coupling\_contract
!
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-08-01 - T.Craig - abstracted from cpl_interface and cpl_msg
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_contract_send(contract,mypid,comm,otherpid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_contract),  intent(inout):: contract ! contract
   integer(IN),         intent(in)   :: mypid    ! my mpi process ID
   integer(IN),         intent(in)   :: comm     ! local communicator group
   integer(IN),         intent(in)   :: otherpid ! mpi process ID to send to

!EOP

   !--- formats ---
   character(*),parameter :: subName = "(cpl_contract_send) "
   character(*),parameter :: F00 = "('(cpl_contract_send) ',8a)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   call cpl_contract_execute('send',contract,mypid,comm,otherpid)
                                         
end subroutine cpl_contract_send

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_contract_recv -- receive data/msg from component.
!
! !DESCRIPTION:
!    Receive msg containing integer parameters
!    Receive msg containing time-variant data such as state and forcing fields.
!    Calls cpl\_coupling\_contract
!
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-08-01 - T.Craig - abstracted from cpl_interface and cpl_msg
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_contract_recv(contract,mypid,comm,otherpid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_contract),  intent(inout):: contract  ! contract
   integer(IN),         intent(in)   :: mypid     ! my mpi process ID
   integer(IN),         intent(in)   :: comm      ! local communicator group
   integer(IN),         intent(in)   :: otherpid  ! mpi process ID to recv from

!EOP

   !--- formats ---
   character(*),parameter :: subName = "(cpl_contract_recv) "
   character(*),parameter :: F00 = "('(cpl_contract_recv) ',8a)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

    call cpl_contract_execute('recv',contract,mypid,comm,otherpid)

end subroutine cpl_contract_recv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_contract_init -- initialize contract
!
! !DESCRIPTION:
!     Initialize a contract
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-Jul-30 - T.Craig -- prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_contract_init(srtype,contract,my_name,other_name,buf,decomp)

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in) :: srtype        ! 'send' or 'recv'
   character(*)        ,intent(in) :: my_name       ! component name (me)
   character(*)        ,intent(in) :: other_name    ! component name (other)
   type(cpl_contract)  ,intent(out):: contract      ! contract
   real(R8),optional   ,intent(in) ::  buf(:,:)     ! data buffer
   integer(IN),optional,intent(in) :: decomp        ! recv side decomp type
                                                    ! 1 = 1d in lat
                                                    ! 2 = 1d in lon

!EOP

   !--- local ---
   integer(IN)             :: i,j,k,n      ! generic indicies
   integer(IN)             :: decomp_type  ! local decomposition type value
   integer(IN)             :: nseg         ! counts number of segments for gsMap
   integer(IN),parameter   :: pid0=0       ! mpi process id for root pe
   integer(IN)             :: ierr         ! routine return code
   integer(IN),pointer     :: indx(:)     ! used to init gsMap 
   integer(IN),allocatable :: start(:)     ! used to init gsMap 
   integer(IN),allocatable :: count(:)     ! used to init gsMap 
   integer(IN)             :: lSize        ! local  grid size
   integer(IN)             :: gSize        ! global grid size
   integer(IN)             :: giSize       ! global grid size in x
   integer(IN)             :: gjSize       ! global grid size in y
   integer(IN)             :: cid_me       ! mph component id (mine)
   integer(IN)             :: cid_other    ! mph component id (other component)
   integer(IN)             :: pid_me       ! root processor id in comm_world (mine)
   integer(IN)             :: pid_other    ! root processor id in comm_world (other component)
   integer(IN)             :: comm_me      ! mpi communicator group (mine)
   character*1             :: suffix       ! suffix
   integer(IN),parameter   :: tag=1001     ! mpi msg tag

   !--- formats ---
   character(*),parameter :: subName = "(cpl_contract_init) "
   character(*),parameter :: F00 = "('(cpl_contract_init) ',8a)"
   character(*),parameter :: F01 = "('(cpl_contract_init) ',2(a,i6),2(a,2i6))"
   character(*),parameter :: F02 = "('(cpl_contract_init) ',2(a,2i6))"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) trim(my_name),'-',trim(srtype),'-',trim(other_name)

   decomp_type = 1
   if (present(decomp)) decomp_type = decomp

   !----------------------------------------------------------------------------
   !  initialize communicator information
   !----------------------------------------------------------------------------
   comm_me = cpl_comm_comp
   cid_me  = cpl_comm_mph_cid
   pid_me  = cpl_comm_wrld_pe0

   if      (other_name == cpl_fields_atmname) then
      cid_other = cpl_comm_mph_cid_atm
      pid_other = cpl_comm_wrld_pe0_atm
      suffix = "a "
   else if (other_name == cpl_fields_lndname) then
      cid_other = cpl_comm_mph_cid_lnd
      pid_other = cpl_comm_wrld_pe0_lnd
      suffix = "l "
   else if (other_name == cpl_fields_rtmname) then
      cid_other = cpl_comm_mph_cid_lnd
      pid_other = cpl_comm_wrld_pe0_lnd
      suffix = "r "
   else if (other_name == cpl_fields_ocnname) then
      cid_other = cpl_comm_mph_cid_ocn
      pid_other = cpl_comm_wrld_pe0_ocn
      suffix = "o "
   else if (other_name == cpl_fields_icename) then
      cid_other = cpl_comm_mph_cid_ice
      pid_other = cpl_comm_wrld_pe0_ice
      suffix = "i "
   else if (other_name == cpl_fields_cplname) then
      cid_other = cpl_comm_mph_cid_cpl
      pid_other = cpl_comm_wrld_pe0_cpl
      if      (my_name == cpl_fields_atmname) then
         suffix = "a "
      else if (my_name == cpl_fields_icename) then
         suffix = "i "
      else if (my_name == cpl_fields_lndname) then
         suffix = "l "
      else if (my_name == cpl_fields_ocnname) then
         suffix = "o "
      else
         write(6,F00) 'ERROR: this should never happen'
         write(6,F00) 'unrecognized my_name = ',my_name
      endif
   else
      write(6,F00) 'ERROR: this should never happen'
      write(6,F00) 'unrecognized other_name = ',other_name
   endif
   write(6,F00) 'OTHER IDS',other_name
   write(6,*) subName,'OTHER IDS',cid_other,pid_other

   !----------------------------------------------------------------------------
   ! send/recv infobuf 
   !----------------------------------------------------------------------------
   if (srtype == 'send') then
     if (cpl_comm_comp_pid == 0) then
        call cpl_infobuf_send(contract%infobuf,pid_other,tag,cpl_comm_wrld)
     endif
     write(6,F00) 'finished infobuf send'
   endif
   if (srtype == 'recv') then
     if (cpl_comm_comp_pid == 0) then
        call cpl_infobuf_recv(contract%infobuf,pid_other,tag,cpl_comm_wrld)
     endif
     write(6,F00) 'finished infobuf recv'
     call cpl_infobuf_bcast(contract%infobuf,pid0,comm_me)
     write(6,F00) 'finished infobuf bcast'
   endif
   write(6,F00) 'finished infobuf transfer'

   !----------------------------------------------------------------------------
   ! get/set local index values
   !----------------------------------------------------------------------------
   gSize = contract%infobuf%ibuf(cpl_fields_ibuf_gsize)
   giSize= contract%infobuf%ibuf(cpl_fields_ibuf_gisize)
   gjSize= contract%infobuf%ibuf(cpl_fields_ibuf_gjsize)
   if (srtype == 'send') then
     lSize = contract%infobuf%ibuf(cpl_fields_ibuf_lsize)
     allocate(indx(lSize))
     indx(:)=buf(:,cpl_fields_grid_index)
   endif

   if (srtype == 'recv') then
#if (OLD_STUFF)
     ! assume 1d decomposition of the recv side, it doesn't know the size of
     ! the grid until it recieves the infobuf.  compute lSize.
     lSize = gSize/cpl_comm_comp_npe
     if (cpl_comm_comp_pid < gsize-lSize*cpl_comm_comp_npe) then
       lSize = lSize + 1
       write(6,*) subName,' adjusting lSize ',cpl_comm_comp_pid,gSize,cpl_comm_comp_npe,lSize
     endif
     allocate(indx(lSize))
     do n=1,lSize
        indx(n) = n + cpl_comm_comp_pid*lSize 
     end do
#else
   write(6,F00) 'cpl_contract_decomp call'
   write(6,*) subName,decomp_type,giSize,gjSize,cpl_comm_comp_pid, &
             cpl_comm_comp_npe
     call cpl_contract_decomp(decomp_type,giSize,gjSize,cpl_comm_comp_pid,cpl_comm_comp_npe,lsize,indx)
   write(6,F00) 'finished cpl_contract_decomp call'
   write(6,*) subName,lsize
#endif
   endif

   !----------------------------------------------------------------------------
   ! initialize contract' gsMap (based on index data) ---
   !----------------------------------------------------------------------------
   if (lSize == 0) then
     !--- this process has a null/empty segment ---
      nseg = 0
      allocate(start(nseg),count(nseg))
      if (dbug>1) then
         write(6,F01) "gsMap_init, nSeg =",nSeg,", gSize =",gsize
         call shr_sys_flush(6)
      end if
   else
      !--- compute segment's start indicies and length counts ---
     nseg=1
     do n=2,lSize
        i = indx(n-1)
        j = indx(n)
        if ( j-i /= 1) nseg=nseg+1
     enddo

     allocate(start(nseg),count(nseg))

     nseg=1
     start(nseg)=indx(1)
     count(nseg)=1
     do n=2,lSize
        i = indx(n-1)
        j = indx(n)
        if ( j-i /= 1) then
            nseg=nseg+1
            start(nseg)=indx(n)
            count(nseg)=1
         else
            count(nseg)=count(nseg)+1
         endif
     enddo
      if (dbug>1) then
         write(6,F01) "gsMap_init, nSeg =",nSeg,", gSize =",gsize
         write(6,F02) "gsMap start(1),count(1)="      ,start(1)   ,count(1),  &
                          ", start(nSeg),count(nSeg)=",start(nSeg),count(nSeg)
         call shr_sys_flush(6)
      endif
   endif

   call cpl_mct_gsMap_init(contract%domain%gsMap,start,count,pid0,comm_me,cid_me,gsize=gsize)

   deallocate(start,count)
   deallocate(indx)

   !----------------------------------------------------------------------------
   ! init router
   !----------------------------------------------------------------------------
   write(6,*) "initializing router to",cid_other,"!"
   call cpl_mct_router_init(cid_other,contract%domain%gsMap,comm_me,contract%rtr)
   write(6,*) "done initializing router to",cid_other,"!"

   !----------------------------------------------------------------------------
   ! init contracts, setup lGrid
   !   on the send side lGrid exists
   !   on the recv side, need to recv lGrid
   !----------------------------------------------------------------------------
   if (srtype == 'send') contract%domain%name = trim(my_name)//' contract domain'
   if (srtype == 'recv') contract%domain%name = trim(other_name)//' contract domain'
   contract%domain%suffix = trim(suffix)
   contract%domain%n  = contract%infobuf%ibuf(cpl_fields_ibuf_gSize)
   contract%domain%ni = contract%infobuf%ibuf(cpl_fields_ibuf_giSize)
   contract%domain%nj = contract%infobuf%ibuf(cpl_fields_ibuf_gjSize)

   call cpl_mct_aVect_init(contract%domain%lGrid," ",cpl_fields_grid_fields,lSize)

   contract%domain%lGrid%rAttr(:,:) = -9999.0       ! special value

   if (srtype == 'send') then
      call cpl_mct_aVect_putRAttr(contract%domain%lGrid,"lon"  ,buf(:,cpl_fields_grid_lon  ),ierr)
      call cpl_mct_aVect_putRAttr(contract%domain%lGrid,"lat"  ,buf(:,cpl_fields_grid_lat  ),ierr)
      call cpl_mct_aVect_putRAttr(contract%domain%lGrid,"mask" ,buf(:,cpl_fields_grid_mask ),ierr)
      call cpl_mct_aVect_putRAttr(contract%domain%lGrid,"area" ,buf(:,cpl_fields_grid_area),ierr)
      call cpl_mct_aVect_putRAttr(contract%domain%lGrid,"index",buf(:,cpl_fields_grid_index),ierr)
!     call cpl_mct_aVect_putRAttr(contract%domain%lGrid,"mpad" ,buf(:,cpl_fields_grid_mpad ),ierr)
!     call cpl_mct_aVect_putRAttr(contract%domain%lGrid,"mcps" ,buf(:,cpl_fields_grid_mcps ),ierr)
      k = cpl_mct_aVect_indexRA(contract%domain%lGrid,"pid"  ,perrWith=subName)
      contract%domain%lGrid%rAttr(k,:) = cpl_comm_wrld_pid
   endif

   !----------------------------------------------------------------------------
   ! send/recv lGrid (requires a router)
   !----------------------------------------------------------------------------
   if (srtype == 'send') then
      call cpl_mct_send(contract%domain%lGrid,contract%rtr)
   endif
   if (srtype == 'recv') then
      call cpl_mct_recv(contract%domain%lGrid,contract%rtr)
   endif

   !--- write some dbug/sanity-check info to stdout ---
   call cpl_domain_info(contract%domain)

end subroutine cpl_contract_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_contract_initSend -- initialize contract, send side
!
! !DESCRIPTION:
!     Initialize a contract on the send side
!     Calls cpl\_coupling\_contract\_init
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-Sep-10 - T.Craig -- prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_contract_initSend(contract,my_name,other_name,buf)

! !INPUT/OUTPUT PARAMETERS:

   character(*)      ,intent(in)  :: my_name       ! component name (me)
   character(*)      ,intent(in)  :: other_name    ! component name (other)
   type(cpl_contract),intent(out) :: contract      ! contract
   real(R8)          ,intent(in)  ::  buf(:,:)     ! data buffer

!EOP

   !--- formats ---
   character(*),parameter :: subName = "(cpl_contract_initSend) "
   character(*),parameter :: F00 = "('(cpl_contract_initSend) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_contract_init('send',contract,my_name,other_name,buf)

end subroutine cpl_contract_initSend

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_contract_initRecv -- initialize contract, recv side
!
! !DESCRIPTION:
!     Initialize a contract on the recv side
!     Calls cpl\_coupling\_contract\_init
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-Sep-10 - T.Craig -- prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_contract_initRecv(contract,my_name,other_name)

! !INPUT/OUTPUT PARAMETERS:

   character(*)      ,intent(in) :: my_name       ! component name (me)
   character(*)      ,intent(in) :: other_name    ! component name (other)
   type(cpl_contract),intent(out):: contract      ! contract

!EOP

   !----- formats -----
   character(*),parameter :: subName = "(cpl_contract_initRecv) "
   character(*),parameter :: F00 = "('(cpl_contract_initRecv) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_contract_init('recv',contract,my_name,other_name)

end subroutine cpl_contract_initRecv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_contract_decomp -- create a decomp based on global info
!
! !DESCRIPTION:
!     Create a decomposition given a global grid and a decomp type.
!     This is used for contract initialization on the recv side.
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2003-Jun-13 - T.Craig -- first implemention
!
! !INTERFACE: ------------------------------------------------------------------

! Helper functions for the QuickSort algorithm:
recursive subroutine quicksort( A, p, r )
  integer(IN), pointer :: A(:) ! Array to sort
  integer(IN), intent(in) :: p, r ! lower , upper bounds of A
  integer(IN),  q ! pivot point of the sort
  ! This "if" causes the recursion to stop when the
  ! subarray is size 1.
  if ( p < r ) then
     call qs_partition( A, p, r, q )
     ! sort the two subarrays
     call quicksort( A, p, q - 1 )
     call quicksort( A, q + 1, r )
  end if
end subroutine quicksort

subroutine qs_partition( A, p, r, q )
  integer(IN), pointer :: A(:) ! Array to sort
  integer(IN), intent(in) :: p, r ! lower , upper bounds of A
  integer(IN), intent(out) :: q ! return value (midpoint of partially sorted A)
  integer(IN)  x ! pivot point of the sort
  integer(IN)  j, tmp ! loop index, tmp value for swap
  integer(IN)  i ! midpoint of the sorted A
  ! partition A into two subarrays, one greater then x
  ! the other less then x.
  x = A(r)
  i = p - 1
  do j = p, r
      if (A(j) <= x) then
        ! swap values:
        i = i + 1
        tmp = A(i)
        A(i) = A(j)
        A(j) = tmp
     endif
  end do
  ! Place the pivot point in the middle
  ! A(i+1) is sorted, the subarrays above/below it are not.
  tmp = A(r)
  A(r) = A(i+1)
  A(i+1) = tmp
  ! Return midpoint
  q = i+1
end subroutine qs_partition


subroutine cpl_contract_decomp(decomp,gi,gj,myid,npes,lsize,indx)

! !INPUT/OUTPUT PARAMETERS:

   integer(IN)      ,intent(in) :: decomp            ! decomp type
   integer(IN)      ,intent(in) :: gi,gj             ! global i and j sizes
   integer(IN)      ,intent(in) :: myid              ! my pe number
   integer(IN)      ,intent(in) :: npes              ! total number of pes
   integer(IN)      ,intent(out):: lsize             ! local size of decomp
   integer(IN),pointer          :: indx(:)           ! global index of decomp

!EOP

   !--- local ---
   integer(IN)             :: gsize,ier
   integer(IN)             :: lsize0
   integer(IN)             :: n,i,j,nl,nr,n1,itmp1,itmp2
   integer(IN),allocatable :: tmparr(:)

   !----- formats -----
   character(*),parameter :: subName = "(cpl_contract_decomp) "
   character(*),parameter :: F00 = "('(cpl_contract_decomp) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   gsize = gi*gj

   allocate(tmparr(gsize),stat=ier)
   if(ier /=0 ) write(0,*)'OUT OF MEMORY?',ier
   write(6,*) subName,' allocated memory for tmparr',gsize

   if (decomp == 1) then
     write(6,*) subName,' decomp 1'
      ! assume 1d decomposition in i direction
      do n=1,gsize
        tmparr(n) = n
      enddo
   elseif (decomp == 2) then
      ! assume 1d decomposition in j direction
      do j=1,gj
      do i=1,gi
        nr = (j-1)*gi + i
        nl = (i-1)*gj + j
        tmparr(nl) = nr
      enddo
      enddo
   elseif (decomp == 901) then
      ! test decomp, do not use
      do n=1,gsize
        tmparr(n) = n
      enddo
   else
      write(6,*) subName,' decomp not available, stop ',decomp
      call shr_sys_abort(trim(subName)//' decomp error ')
   endif

   lSize0 = gSize/npes
   lSize = lSize0
   write(6,*) subName,' lsizes ',lSize0,lSize
   if (myid < mod(gsize,npes)) then
      lSize = lSize + 1
      write(6,*) subName,' adjusting lSize ',myid,gSize,npes,lSize
   endif
   allocate(indx(lSize),stat=ier)
   if(ier /=0 ) write(0,*)'OUT OF MEMORY?',ier

   do n=1,lSize
      indx(n) = tmparr(n + myid*lSize0 +  min(myid,mod(gsize,npes)))
   end do
   write(6,*) subName,' Set index'

   !--- sort for mapping performance
!  do n=1,lSize-1
!  itmp1 = n
!  do n1 = n+1,lSize
!    if (indx(n1) < indx(itmp1)) then
!      itmp1 = n1
!    endif
!  enddo
!  itmp2 = indx(n)
!  indx(n) = indx(itmp1)
!  indx(itmp1) = itmp2
!  enddo
   !write(6,*) subName,"Sorting for mapping performance:"
   !call quicksort(indx,1,lSize)
   !write(6,*) subName,'Done sorting'

   deallocate(tmparr)

end subroutine cpl_contract_decomp

!===============================================================================

end module cpl_contract_mod

