!BOP ===========================================================================
!
! !MODULE: data_mod -- data declaration
!
! !DESCRIPTION:
!  Does data that's shared across the init, run, and finalize methods
!
! !REVISION HISTORY:
!     2003-Apr-20 - T. Craig - subroutinize dead model, add data module
!
! !INTERFACE: ------------------------------------------------------------------

MODULE data_mod

! !USES:

   use shr_kind_mod
   use cpl_fields_mod
   use cpl_contract_mod

   implicit none

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   ! no public functions

! !PUBLIC DATA MEMBERS:

   integer(SHR_KIND_IN),parameter :: max_con = 2
   integer(SHR_KIND_IN)           :: n_con 

   type(cpl_contract):: contractS(max_con)           ! send contract
   type(cpl_contract):: contractR(max_con)           ! recv contract
   integer(SHR_KIND_IN)           :: nx(max_con)       ! local i dimension
   integer(SHR_KIND_IN)           :: ny(max_con)       ! local j dimension

   real(SHR_KIND_R8), allocatable :: fields_m2c(:,:,:) ! component fields buffer
   real(SHR_KIND_R8), allocatable :: fields_c2m(:,:,:) ! component fields buffer
   real(SHR_KIND_R8), allocatable :: gbuf(:,:)         ! cpl buffer for grid info
   real(SHR_KIND_R8), allocatable :: sbuf(:,:)         ! buffer for send contract
   real(SHR_KIND_R8), allocatable :: rbuf(:,:)         ! buffer for recv contract
   integer(SHR_KIND_IN)           :: isbuf(cpl_fields_ibuf_total)
                                                       ! info-buffer to send
   integer(SHR_KIND_IN)           :: irbuf(cpl_fields_ibuf_total)
                                                       ! info-buffer recieved

   integer(SHR_KIND_IN)           :: fields_m2c_total  ! # of flds m2c
   integer(SHR_KIND_IN)           :: fields_c2m_total  ! # of flds c2m
   character(SHR_KIND_CL)         :: fields_m2c_list   ! list of flds m2c  
   character(SHR_KIND_CL)         :: fields_c2m_list   ! list of flds c2m

   !--- runoff domain/contract info ---
   integer(SHR_KIND_IN)           :: fields_r2c_total
   character(SHR_KIND_CL)         :: fields_r2c_list 
   real(SHR_KIND_R8), allocatable :: fields_r2c(:,:,:)

   integer(SHR_KIND_IN)           :: isbuf_r(cpl_fields_ibuf_total)
   real(SHR_KIND_R8), allocatable :: gbuf_r(:,:)
   real(SHR_KIND_R8), allocatable :: sbuf_r(:,:)

   integer(SHR_KIND_IN)           :: irbuf_g(cpl_fields_ibuf_total)
   real(SHR_KIND_R8), allocatable :: rbuf_g(:,:)

   !--- timers ---
   integer(SHR_KIND_IN)   :: timer01,timer02,timer03,timer04,timer05

   !--- stdin input stuff ---
   character(SHR_KIND_CS) :: str               ! cpp  defined model name
   character(SHR_KIND_CS) :: myModelName       ! user defined model name
   integer(SHR_KIND_IN)   :: nxg    (max_con)  ! global dim i-direction
   integer(SHR_KIND_IN)   :: nyg    (max_con)  ! global dim j-direction
   integer(SHR_KIND_IN)   :: nproc_x(max_con)  ! num of i pes (type 3)
   integer(SHR_KIND_IN)   :: seg_len(max_con)  ! length of segs (type 4)
   integer(SHR_KIND_IN)   :: ncpl   (max_con)  ! cpling freq (# per day)
   integer(SHR_KIND_IN)   :: ncpl_v1(max_con)  ! relative coupling freq 
   integer(SHR_KIND_IN)   :: decomp_type(max_con) ! data decomp type:
                                               !  1 = 1d decomp by lat, eq chunks 
                                               !  2 = 1d decomp by lon, eq chunks
                                               !  3 = 2d decomp, eq chunks
                                               !  4 = segmented decomp, repetitive
   real(SHR_KIND_R8)      ::   simTime         !  simulation time proxy (secs)   

   !--- other ---
   integer(SHR_KIND_IN)   :: ncomp             ! component index
   integer(SHR_KIND_IN)   :: dbug = 0          ! debug level (higher is more)

!EOP

END MODULE data_mod
