#include <misc.h>
#include <preproc.h>

module clm_csmMod

#if (defined COUP_CSM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_csmMod
!
! !DESCRIPTION:
! Set of routines that define communication between the
! land model and flux coupler. The order of sends/receives is:
! 1) receive orbital data from coupler
! 2) send control data (grids and masks) to coupler
!    land grid does not have valid data, runoff grid does
! 3) receive valid land grid from flux coupler
! 4) send compressed runoff information to flux coupler
! 5) start normal send/recv communication patterm
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clm_varpar
#if (defined SPMD)
  use spmdMod        , only : masterproc, mpicom
  use spmdGathScatMod, only : gather_data_to_master
#else
  use spmdMod        , only : masterproc
#endif
  use mpiinc
  use cpl_fields_mod
  use cpl_contract_mod
  use cpl_interface_mod
  use cpl_comm_mod
  use RunoffMod        , only : runoff
  use shr_sys_mod      , only : shr_sys_irtc, shr_sys_flush ! csm_share system utility routines
  use system_messages  , only : allocation_err              ! allocation error output
  use abortutils       , only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: csm_setup          ! Setup, mpi_init
  public :: csm_shutdown       ! Shutdown, mpi_finalize
  public :: csm_initialize     ! Initialize contracts, etc
  public :: csm_recvgrid       ! Receive grid and land mask
  public :: csm_dosndrcv       ! Logic for determining if send/recv
  public :: csm_recv           ! Receive data from flux coupler
  public :: csm_send           ! Send data to flux coupler
  public :: csm_sendalb        ! Send initial albedos, surface temp and snow data
  public :: csm_flxave         ! Flux averaging rougine
  public :: restart_coupler    ! Restart code
  public :: compat_check_spval ! Checks that data sent from the coupler is valid
  public :: csm_compat         ! Checks compatibility of messages send/received

! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!  03.01.15 T. Craig Update for cpl6
!  03.04.27 M. Vertenstein, added qref_2m to communication and
!           generalized global sums to include all fields
!
!EOP
!
  private
!
! PRIVATE MEMBER FUNCTIONS:
!
! PRIVATE TYPES:
!
  integer   :: ibuffr(cpl_fields_ibuf_total)  ! Integer buffer from cpl
  integer   :: ibuffs(cpl_fields_ibuf_total)  ! Integer buffer to   cpl
  real(r8)  :: rbuffr(cpl_fields_rbuf_total)  ! Real    buffer from cpl
  real(r8)  :: rbuffs(cpl_fields_rbuf_total)  ! Real    buffer to   cpl
  type(cpl_contract)    :: contractRg         ! Contract for grid recvs from cpl
  type(cpl_contract)    :: contractR          ! Contract for recvs from cpl
  type(cpl_contract)    :: contractS          ! Contract for sends to   cpl
  type(cpl_contract)    :: contractSr         ! Contract for runoff sends to cpl
  real(r8), allocatable :: Gbuf(:,:)          ! Temporary generic buffer
  real(r8), allocatable :: bufS(:,:)          ! Send buffer for land
  real(r8), allocatable :: bufR(:,:)          ! Recv buffer for land
  real(r8), allocatable :: bufSr(:,:)         ! Send buffer for runoff
  real(r8), pointer     :: bufSglob(:,:)      ! Send global sum buffer for land
  real(r8), pointer     :: bufRglob(:,:)      ! Recv global sum buffer for land
  real(r8), pointer     :: bufSloc(:,:)       ! Send local sum buffer for land
  real(r8), pointer     :: bufRloc(:,:)       ! Recv local sum buffer for land
  real(r8), pointer     :: fieldS(:)          ! Global sum send field
  real(r8), pointer     :: fieldR(:)          ! Global sum receive field

  integer :: csm_nptg                         ! Loc sizes, grid coupling buffers
  integer :: csm_nptr                         ! Loc sizes, roff coupling buffers
  integer :: begp, endp                       ! per-proc beginning and ending pft indices
  integer :: begc, endc                       ! per-proc beginning and ending column indices
  integer :: begl, endl                       ! per-proc beginning and ending landunit indices
  integer :: begg, endg                       ! per-proc gridcell ending gridcell indices
  integer :: numg                             ! total number of gridcells across all processors
  integer :: numl                             ! total number of landunits across all processors
  integer :: numc                             ! total number of columns across all processors
  integer :: nump                             ! total number of pfts across all processors
  integer :: beg_lnd_rof,end_lnd_rof          ! beginning,ending landrunoff points
  integer :: beg_ocn_rof,end_ocn_rof          ! beginning,ending oceaan
!
! Send/recv buffers
!
  integer, parameter :: nsend_csm = cpl_fields_l2c_total
  integer, parameter :: nrecv_csm = cpl_fields_c2l_total
!
! Flux averaging arrays and counters
!
  integer  :: icnt                         ! step counter for flux averager
  integer  :: ncnt                         ! number of steps over which to average output fluxes
  real(r8) :: rncnt                        ! reciprocal of ncnt

  real(r8), allocatable :: taux_ave(:)     ! averaged array
  real(r8), allocatable :: tauy_ave(:)     ! averaged array
  real(r8), allocatable :: lhflx_ave(:)    ! averaged array
  real(r8), allocatable :: shflx_ave(:)    ! averaged array
  real(r8), allocatable :: lwup_ave(:)     ! averaged array
  real(r8), allocatable :: qflx_ave(:)     ! averaged array
  real(r8), allocatable :: swabs_ave(:)    ! averaged array
!
! When to send/receive messages to coupler and when to make restart and stop
!
  integer, private:: ncpday         ! number of send/recv calls per day
  logical, public :: dorecv         ! receive data from coupler this step
  logical, public :: dosend         ! send data to coupler this step
  logical, public :: csmstop_next   ! received stop at eod signal and will stop on next ts
  logical, public :: csmstop_now    ! received stop now signal from coupler
  logical, public :: csmrstrt       ! restart write signal received from coupler
!
! Indices for send/recv fields
!
  integer, parameter :: irecv_hgt    = cpl_fields_c2l_z
  integer, parameter :: irecv_u      = cpl_fields_c2l_u
  integer, parameter :: irecv_v      = cpl_fields_c2l_v
  integer, parameter :: irecv_th     = cpl_fields_c2l_ptem
  integer, parameter :: irecv_q      = cpl_fields_c2l_shum
  integer, parameter :: irecv_pbot   = cpl_fields_c2l_pbot
  integer, parameter :: irecv_t      = cpl_fields_c2l_tbot
  integer, parameter :: irecv_lwrad  = cpl_fields_c2l_lwdn
  integer, parameter :: irecv_rainc  = cpl_fields_c2l_rainc
  integer, parameter :: irecv_rainl  = cpl_fields_c2l_rainl
  integer, parameter :: irecv_snowc  = cpl_fields_c2l_snowc
  integer, parameter :: irecv_snowl  = cpl_fields_c2l_snowl
  integer, parameter :: irecv_soll   = cpl_fields_c2l_swndr
  integer, parameter :: irecv_sols   = cpl_fields_c2l_swvdr
  integer, parameter :: irecv_solld  = cpl_fields_c2l_swndf
  integer, parameter :: irecv_solsd  = cpl_fields_c2l_swvdf

  integer, parameter :: isend_trad   = cpl_fields_l2c_t
  integer, parameter :: isend_asdir  = cpl_fields_l2c_avsdr
  integer, parameter :: isend_aldir  = cpl_fields_l2c_anidr
  integer, parameter :: isend_asdif  = cpl_fields_l2c_avsdf
  integer, parameter :: isend_aldif  = cpl_fields_l2c_anidf
  integer, parameter :: isend_sno    = cpl_fields_l2c_snowh
  integer, parameter :: isend_taux   = cpl_fields_l2c_taux
  integer, parameter :: isend_tauy   = cpl_fields_l2c_tauy
  integer, parameter :: isend_lhflx  = cpl_fields_l2c_lat
  integer, parameter :: isend_shflx  = cpl_fields_l2c_sen
  integer, parameter :: isend_lwup   = cpl_fields_l2c_lwup
  integer, parameter :: isend_qflx   = cpl_fields_l2c_evap
  integer, parameter :: isend_tref2m = cpl_fields_l2c_tref
  integer, parameter :: isend_qref2m = cpl_fields_l2c_qref
  integer, parameter :: isend_swabs  = cpl_fields_l2c_swnet
!
! CCSM timers
!
  logical  :: timer_lnd_sendrecv = .false. ! true => timer is on
  logical  :: timer_lnd_recvsend = .false. ! true => timer is on
!
! Input dtime
!
  integer, public :: csm_dtime   ! value passed to coupler on initialization set to namelist input

! !PRIVATE MEMBER FUNCTIONS:
  private :: global_sum_fld2d    ! global sum of 2d grid fields
  private :: global_sum_fld1d    ! global sum of 1d fluxes
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_setup
!
! !INTERFACE:
  subroutine csm_setup(mpicom)
!
! !DESCRIPTION:
!  Initialize csm coupling, return the communicator group to the
!  application.
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: mpicom  !MPI group communicator
!
! !REVISION HISTORY:
!  03.01.15 T. Craig: initial version
!
!EOP
!
! !LOCAL VARIABLES:
!------------------------------------------------------------------------

   !call cpl_interface_init(cpl_fields_lndname,mpicom)
   call cpl_comm_init_wo_mph(cpl_fields_lndname, mpicom)

 end subroutine csm_setup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_shutdown
!
! !INTERFACE:
  subroutine csm_shutdown
!
! !DESCRIPTION:
!  Finalize csm coupling
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  03.01.15 T. Craig: initial version
!
!EOP
!
! !LOCAL VARIABLES:
!------------------------------------------------------------------------

   call cpl_interface_finalize(cpl_fields_lndname)

 end subroutine csm_shutdown

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_initialize
!
! !INTERFACE:
 subroutine csm_initialize(irad, eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION:
!  Send first control data to flux coupler and "invalid" grid
!  containing special value data.
!  The coupler treats points where the mask is nonzero as points where
!  you could possibly do a calculation (in the case of the runoff, this
!  corresponds to all rtm ocean points). The coupler then defines a "key"
!  as points where the model can give you valid data (in the case of runoff,
!  this corresponds to points where the land model will give you valid
!  compressed data points). The key can be 0 where the mask is 1. However,
!  the key cannot be 1 where the mask is 0 unless the data is also zero.
!  In the case of runoff, the key the coupler builds is time invariant.
!  Send first control data to flux coupler and "invalid" grid
!  containing special value data
!
! !USES:
    use clmtype
    use decompMod    , only : get_proc_bounds, get_proc_global
    use RunoffMod    , only : get_proc_rof_bounds, runoff
    use clm_varctl   , only : csm_doflxave, nsrest
    use RtmMod       , only : area_r, longxy_r, latixy_r, mask_r
    use clm_varcon   , only : re
    use time_manager , only : get_step_size
    use shr_const_mod, only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: irad   ! frequency of radiation computation
    real(r8), intent(out) :: eccen  ! Earth's eccentricity of orbit
    real(r8), intent(out) :: obliqr ! Earth's obliquity in radians
    real(r8), intent(out) :: lambm0 ! Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(out) :: mvelpp ! Earth's moving vernal equinox long of perihelion plus pi (radians)
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!  03.01.15 T.Craig Update for cpl6
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni              ! indices
    real(r8):: dtime                      ! step size
    real(r8):: spval                      ! special value
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

    ! Determine processor bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
    call get_proc_rof_bounds(beg_lnd_rof, end_lnd_rof, beg_ocn_rof, end_ocn_rof)

    ! Determine number of send/recv calls steps per day to flux coupler

    if (nsrest == 0) then
       dtime = get_step_size()
    else
       dtime = csm_dtime
    endif
    if (csm_doflxave) then
       ncpday = nint(SHR_CONST_CDAY/dtime)/irad
    else
       ncpday = nint(SHR_CONST_CDAY/dtime)
    endif

    ! Setup contracts for grid communication

    ibuffs(:)  = 0                                   ! initialize ibuffs
    ibuffs(cpl_fields_ibuf_gsize  ) = lsmlon*lsmlat  ! global array size
    ibuffs(cpl_fields_ibuf_gisize ) = lsmlon         ! global number of lons
    ibuffs(cpl_fields_ibuf_gjsize ) = lsmlat         ! global number of lats
    ibuffs(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
    ibuffs(cpl_fields_ibuf_ncpl   ) = ncpday         ! number of land send/recv calls per day
    ibuffs(cpl_fields_ibuf_inimask) = 1              ! T(or 1) => requests cpl to send valid land domain info back

    ibuffs(cpl_fields_ibuf_lsize  ) = endg-begg+1
    ibuffs(cpl_fields_ibuf_lisize ) = endg-begg+1
    ibuffs(cpl_fields_ibuf_ljsize ) = 1

    allocate(Gbuf(begg:endg,cpl_fields_grid_total))

    do g = begg,endg
       Gbuf(g,cpl_fields_grid_lon  ) = 1.0e30
       Gbuf(g,cpl_fields_grid_lat  ) = 1.0e30
       Gbuf(g,cpl_fields_grid_area ) = 1.0e30
       Gbuf(g,cpl_fields_grid_mask ) = 1.0e30
       gi = (gptr%jxy(g)-1)*lsmlon + gptr%ixy(g)
       Gbuf(g,cpl_fields_grid_index) = gi
    end do

    call cpl_interface_contractInit(contractS, cpl_fields_lndname,cpl_fields_cplname,cpl_fields_l2c_fields, ibuffs,Gbuf)
    call cpl_interface_contractInit(contractR, cpl_fields_lndname,cpl_fields_cplname,cpl_fields_c2l_fields, ibuffs,Gbuf)

    deallocate(Gbuf)

    ! Setup contracts for  runoff communication

    ibuffs(:) = 0
    ibuffs(cpl_fields_ibuf_gsize  ) = rtmlon*rtmlat                 ! global array size
    ibuffs(cpl_fields_ibuf_gisize ) = rtmlon                        ! global number of lons
    ibuffs(cpl_fields_ibuf_gjsize ) = rtmlat                        ! global number of lats
    ibuffs(cpl_fields_ibuf_ncpl   ) = ncpday                        ! number of land send/recv calls per day
    ibuffs(cpl_fields_ibuf_lsize  ) = end_ocn_rof - beg_ocn_rof + 1 ! local array size
    ibuffs(cpl_fields_ibuf_lisize ) = end_ocn_rof - beg_ocn_rof + 1 ! local array size
    ibuffs(cpl_fields_ibuf_ljsize ) = 1                             ! local array size

    csm_nptr = ibuffs(cpl_fields_ibuf_lsize)

    allocate(Gbuf(csm_nptr,cpl_fields_grid_total))

    do n = beg_ocn_rof,end_ocn_rof
       ni = n - beg_ocn_rof+ 1
       gi = (runoff%ocn_jxy(n)-1)*rtmlon + runoff%ocn_ixy(n)
       j = (gi-1) / rtmlon + 1
       i = mod(gi-1,rtmlon) + 1
       Gbuf(ni,cpl_fields_grid_lon  ) = longxy_r(i,j)
       Gbuf(ni,cpl_fields_grid_lat  ) = latixy_r(i,j)
       Gbuf(ni,cpl_fields_grid_area ) = area_r(i,j)/(re*re)
       Gbuf(ni,cpl_fields_grid_mask ) = 1.0 - float(mask_r(i,j))
       Gbuf(ni,cpl_fields_grid_index) = gi
    end do

    call cpl_interface_contractInit(contractSr,cpl_fields_lndname,cpl_fields_cplname,cpl_fields_r2c_fields, ibuffs,Gbuf)

    deallocate(Gbuf)

    ! Receive initial ibuf message

    call cpl_interface_ibufRecv(cpl_fields_cplname,ibuffr,rbuffr)

    spval  = rbuffr(cpl_fields_rbuf_spval)
    eccen  = rbuffr(cpl_fields_rbuf_eccen)
    obliqr = rbuffr(cpl_fields_rbuf_obliqr)
    lambm0 = rbuffr(cpl_fields_rbuf_lambm0)
    mvelpp = rbuffr(cpl_fields_rbuf_mvelpp)

    ! Check that data is good data and not the special value

    if (masterproc) then
       call compat_check_spval(spval, eccen ,'Eccentricity'     )
       call compat_check_spval(spval, obliqr,'Obliquity'        )
       call compat_check_spval(spval, lambm0,'Long of perhelion')
       call compat_check_spval(spval, mvelpp,'Move long of perh')

       write(6,*)'(CSM_INITIALIZE): eccen:  ', eccen
       write(6,*)'(CSM_INITIALIZE): obliqr: ', obliqr
       write(6,*)'(CSM_INITIALIZE): lambm0: ', lambm0
       write(6,*)'(CSM_INITIALIZE): mvelpp: ', mvelpp
    end if

    write(6,*)'(CSM_INITIALIZE): there will be ',ncpday, &
         ' send/recv calls per day from the land model to the flux coupler'
    write(6,*)'(CSM_INITIALIZE):sent l->d control data '

    ! Allocate memory for grid and runoff communication

    allocate(bufS(begg:endg,nsend_csm))
    allocate(bufR(begg:endg,nrecv_csm))
    allocate(bufSr(csm_nptr,cpl_fields_r2c_total))

  end subroutine csm_initialize

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_recvgrid
!
! !INTERFACE:
  subroutine csm_recvgrid(cam_longxy, cam_latixy, cam_numlon, &
                          cam_landfrac, cam_landmask)
!
! !DESCRIPTION:
!  Receive valid land grid and land mask from coupler
!
! !ARGUMENTS:
    implicit none
    integer , intent(out) :: cam_numlon(lsmlat)           !cam number of longitudes
    real(r8), intent(out) :: cam_longxy(lsmlon,lsmlat)    !cam lon values
    real(r8), intent(out) :: cam_latixy(lsmlon,lsmlat)    !cam lat values
    real(r8), intent(out) :: cam_landfrac(lsmlon,lsmlat)  !cam fractional land
    integer , intent(out) :: cam_landmask(lsmlon,lsmlat)  !cam land mask
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: i,j,n                    ! indices
    real(r8) :: area_a(lsmlon,lsmlat)    ! coupler atm grid areas
    integer  :: mask_a(lsmlon,lsmlat)    ! coupler atm valid grid mask
    real(r8) :: tmp(lsmlon,lsmlat)       ! temporary
!------------------------------------------------------------------------

    ! Set integer control information and setup contracts

    ibuffs(:)  = 0                                   !initialize ibuffs
    ibuffs(cpl_fields_ibuf_gsize  ) = lsmlon*lsmlat  !global array size
    ibuffs(cpl_fields_ibuf_gisize ) = lsmlon         !global number of lons
    ibuffs(cpl_fields_ibuf_gjsize ) = lsmlat         !global number of lats
    ibuffs(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
    ibuffs(cpl_fields_ibuf_inimask) = 1              !T(or 1) => requests cpl to send valid land domain info back
    ibuffs(cpl_fields_ibuf_lsize  ) = lsmlon*lsmlat  !local array size
    ibuffs(cpl_fields_ibuf_lisize ) = lsmlon         !local array size
    ibuffs(cpl_fields_ibuf_ljsize ) = lsmlat         !local array size

    csm_nptg = ibuffs(cpl_fields_ibuf_lsize)

    allocate(Gbuf(csm_nptg,cpl_fields_grid_total))
    do n = 1,csm_nptg
       Gbuf(n,cpl_fields_grid_lon  ) = 1.0e30
       Gbuf(n,cpl_fields_grid_lat  ) = 1.0e30
       Gbuf(n,cpl_fields_grid_area ) = 1.0e30
       Gbuf(n,cpl_fields_grid_mask ) = 1.0e30
       Gbuf(n,cpl_fields_grid_index) = n
    enddo
    call cpl_interface_contractInit(contractRg,cpl_fields_lndname,cpl_fields_cplname,cpl_fields_c2lg_fields,ibuffs,Gbuf)
    deallocate(Gbuf)

    allocate(Gbuf(csm_nptg,cpl_fields_c2lg_total))
    call cpl_interface_contractRecv(cpl_fields_cplname,contractRg,ibuffr,Gbuf)
    cam_numlon(:) = 0
    do n = 1,csm_nptg
       j = (n-1) / lsmlon + 1
       i = mod(n-1,lsmlon) + 1
       cam_longxy(i,j)   = Gbuf(n,cpl_fields_c2lg_alon)
       cam_latixy(i,j)   = Gbuf(n,cpl_fields_c2lg_alat)
       write(6,*) n,i,j,cam_latixy(i,j)
       area_a(i,j)       = Gbuf(n,cpl_fields_c2lg_aarea)
       cam_landfrac(i,j) = Gbuf(n,cpl_fields_c2lg_lfrac)
       cam_landmask(i,j) = nint(Gbuf(n,cpl_fields_c2lg_lmask))
       mask_a(i,j)       = nint(Gbuf(n,cpl_fields_c2lg_amask))
       if (mask_a(i,j) /= 0) cam_numlon(j) = cam_numlon(j)+1
    end do
    deallocate(Gbuf)

    if (masterproc) then
       write(6,*)'(CSM_SENDGRID):recd d->l land grid '
       write(6,100) 'lnd','recv', cpl_fields_c2lg_alon, global_sum_fld2d(cam_longxy  , 1.e30), ' lon'
       write(6,100) 'lnd','recv', cpl_fields_c2lg_alat, global_sum_fld2d(cam_latixy  , 1.e30), ' lat'
       write(6,100) 'lnd','recv', cpl_fields_c2lg_aarea,global_sum_fld2d(area_a      , 1.e30), ' aarea'
       write(6,100) 'lnd','recv', cpl_fields_c2lg_lfrac,global_sum_fld2d(cam_landfrac, 1.e30), ' lfrac'
       tmp = float(cam_landmask)
       write(6,100) 'lnd','recv', cpl_fields_c2lg_lmask,global_sum_fld2d(tmp, 1.e30), ' lmask'
       tmp = float(mask_a)
       write(6,100) 'lnd','recv', cpl_fields_c2lg_amask,global_sum_fld2d(tmp, 1.e30), ' amask'
100    format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
    endif

  end subroutine csm_recvgrid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_sendalb
!
! !INTERFACE:
  subroutine csm_sendalb()
!
! !DESCRIPTION:
! Send initial albedos, surface temperature and snow data to the
! flux coupler
!
! !USES:
    use clmtype
    use clm_varsur
    use clm_varctl  , only : csm_doflxave, nsrest
    use clm_varcon  , only : sb
    use time_manager, only : get_curr_date, get_prev_date, get_nstep
    use lnd2atmMod  , only : lnd2atm
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni ! generic indices
    integer :: yr            ! current year
    integer :: mon           ! current month
    integer :: day           ! current day (0, 1, ...)
    integer :: ncsec         ! current seconds of current date (0, ..., 86400)
    integer :: ncdate        ! current date (yymmdd format) (e.g., 021105)
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
! -----------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

    ! Fill send buffer for grid data

    bufS(:,:) = 0.0

    if (nsrest == 0 ) then   !initial run

       ! On initial timestep ONLY: determine 1d vector of states that will be sent
       ! to coupler and map fields from 1d subgrid vector to 2d [lsmlon]x[lsmlat] grid.

       call lnd2atm(init=.true.)

       do g = begg,endg
          bufS(g,isend_trad ) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb))
          bufS(g,isend_sno  ) = gptr%l2as%h2osno(g)
          bufS(g,isend_asdir) = gptr%l2as%albd(g,1)
          bufS(g,isend_aldir) = gptr%l2as%albd(g,2)
          bufS(g,isend_asdif) = gptr%l2as%albi(g,1)
          bufS(g,isend_aldif) = gptr%l2as%albi(g,2)
       end do

    else  ! restart run

       ! On a restart run, no meaningful data is sent to the flux coupler -
       ! this includes ocean runoff (which should only contain zero values)
       ! since the runoff code (riverfluxrtm) has not been called yet

       bufS(:,:) = 1.e30

    endif

    ! Determine time index to send to coupler. Note that for a restart run,
    ! the next time step is nstep+1. But must send current time step to flux couper here.

    if (nsrest == 0) then
       call get_curr_date (yr, mon, day, ncsec)
    else
       call get_prev_date (yr, mon, day, ncsec)
    endif

    ncdate = yr*10000 + mon*100 + day

    ibuffs(:)  = 0
    ibuffs(cpl_fields_ibuf_cdate)  = ncdate      !model date (yyyymmdd)
    ibuffs(cpl_fields_ibuf_sec  )  = ncsec       !elapsed seconds in current date

    ! Send grid data to coupler

    call cpl_interface_contractSend(cpl_fields_cplname,contractS,ibuffs,bufS)

    ! Send runoff data to coupler
    ! Must convert runoff to units of kg/m^2/s from m^3/s

    bufSr(:,:) = 0.
    do n = beg_ocn_rof, end_ocn_rof
       ni = n - beg_ocn_rof+ 1
       bufSr(ni,cpl_fields_r2c_runoff) = runoff%ocn(n) / (runoff%ocn_area(n)*1000.)
    end do
    call cpl_interface_contractSend(cpl_fields_cplname,contractSr,ibuffs,bufSr)

  end subroutine csm_sendalb

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_dosndrcv
!
! !INTERFACE:
  subroutine csm_dosndrcv(doalb)
!
! !DESCRIPTION:
! Determine when to send and receive messages to/from the
! flux coupler on this time-step.
! Determine if send/receive information to/from flux coupler
! Send msgs (land state and fluxes) to the flux coupler only when
! doalb is true (i.e. on time steps before the atm does a solar
! radiation computation). Receive msgs (atm state) from the
! flux coupler only when dorad is true (i.e. on time steps
! when the atm does a solar radiation computation).
! The fluxes are then averaged between the send and receive calls.
!
! !USES:
    use clm_varctl   , only : csm_doflxave
    use time_manager , only : get_step_size, get_nstep
    use shr_const_mod, only : SHR_CONST_CDAY
!
! !ARGUMENTS:
    implicit none
    logical, intent(in) :: doalb  !true=>next timestep a radiation time step
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ntspday           !model steps per day
    real(r8) :: dtime             !step size (seconds)
    integer  :: nstep             !time step
!-----------------------------------------------------------------------

    ! Determine if send/receive information to/from flux coupler

    nstep = get_nstep()
    if (csm_doflxave) then
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = doalb
       else
          dorecv = dosend
          dosend = doalb
       endif
    else
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = .true.
       else
          dorecv = .true.
          dosend = .true.
       endif
    endif

    ! If at end of day: check if should write restart file or stop
    ! at next time step. Note, these statements must appear here since
    ! ibuffr is not received at every time step when flux averaging occurs.

    csmstop_next = .false.
    csmrstrt     = .false.
    dtime        = get_step_size()
    ntspday      = nint(SHR_CONST_CDAY/dtime)
    if (mod(nstep,ntspday) == 0) then
       if (ibuffr(cpl_fields_ibuf_stopeod) /= 0) then  !stop at end of day
          csmstop_next = .true.  !will stop on next time step
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
       if (ibuffr(cpl_fields_ibuf_resteod) /= 0) then !write restart at end of day
          csmrstrt = .true.      !will write restart now
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
    endif

  end subroutine csm_dosndrcv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_recv
!
! !INTERFACE:
  subroutine csm_recv()
!
! !DESCRIPTION:
!  Receive and map data from flux coupler
!
! !USES:
    use clmtype
    use clm_varcon, only : rair, po2, pco2
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni ! generic indices
    real(r8):: forc_rainc    ! rainxy Atm flux mm/s
    real(r8):: forc_rainl    ! rainxy Atm flux mm/s
    real(r8):: forc_snowc    ! snowfxy Atm flux  mm/s
    real(r8):: forc_snowl    ! snowfxl Atm flux  mm/s
    integer :: ier           ! return error code
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

    ! Start timers

     if (timer_lnd_sendrecv) then
        call t_stopf ('lnd_sendrecv') ; timer_lnd_sendrecv = .false.
     endif

     call t_startf('lnd_recv')

     ibuffr(:) = 0

     call cpl_interface_contractRecv(cpl_fields_cplname,contractR,ibuffr,bufR)

     ! Do global integrals of fluxes if flagged

     if (ibuffr(cpl_fields_ibuf_infobug) >= 2) then
        if (masterproc) then
           if (.not. associated(bufRglob)) then
              allocate(bufRglob(nrecv_csm,numg), stat=ier)
              if (ier /= 0) then
                 write(6,*)'clm_csmMod: allocation error for bufRglob'; call endrun
              end if
              if (.not. associated(fieldR)) then
                 allocate(fieldR(numg), stat=ier)
                 if (ier /= 0) then
                    write(6,*)'clm_csmMod: allocation error for fieldR'; call endrun
                 end if
              endif
           end if
        end if
        if (.not. associated(bufRloc)) then
           allocate(bufRloc(nrecv_csm,begg:endg), stat=ier)
           if (ier /= 0) then
              write(6,*)'clm_csmMod: allocation error for bufRloc'; call endrun
           end if
        end if
        do g = begg,endg
           do n = 1,nrecv_csm
              bufRloc(n,g) = bufR(g,n)
           end do
        end do
#if (defined SPMD)
        call gather_data_to_master(bufRloc, bufRglob, clmlevel=nameg)
#else
        bufRglob(:,:) = bufRloc(:,:)
#endif
        if (masterproc) then
           write(6,*)
           fieldR(:) = bufRglob(irecv_hgt,:)
           write(6,100) 'lnd','recv', irecv_hgt  , global_sum_fld1d(fieldR), ' hgt'
           fieldR(:) = bufRglob(irecv_u,:)
           write(6,100) 'lnd','recv', irecv_u    , global_sum_fld1d(fieldR), ' u'
           fieldR(:) = bufRglob(irecv_v,:)
           write(6,100) 'lnd','recv', irecv_v    , global_sum_fld1d(fieldR), ' v'
           fieldR(:) = bufRglob(irecv_th,:)
           write(6,100) 'lnd','recv', irecv_th   , global_sum_fld1d(fieldR), ' th'
           fieldR(:) = bufRglob(irecv_q,:)
           write(6,100) 'lnd','recv', irecv_q    , global_sum_fld1d(fieldR), ' q'
           fieldR(:) = bufRglob(irecv_pbot,:)
           write(6,100) 'lnd','recv', irecv_pbot , global_sum_fld1d(fieldR), ' pbot'
           fieldR(:) = bufRglob(irecv_t,:)
           write(6,100) 'lnd','recv', irecv_t    , global_sum_fld1d(fieldR), ' t'
           fieldR(:) = bufRglob(irecv_lwrad,:)
           write(6,100) 'lnd','recv', irecv_lwrad, global_sum_fld1d(fieldR), ' lwrad'
           fieldR(:) = bufRglob(irecv_rainc,:)
           write(6,100) 'lnd','recv', irecv_rainc, global_sum_fld1d(fieldR), ' rainc'
           fieldR(:) = bufRglob(irecv_rainl,:)
           write(6,100) 'lnd','recv', irecv_rainl, global_sum_fld1d(fieldR), ' rainl'
           fieldR(:) = bufRglob(irecv_snowc,:)
           write(6,100) 'lnd','recv', irecv_snowc, global_sum_fld1d(fieldR), ' snowc'
           fieldR(:) = bufRglob(irecv_snowl,:)
           write(6,100) 'lnd','recv', irecv_snowl, global_sum_fld1d(fieldR), ' snowl'
           fieldR(:) = bufRglob(irecv_soll,:)
           write(6,100) 'lnd','recv', irecv_soll , global_sum_fld1d(fieldR), ' soll '
           fieldR(:) = bufRglob(irecv_sols,:)
           write(6,100) 'lnd','recv', irecv_sols , global_sum_fld1d(fieldR), ' sols '
           fieldR(:) = bufRglob(irecv_solld,:)
           write(6,100) 'lnd','recv', irecv_solld, global_sum_fld1d(fieldR), ' solld'
           fieldR(:) = bufRglob(irecv_solsd,:)
           write(6,100) 'lnd','recv', irecv_solsd, global_sum_fld1d(fieldR), ' solsd'
100        format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
           write(6,*)
        endif
     endif

     ! Stop timer

     call t_stopf('lnd_recv')

     ! Check if end of run now, if so stop (each processor does this)

     csmstop_now = .false.
     if (ibuffr(cpl_fields_ibuf_stopnow) /= 0) then
        csmstop_now = .true.
        if (timer_lnd_recvsend) call t_stopf('lnd_recvsend')
        if (timer_lnd_sendrecv) call t_stopf('lnd_sendrecv')
        write(6,*)'(CSM_RECV) stop now signal from flux coupler'
        write(6,*)'(CSM_RECV) ibuffr(cpl_fields_ibuf_stopnow) = ',ibuffr(cpl_fields_ibuf_stopnow)
        if (masterproc) then
           write(6,9001)
           write(6,9002) ibuffr(cpl_fields_ibuf_cdate)
           write(6,9003)
9001       format(/////' ===========> Terminating CLM Model')
9002       format(     '      Date: ',i8)
9003       format(/////' <=========== CLM Model Terminated')
        endif
        RETURN
     endif

     ! More timer logic

     if (.not. timer_lnd_recvsend) then
        call t_startf('lnd_recvsend') ; timer_lnd_recvsend = .true.
     endif

    ! Split data from coupler into component arrays.
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.


!$OMP PARALLEL DO PRIVATE (g,forc_rainc,forc_rainl,forc_snowc,forc_snowl)
!CSD$ PARALLEL DO PRIVATE (g,forc_rainc,forc_rainl,forc_snowc,forc_snowl)
    do g = begg,endg
       gptr%a2ls%forc_hgt(g)     = bufR(g,irecv_hgt  )  !zgcmxy  Atm state m
       gptr%a2ls%forc_u(g)       = bufR(g,irecv_u    )  !forc_uxy  Atm state m/s
       gptr%a2ls%forc_v(g)       = bufR(g,irecv_v    )  !forc_vxy  Atm state m/s
       gptr%a2ls%forc_th(g)      = bufR(g,irecv_th   )  !forc_thxy Atm state K
       gptr%a2ls%forc_q(g)       = bufR(g,irecv_q    )  !forc_qxy  Atm state kg/kg
       gptr%a2ls%forc_pbot(g)    = bufR(g,irecv_pbot )  !ptcmxy  Atm state Pa
       gptr%a2ls%forc_t(g)       = bufR(g,irecv_t    )  !forc_txy  Atm state K
       gptr%a2lf%forc_lwrad(g)   = bufR(g,irecv_lwrad)  !flwdsxy Atm flux  W/m^2
       forc_rainc                = bufR(g,irecv_rainc)  !mm/s
       forc_rainl                = bufR(g,irecv_rainl)  !mm/s
       forc_snowc                = bufR(g,irecv_snowc)  !mm/s
       forc_snowl                = bufR(g,irecv_snowl)  !mm/s
       gptr%a2lf%forc_solad(g,2) = bufR(g,irecv_soll )  !forc_sollxy  Atm flux  W/m^2
       gptr%a2lf%forc_solad(g,1) = bufR(g,irecv_sols )  !forc_solsxy  Atm flux  W/m^2
       gptr%a2lf%forc_solai(g,2) = bufR(g,irecv_solld)  !forc_solldxy Atm flux  W/m^2
       gptr%a2lf%forc_solai(g,1) = bufR(g,irecv_solsd)  !forc_solsdxy Atm flux  W/m^2

       ! Determine derived quantities

       gptr%a2ls%forc_hgt_u(g) = gptr%a2ls%forc_hgt(g)    !observational height of wind [m]
       gptr%a2ls%forc_hgt_t(g) = gptr%a2ls%forc_hgt(g)    !observational height of temperature [m]
       gptr%a2ls%forc_hgt_q(g) = gptr%a2ls%forc_hgt(g)    !observational height of humidity [m]
       gptr%a2ls%forc_vp(g)    = gptr%a2ls%forc_q(g) * gptr%a2ls%forc_pbot(g) &
                                 / (0.622 + 0.378 * gptr%a2ls%forc_q(g))
       gptr%a2ls%forc_rho(g)   = (gptr%a2ls%forc_pbot(g) - 0.378 * gptr%a2ls%forc_vp(g)) &
                                 / (rair * gptr%a2ls%forc_t(g))
       gptr%a2ls%forc_co2(g)   = pco2 * gptr%a2ls%forc_pbot(g)
       gptr%a2ls%forc_o2(g)    = po2 * gptr%a2ls%forc_pbot(g)
       gptr%a2ls%forc_wind(g)  = sqrt(gptr%a2ls%forc_u(g)**2 + gptr%a2ls%forc_v(g)**2)
       gptr%a2lf%forc_solar(g) = gptr%a2lf%forc_solad(g,1) + gptr%a2lf%forc_solai(g,1) + &
                                 gptr%a2lf%forc_solad(g,2) + gptr%a2lf%forc_solai(g,2)

       ! Determine precipitation needed by clm

       gptr%a2lf%forc_rain(g) = forc_rainc + forc_rainl
       gptr%a2lf%forc_snow(g) = forc_snowc + forc_snowl

    end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

  end subroutine csm_recv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_send
!
! !INTERFACE:
  subroutine csm_send()
!
! !DESCRIPTION:
! Send data to the flux coupler
!
! !USES:
    use clmtype
    use clm_varctl  , only : csm_doflxave
    use clm_varsur  , only : landmask
    use clm_varcon  , only : sb
    use time_manager, only : get_curr_date, get_nstep
    use lnd2atmMod  , only : lnd2atm
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,n,ni,g,c,p  ! indices
    integer :: yr              ! current year
    integer :: mon             ! current month
    integer :: day             ! current day (0, 1, ...)
    integer :: ncsec           ! current seconds of current date (0, ..., 86400)
    integer :: ncdate          ! current date (yymmdd format) (e.g., 021105)
    integer :: ier             ! error status
    real(r8):: wt              ! weight
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived type
    type(pft_type)     , pointer :: pptr  ! pointer to column   derived type
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    pptr => clm3%g%l%c%p

    ! Send data to the flux coupler

    if (timer_lnd_recvsend) then
       call t_stopf ('lnd_recvsend') ; timer_lnd_recvsend = .false.
    endif

    ! Start timer

    call t_startf('lnd_send')

    ! Determine 1d vector of fields that will be sent to coupler.
    ! Coupler has convention that fluxes are positive downward.

    bufS(:,:) = 0.0

    ! Calculate fluxes and states to send to atm

    call lnd2atm()

    ! Recalculate fluxes to send to atm in flux averaged case
    ! Also recalculate a new t_rad based on flux averged lwup

    if (csm_doflxave) then
       do g = begg,endg
          gptr%l2af%taux(g) = 0.
          gptr%l2af%tauy(g) = 0.
          gptr%l2af%eflx_lh_tot(g) = 0.
          gptr%l2af%eflx_sh_tot(g) = 0.
          gptr%l2af%eflx_lwrad_out(g) = 0.
          gptr%l2af%qflx_evap_tot(g) = 0.
          gptr%l2af%fsa(g) = 0.
       end do
       do p = begp,endp
          g  = pptr%gridcell(p)
          wt = pptr%wtgcell(p)
          gptr%l2af%taux(g)           = gptr%l2af%taux(g)           + taux_ave(p)  * wt
          gptr%l2af%tauy(g)           = gptr%l2af%tauy(g)           + tauy_ave(p)  * wt
          gptr%l2af%eflx_lh_tot(g)    = gptr%l2af%eflx_lh_tot(g)    + lhflx_ave(p) * wt
          gptr%l2af%eflx_sh_tot(g)    = gptr%l2af%eflx_sh_tot(g)    + shflx_ave(p) * wt
          gptr%l2af%eflx_lwrad_out(g) = gptr%l2af%eflx_lwrad_out(g) + lwup_ave(p)  * wt
          gptr%l2af%qflx_evap_tot(g)  = gptr%l2af%qflx_evap_tot(g)  + qflx_ave(p)  * wt
          gptr%l2af%fsa(g)            = gptr%l2af%fsa(g)            + swabs_ave(p) * wt
       end do
       do g = begg,endg
          gptr%l2as%t_rad(g) = (abs(gptr%l2af%eflx_lwrad_out(g)/sb))**0.25
       end do
    endif

    do g = begg,endg
       bufS(g,isend_trad  ) =  gptr%l2as%t_rad(g)
       bufS(g,isend_sno   ) =  gptr%l2as%h2osno(g)
       bufS(g,isend_asdir ) =  gptr%l2as%albd(g,1)
       bufS(g,isend_aldir ) =  gptr%l2as%albd(g,2)
       bufS(g,isend_asdif ) =  gptr%l2as%albi(g,1)
       bufS(g,isend_aldif ) =  gptr%l2as%albi(g,2)
       bufS(g,isend_tref2m) =  gptr%l2as%t_ref2m(g)
       bufS(g,isend_qref2m) =  gptr%l2as%q_ref2m(g)
       bufS(g,isend_taux  ) = -gptr%l2af%taux(g)
       bufS(g,isend_tauy  ) = -gptr%l2af%tauy(g)
       bufS(g,isend_lhflx ) = -gptr%l2af%eflx_lh_tot(g)
       bufS(g,isend_shflx ) = -gptr%l2af%eflx_sh_tot(g)
       bufS(g,isend_lwup  ) = -gptr%l2af%eflx_lwrad_out(g)
       bufS(g,isend_qflx  ) = -gptr%l2af%qflx_evap_tot(g)
       bufS(g,isend_swabs ) = -gptr%l2af%fsa(g)
    end do

    ! Map fields from 1d-grid vector to 2d-grid
    ! NOTE: snow is sent as zero over non-land because currently
    ! the ocn and sea-ice send no snow cover to cpl and so the cpl
    ! sends back zero snow over non-land to  the atm (the atm and
    ! land grid are currently assumed to be identical)

    call get_curr_date (yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day

    ibuffs(:)  = 0
    ibuffs(cpl_fields_ibuf_cdate)  = ncdate      ! model date (yyyymmdd)
    ibuffs(cpl_fields_ibuf_sec  )  = ncsec       ! elapsed seconds in current date

    call cpl_interface_contractSend(cpl_fields_cplname,contractS ,ibuffs,bufS)

    ! Must convert runoff to units of kg/m^2/s from m^3/s

    bufSr(:,:) = 0.
    do n = beg_ocn_rof, end_ocn_rof
       ni = n - beg_ocn_rof+ 1
       bufSr(ni,cpl_fields_r2c_runoff) = runoff%ocn(n) / (runoff%ocn_area(n)*1000.)
    enddo
    call cpl_interface_contractSend(cpl_fields_cplname,contractSr,ibuffs,bufSr)

    ! Do global integrals if flag is set

    if (ibuffr(cpl_fields_ibuf_infobug) >= 2) then
       if (masterproc) then
          if (.not. associated(bufSglob)) then
             allocate(bufSglob(nsend_csm,numg), stat=ier)
             if (ier /= 0) then
                write(6,*)'clm_csmMod: allocation error for bufSglob'; call endrun
             end if
          end if
          if (.not. associated(fieldS)) then
             allocate(fieldS(numg), stat=ier)
             if (ier /= 0) then
                write(6,*)'clm_csmMod: allocation error for fieldS'; call endrun
             end if
          endif
       end if
       if (.not. associated(bufSloc)) then
          allocate(bufSloc(nsend_csm,begg:endg), stat=ier)
          if (ier /= 0) then
             write(6,*)'clm_csmMod: allocation error for bufSloc'; call endrun
          end if
       end if
       do g = begg,endg
          do n = 1,nsend_csm
             bufSloc(n,g) = bufS(g,n)
          end do
       end do
#if (defined SPMD)
       call gather_data_to_master(bufSloc, bufSglob, clmlevel=nameg)
#else
       bufSglob(:,:) = bufSloc(:,:)
#endif
       if (masterproc) then
          write(6,*)
          fieldS(:) = bufSglob(isend_trad,:)
          write(6,100) 'lnd','send', isend_trad , global_sum_fld1d(fieldS), ' trad'
          fieldS(:) = bufSglob(isend_asdir,:)
          write(6,100) 'lnd','send', isend_asdir, global_sum_fld1d(fieldS),' asdir'
          fieldS(:) = bufSglob(isend_aldir,:)
          write(6,100) 'lnd','send', isend_aldir, global_sum_fld1d(fieldS),' aldir'
          fieldS(:) = bufSglob(isend_asdif,:)
          write(6,100) 'lnd','send', isend_asdif, global_sum_fld1d(fieldS),' asdif'
          fieldS(:) = bufSglob(isend_aldif,:)
          write(6,100) 'lnd','send', isend_aldif, global_sum_fld1d(fieldS),' aldif'
          fieldS(:) = bufSglob(isend_taux,:)
          write(6,100) 'lnd','send', isend_taux , global_sum_fld1d(fieldS), ' taux'
          fieldS(:) = bufSglob(isend_tauy,:)
          write(6,100) 'lnd','send', isend_tauy , global_sum_fld1d(fieldS), ' tauy'
          fieldS(:) = bufSglob(isend_lhflx,:)
          write(6,100) 'lnd','send', isend_lhflx, global_sum_fld1d(fieldS), ' lhflx'
          fieldS(:) = bufSglob(isend_shflx,:)
          write(6,100) 'lnd','send', isend_shflx, global_sum_fld1d(fieldS), ' shflx'
          fieldS(:) = bufSglob(isend_lwup,:)
          write(6,100) 'lnd','send', isend_lwup , global_sum_fld1d(fieldS), ' lwup'
          fieldS(:) = bufSglob(isend_qflx,:)
          write(6,100) 'lnd','send', isend_qflx , global_sum_fld1d(fieldS), ' qflx'
          fieldS(:) = bufSglob(isend_swabs,:)
          write(6,100) 'lnd','send', isend_swabs, global_sum_fld1d(fieldS), ' swabs'
          write(6,*)
100       format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
       endif
    endif

    ! Stop timers

    call t_stopf('lnd_send')

    if (.not. timer_lnd_recvsend) then
       call t_startf('lnd_sendrecv') ; timer_lnd_sendrecv = .true.
    endif

  end subroutine csm_send

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_flxave
!
! !INTERFACE:
  subroutine csm_flxave()
!
! !DESCRIPTION:
! Average output fluxes for flux coupler
! Add land surface model output fluxes to accumulators every time step.
! When icnt==ncnt, compute the average flux over the time interval.
!
! !USES:
!
    use clmtype
    use clm_varctl  , only : irad
    use time_manager, only : get_nstep
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: nstep        ! model time step
    type(pft_type), pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    pptr => clm3%g%l%c%p

    ! Allocate dynamic memory if necessary
    if (.not. allocated(taux_ave)) then
       allocate (taux_ave(begp:endp)) ; taux_ave(:) = nan
    endif
    if (.not. allocated(tauy_ave)) then
       allocate (tauy_ave(begp:endp)) ; tauy_ave(:) = nan
    endif
    if (.not. allocated(lhflx_ave)) then
       allocate (lhflx_ave(begp:endp)); lhflx_ave(:) = nan
    endif
    if (.not. allocated(shflx_ave)) then
       allocate (shflx_ave(begp:endp)); shflx_ave(:) = nan
    endif
    if (.not. allocated(lwup_ave)) then
       allocate (lwup_ave(begp:endp)) ; lwup_ave(:) = nan
    endif
    if (.not. allocated(qflx_ave)) then
       allocate (qflx_ave(begp:endp)) ; qflx_ave(:) = nan
    endif
    if (.not. allocated(swabs_ave)) then
       allocate (swabs_ave(begp:endp)) ; swabs_ave(:) = nan
    endif

    ! Determine output flux averaging interval

    nstep = get_nstep()
    if (dorecv) then
       icnt = 1
       if ( nstep==0 ) then
          ncnt = irad + 1
       else
          ncnt = irad
       endif
       rncnt = 1./ncnt
    endif

    if (icnt == 1) then

       ! Initial call of averaging interval, copy data to accumulators

       do p = begp,endp
          taux_ave(p)  = pptr%pmf%taux(p)
          tauy_ave(p)  = pptr%pmf%tauy(p)
          lhflx_ave(p) = pptr%pef%eflx_lh_tot(p)
          shflx_ave(p) = pptr%pef%eflx_sh_tot(p)
          lwup_ave(p)  = pptr%pef%eflx_lwrad_out(p)
          qflx_ave(p)  = pptr%pwf%qflx_evap_tot(p)
          swabs_ave(p) = pptr%pef%fsa(p)
       end do

    else if (icnt == ncnt) then

       ! Final call of averaging interval, complete averaging

       do p = begp,endp
          taux_ave(p)  = rncnt * (taux_ave(p)  + pptr%pmf%taux(p))
          tauy_ave(p)  = rncnt * (tauy_ave(p)  + pptr%pmf%tauy(p))
          lhflx_ave(p) = rncnt * (lhflx_ave(p) + pptr%pef%eflx_lh_tot(p))
          shflx_ave(p) = rncnt * (shflx_ave(p) + pptr%pef%eflx_sh_tot(p))
          lwup_ave(p)  = rncnt * (lwup_ave(p)  + pptr%pef%eflx_lwrad_out(p))
          qflx_ave(p)  = rncnt * (qflx_ave(p)  + pptr%pwf%qflx_evap_tot(p))
          swabs_ave(p) = rncnt * (swabs_ave(p) + pptr%pef%fsa(p))
       end do

    else

       ! Intermediate call, add data to accumulators

       do p = begp,endp
          taux_ave(p)  = (taux_ave(p)  + pptr%pmf%taux(p))
          tauy_ave(p)  = (tauy_ave(p)  + pptr%pmf%tauy(p))
          lhflx_ave(p) = (lhflx_ave(p) + pptr%pef%eflx_lh_tot(p))
          shflx_ave(p) = (shflx_ave(p) + pptr%pef%eflx_sh_tot(p))
          lwup_ave(p)  = (lwup_ave(p)  + pptr%pef%eflx_lwrad_out(p))
          qflx_ave(p)  = (qflx_ave(p)  + pptr%pwf%qflx_evap_tot(p))
          swabs_ave(p) = (swabs_ave(p) + pptr%pef%fsa(p))
       end do

    end if

    ! Increment counter

    icnt = icnt + 1

  end subroutine csm_flxave

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: compat_check
!
! !INTERFACE:
  subroutine compat_check_spval(spval, data, string)
!
! !DESCRIPTION:
! Check that the given piece of real data sent from the coupler is valid
! data and not the couplers special data flag.  This ensures that the data
! you expect is actually being sent by the coupler.
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: spval
    real(r8), intent(in) :: data
    character(len=*), intent(in) :: string
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    if ( spval == data )then
       write(6,*)'ERROR:(compat_check_spval) msg incompatibility'
       write(6,*)'ERROR: I expect to recieve the data type: ',string
       write(6,*)'from CPL, but all I got was the special data flag'
       write(6,*)'coupler must not be sending this data, you are'
       write(6,*)'running with an incompatable version of the coupler'
       call endrun
    end if

  end subroutine compat_check_spval

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_compat
!
! !INTERFACE:
  subroutine csm_compat(cpl_maj_vers, cpl_min_vers, expect_maj_vers, &
                        expect_min_vers)
!
! !DESCRIPTION:
! Checks that the message recieved from the coupler is compatable
! with the type of message that I expect to recieve.  If the minor
! version numbers differ I print a warning message.  If the major
! numbers differ I abort since that means that the change is
! drastic enough that I can't run with the differences.
! Original Author: Erik Kluzek Dec/97
!
! !PARAMETERS:
    implicit none
    integer, intent(in) :: cpl_maj_vers    ! major version from coupler initial ibuffr array
    integer, intent(in) :: cpl_min_vers    ! minor version from coupler initial ibuffr array
    integer, intent(in) :: expect_maj_vers ! major version of the coupler I'm expecting
    integer, intent(in) :: expect_min_vers ! minor version of the coupler I'm expecting
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    write(6,*)'(cpl_COMPAT): This is revision: $Revision: 1.11.4.37 $'
    write(6,*)'              Tag: $Name: ccsm3_0_rel04 $'
    write(6,*)'              of the message compatability interface:'

    if ( cpl_min_vers /= expect_min_vers )then
       write(6,*) 'WARNING(cpl_compat):: Minor version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_min_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_min_vers
    end if

    if ( cpl_maj_vers /= expect_maj_vers )then
       write(6,*) 'ERROR(cpl_compat):: Major version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_maj_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_maj_vers
       call endrun
    end if

  end subroutine csm_compat

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_coupler
!
! !INTERFACE:
  subroutine restart_coupler(nio, flag)
!
! !DESCRIPTION:
!  Read/write restart data needed for running in flux coupled mode
!
! !USES:
    use clm_varctl, only : csm_doflxave
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit
    character(len=*), intent(in) :: flag   !"read" or "write"
!
! !REVISION HISTORY:
!  02.09.17  Mariana Vertenstein: moved code to be part of ccsm module
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: flxave_res   !flux averaging flag read from restart file
    integer :: ier          !mpi return error code
!-----------------------------------------------------------------------

    if (flag == 'read') then
       if (masterproc) then
          read(nio) flxave_res
          read(nio) dosend
          if ((flxave_res .and. .not.csm_doflxave).or.(.not.flxave_res .and. csm_doflxave)) then
             write(6,*)'(RESTART_COUPLER): flxave value from namelist ',csm_doflxave, &
                  ' must be the same as flxave value from restart dataset ',flxave_res
             call endrun
          endif
          if (flxave_res .and. .not. dosend) then
             write(6,*)'(RESTART_COUPLER): assume that current flux coupled model ', &
                  'with flux averaging must stop on a time step where dosend (doalb) is true'
             call endrun
          end if
       endif
#if ( defined SPMD )
       call mpi_bcast(dosend    , 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= MPI_SUCCESS) then
          write(6,*) 'MPI BCAST ERROR: for dosend in restart_coupler'
          call endrun
       end if
       call mpi_bcast(flxave_res, 1, MPI_INTEGER, 0, mpicom, ier)
       if (ier /= MPI_SUCCESS) then
          write(6,*) 'MPI BCAST ERROR: for flxave_res in restart_coupler'
          call endrun
       end if
#endif
    endif

    if (flag == 'write') then
       if (masterproc) then
          write(nio) csm_doflxave
          write(nio) dosend
       endif
    end if

  end subroutine restart_coupler

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: global_sum_fld2d
!
! !INTERFACE:
  real(r8) function global_sum_fld2d(array, spval)
!
! !DESCRIPTION:
! Performs a global sum on an input 2d grid array
!
! !USES:
    use clm_varsur, only : area                 !km^2
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: array(lsmlon,lsmlat) !W/m2, Kg/m2-s or N/m2
    real(r8), intent(in) :: spval                !points to not include in global sum
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j         !indices
!------------------------------------------------------------------------

    global_sum_fld2d = 0.
    do j = 1,lsmlat
       do i = 1,lsmlon
          if (array(i,j) /= spval) then
             global_sum_fld2d = global_sum_fld2d + array(i,j) * area(i,j) * 1.e6
          endif
       end do
    end do

  end function global_sum_fld2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: global_sum_fld1d
!
! !INTERFACE:
  real(r8) function global_sum_fld1d(array)
!
! !DESCRIPTION:
! Performs a global sum on an input flux array
!
! !USES:
    use clmtype
    use clm_varsur, only : area
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: array(:) !W/m2, Kg/m2-s or N/m2
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g,i,j  ! indices
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived type
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

    ! Note: area is in km^2

    global_sum_fld1d = 0.
    do g = 1,numg
       i = gptr%ixy(g)
       j = gptr%jxy(g)
       global_sum_fld1d = global_sum_fld1d + array(g) * area(i,j) * 1.e6
    end do

  end function global_sum_fld1d

#endif

end module clm_csmMod
