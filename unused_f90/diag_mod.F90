!===============================================================================
! CVS: $Id: diag_mod.F90,v 1.1.1.1 2005/02/03 22:29:00 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/diag_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: diag_mod -- computes spatial \& time averages of fluxed quatities
!
! !DESCRIPTION:
!    The coupler is required to do certain diagnostics, those calculations are
!    located in this module.
!
! !REMARKS:
!    Sign convention:
!       positive value <=> the model is gaining water, heat, momentum, etc.
!    Unit convention:
!       heat flux     ~ W/m^2 
!       momentum flux ~ N/m^2
!       water flux    ~ (kg/s)/m^2
!       salt  flux    ~ (kg/s)/m^2
!
! !REVISION HISTORY:
!    199x-mmm-dd - B. Kauffman - original cpl5 version
!    2002-nov-21 - R. Jacob - initial port to cpl6. Does atm and lnd
!    2002-nov-27 - R. Jacob - add ocean
!    2002-dec-03 - R. Jacob - add solar diagnostics
!    2002-dec-15 - R. Jacob - time average diagnostics
!    2003-Feb-10 - R. Jacob - calculate sums locally
!
! !INTERFACE: ------------------------------------------------------------------

module diag_mod

! !USES:

   use shr_date_mod    ! shared date module
   use shr_sys_mod     ! shared system routines
   use shr_timer_mod   ! shared timers
   use shr_mpi_mod     ! shared mpi layer
   use cpl_kind_mod    ! kinds
   use cpl_const_mod   ! physical constants
   use cpl_mct_mod     ! mct library
   use cpl_comm_mod    ! communicator module
   use cpl_fields_mod  ! index to fields in bundles
   use cpl_domain_mod  ! domain data types
   use cpl_bundle_mod  ! bundle data types
   use cpl_control_mod ! cpl control flags & methods

   implicit none

   private ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: diag_doDiag   ! coordinates all diagnostic subroutines
   public :: diag_solar    ! verifies net-solar coordination

! !PUBLIC DATA MEMBERS:

   !--- note: this partial-sum data needs to be saved in a restart file ---
   real(R8),save,public        :: diag_eday0        ! partial sum: start date
   real(R8),save,public        :: diag_eday1        ! partial sum: end   date
   real(R8),save,public        :: diag_ns           ! partial sum: number of samples
   real(R8),save,public,target :: diag_datas(8,6,3) ! partial sum: the p-sum data

!EOP

   !----- local constants -----
   real(R8),parameter :: HFLXtoWFLX =  &        ! converts heat to water
   &  - (cpl_const_ocn_ref_sal-cpl_const_ice_ref_sal) / &
   &    (cpl_const_ocn_ref_sal*cpl_const_latice)     

   !----- diagnostic data array -----
   real(R8),target :: diag_datai (8,6,3) ! data: instantaneous (n,m,p), local
   real(R8)        :: diag_Gdatai(8,6,3) ! data: instantaneous (n,m,p), global sum

   !----- indicies into data array -----
   integer(IN),parameter :: p_heat    = 1 ! group index: heat fluxes
   integer(IN),parameter :: p_water   = 2 ! group index: water fluxes
   integer(IN),parameter :: p_area    = 3 ! group index: area

   integer(IN),parameter :: m_atm     = 1 ! model index: atm
   integer(IN),parameter :: m_ice_nh  = 2 ! model index: ice, northern
   integer(IN),parameter :: m_ice_sh  = 3 ! model index: ice, southern
   integer(IN),parameter :: m_lnd     = 4 ! model index: lnd
   integer(IN),parameter :: m_ocn     = 5 ! model index: ocn
   integer(IN),parameter :: m_sum     = 6 ! model index: sum of all

   integer(IN),parameter :: n_area    = 1 ! area (wrt to unit sphere)

   integer(IN),parameter :: n_hfrz    = 1 ! heat : latent, freezing 
   integer(IN),parameter :: n_hmelt   = 2 ! heat : latent, melting  
   integer(IN),parameter :: n_hswnet  = 3 ! heat : short wave, net
   integer(IN),parameter :: n_hlwdn   = 4 ! heat : longwave down
   integer(IN),parameter :: n_hlwup   = 5 ! heat : longwave up
   integer(IN),parameter :: n_hlat    = 6 ! heat : latent, vaporization
   integer(IN),parameter :: n_hsen    = 7 ! heat : sensible
   integer(IN),parameter :: n_hnet    = 8 ! heat : sum of all heat

   integer(IN),parameter :: n_wfrz    = 1 ! water: freezing
   integer(IN),parameter :: n_wmelt   = 2 ! water: melting
   integer(IN),parameter :: n_wrain   = 3 ! water: precip, liquid
   integer(IN),parameter :: n_wsnow   = 4 ! water: precip, frozen
   integer(IN),parameter :: n_wevap   = 5 ! water: evaporation
   integer(IN),parameter :: n_wroff   = 6 ! water: runoff
   integer(IN),parameter :: n_wnet    = 7 ! water: sum of all water

   !----- used for short-wave net verification -----
   real(R8),save :: swnet_atm       ! cpl's expected swnet from atm
   real(R8),save :: swnet_lnd       ! cpl's expected swnet from lnd
   real(R8),save :: swnet_ice_nh    ! cpl's expected swnet from ice, northern hemi
   real(R8),save :: swnet_ice_sh    ! cpl's expected swnet from ice, southern hemi
   real(R8),save :: swnet_ocn       ! cpl's expected swnet from ocn

   !--- module variables ---
   integer(IN),parameter :: pid0 = 0        ! root process pid = zero
   integer(IN),save      :: tda,tdl,tdi,tdo ! timers for atm/lnd/ice/ocn

   !----- formats -----
   character(*),parameter :: F10="('(diag) ',a17,'    date    sec    area',&
   &       '    freeze      melt     netsw      lwdn',                     &
   &       '      lwup       lat       sen       net W/m^2')"
   character(*),parameter :: F11="('(diag) ',a17,' -------- ----- -------',&
   &       8(' ---------'))"
   character(*),parameter :: F12="('(diag) ',a17,i9.8,i6,f8.5,8(f10.4))"

   character(*),parameter :: F20="('(diag) ',a17,'    date    sec    area',&
   &       '    freeze      melt      rain      snow',                     &
   &       '      evap   run-off       net 10^-6 (kg/s)/m^2')"
   character(*),parameter :: F21="('(diag) ',a17,' -------- ----- -------',&
   &       7(' ---------'))"
   character(*),parameter :: F22="('(diag) ',a17,i9.8,i6,f8.5,7(f10.5))"

   character(*),parameter :: F30="('(diag) ',a17,'    date  ndays    area',&
   &       '    freeze      melt     netsw      lwdn',                     &
   &       '      lwup       lat       sen       net W/m^2')"
   character(*),parameter :: F40="('(diag) ',a17,'    date  ndays    area',&
   &       '    freeze      melt      rain      snow',                     &
   &       '      evap   run-off       net 10^-6 (kg/s)/m^2')"

   character(*),parameter :: modName = "diag_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_doDiag - coordinates calculation of diagnostic data
!
! !DESCRIPTION:
!    Calculate global diagnostics.  
!
! !REMARKS:
!    if (cpl_control_diagNow  ) then print instantaneous diagnostics
!    if (cpl_control_avDiagNow) then print time-avg diagnostics
!
!    This is hard-coded to print/reset the t-avg data at the end of every year.
!
! !REVISION HISTORY:
!   199x-mmm-dd - B. Kauffman original cpl5 version, called diagnos in diag_mod
!
! !INTERFACE: ------------------------------------------------------------------

subroutine diag_doDiag(date,bun_a2c,bun_c2a ,bun_l2c ,bun_c2l ,bun_r2c  , &
           &                bun_i2c,bun_c2i ,bun_o2c ,bun_c2o ,bun_a2c_o, &
           &                bun_alb,bun_lfrac,bun_ifrac,bun_ofrac)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(shr_date)  ,intent(in) :: date       ! current model date
   type(cpl_bundle),intent(in) :: bun_a2c    ! atm->cpl bundle
   type(cpl_bundle),intent(in) :: bun_c2a    ! cpl->atm bundle
   type(cpl_bundle),intent(in) :: bun_l2c    ! lnd->cpl bundle
   type(cpl_bundle),intent(in) :: bun_c2l    ! cpl->lnd bundle
   type(cpl_bundle),intent(in) :: bun_r2c    ! rof->cpl bundle
   type(cpl_bundle),intent(in) :: bun_i2c    ! ice->cpl bundle
   type(cpl_bundle),intent(in) :: bun_c2i    ! cpl->ice bundle
   type(cpl_bundle),intent(in) :: bun_o2c    ! ocn->cpl bundle
   type(cpl_bundle),intent(in) :: bun_c2o    ! cpl->ocn bundle
   type(cpl_bundle),intent(in) :: bun_a2c_o  ! atm->cpl bundle
   type(cpl_bundle),intent(in) :: bun_alb    ! albedo   bundle

   type(cpl_bundle),intent(in) :: bun_lfrac  ! surface fractions on lnd domain
   type(cpl_bundle),intent(in) :: bun_ifrac  ! surface fractions on ice domain
   type(cpl_bundle),intent(in) :: bun_ofrac  ! surface fractions on ocn domain

!EOP

   !----- local -----
   integer(IN)    :: year,month,day,sec,eday ! date info
   integer(IN)    :: diagsize                ! data size for mpi reduce
   integer(IN)    :: rcode                   ! return code
   logical        :: first_call=.true.       ! flags initialization
   type(shr_date) :: nextDate                ! date after advancing 1 time step

   !----- formats -----
   character(*),parameter :: subName = '(diag_doDiag) '
   character(*),parameter :: F00   = "('(diag_doDiag) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !--- get current elapsed day & seconds ---
   call shr_date_getEDay(date,eday,sec)

   !----------------------------------------------------------------------------
   ! one-time initializations
   !----------------------------------------------------------------------------
   if (first_call) then
      if (cpl_control_restType == 'initial' ) then ! initialize partial sum
         diag_datai = 0.0  ! instantaneous
         diag_datas = 0.0  ! partial sum
         diag_eday0 = eday ! 1st  day in partial sum
         diag_eday1 = -1   ! last day in partial sum
         diag_ns    =  0   ! number of samples in partial sum
      end if
      call shr_timer_get(tda,'diag_atm')
      call shr_timer_get(tdl,'diag_lnd')
      call shr_timer_get(tdi,'diag_ice')
      call shr_timer_get(tdo,'diag_ocn')
   end if

   !----------------------------------------------------------------------------
   ! compute instantaneous diagnostics
   !----------------------------------------------------------------------------
   call shr_timer_start(tda)
   call diag_atm(bun_a2c,bun_c2a)
   call shr_timer_stop(tda)

   call shr_timer_start(tdl)
   call diag_lnd(bun_l2c,bun_c2l,bun_r2c,bun_lfrac)
   call shr_timer_stop(tdl)

   call shr_timer_start(tdi)
   call diag_ice(bun_i2c,bun_c2i,bun_o2c,bun_ifrac)
   call shr_timer_stop(tdi)

   call shr_timer_start(tdo)
   call diag_ocn(bun_o2c,bun_c2o,bun_a2c_o,bun_i2c,bun_alb,bun_ofrac)
   call shr_timer_stop(tdo)

   !--- sum of instantaneous, inter-model data ---
   diag_datai(:,m_sum,:      ) = 0.0
   diag_datai(:,m_sum,p_heat ) = sum(diag_datai(:,:,p_heat ),DIM=2)
   diag_datai(:,m_sum,p_water) = sum(diag_datai(:,:,p_water),DIM=2)
   diag_datai(:,m_sum,p_area ) = sum(diag_datai(:,:,p_area ),DIM=2)

   !--- reduce instantaneous data onto master process ---
   diagsize=size(diag_datai)
   call shr_mpi_sum(diag_datai,diag_Gdatai,cpl_comm_comp,subName//" diag_datai")

   !----------------------------------------------------------------------------
   ! add instantaneous data to partial sum, printout & re-init as necessary
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid == 0) then

      !--- finish the global integral ---
      diag_Gdatai = diag_Gdatai/(4*cpl_const_pi)

      !--- add instantaneous data to partial sum ---
      diag_datas = diag_datas + diag_Gdatai
      diag_ns    = diag_ns + 1
      diag_eday1 = eday + sec/86400.0

      !--- is it the very end of a year? ---
      nextDate = date  ! will be the date at very next time step
      call shr_date_adv1Step(nextDate)
      call shr_date_getYMD(nextDate,year,month,day,sec)

      !--- print data to stdout ---
      if (cpl_control_diagNow  ) call diag_print   (date)
      if (cpl_control_avdiagNow) then
         call diag_printAvg(date)
      else if ( month==1 .and. day==1 .and. sec==0) then ! => end-of-year
         call diag_printAvg(date)
      end if

      !--- re-initialize partial sum data at *end* of year ---
      if ( month==1 .and. day==1 .and. sec==0) then
         diag_datas = 0.0  ! partial sum
         diag_eday0 = eday ! 1st  day in partial sum
         diag_eday1 = -1   ! last day in partial sum
         diag_ns    =  0   ! number of samples in partial sum
      end if

   endif

   first_call = .false.

END subroutine diag_dodiag

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_atm - compute atmosphere diagnostics
!
! !DESCRIPTION:
!    Compute atmosphere diagnostics (instantaneous global averages)
!
! !REMARKS:
!    Area averages are relative to the entire unit sphere, area = 4*pi rad^2
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE diag_atm(bun_a2c,bun_c2a)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(cpl_bundle),intent(in ) :: bun_a2c ! atm->cpl bundle
   type(cpl_bundle),intent(in ) :: bun_c2a ! cpl->atm bundle

!EOP

   !----- local -----
   integer(IN),save          :: nloc = -1      ! size of local data array
   integer(IN)               :: n              ! generic index
   integer(IN)               :: rcode          ! return code
   real(R8),pointer          :: h(:),w(:),a(:) ! heat/water/area data for this model
   real(R8)                  :: da             ! area of one grid cell
   real(R8),allocatable,save :: area(:)        ! area of all grid cells
   logical,save              :: first_call = .true.

   integer(IN) ,save :: arainc,arainl,asnowc,asnowl ! aVect index: rain,snow
   integer(IN) ,save :: anidr,avsdr,anidf,avsdf     ! aVect index: albedoes
   integer(IN) ,save :: swvdr,swndr,swvdf,swndf     ! aVect index: shortwave rad
   integer(IN) ,save :: alwdn                       ! aVect index: longwave down
   integer(IN) ,save :: cmlwup                      ! aVect index: longwave up
   integer(IN) ,save :: cmlat,cmsen,cmevap          ! aVect index: latent,sens,evap

   !----- formats -----
   character(*),parameter :: subName = '(diag_atm) '
   character(*),parameter :: F00   = "('(diag_atm) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! one-time initializations
   !----------------------------------------------------------------------------
   if (first_call) then

      !--- get local size of atmosphere bundle ---
      nloc = cpl_mct_aVect_lsize(bun_a2c%data)

      !--- get indicies to fields in aVect ---
      avsdr  = cpl_mct_aVect_indexRA(bun_c2a%data,"Sx_avsdr") 
      anidr  = cpl_mct_aVect_indexRA(bun_c2a%data,"Sx_anidr") 
      avsdf  = cpl_mct_aVect_indexRA(bun_c2a%data,"Sx_avsdf") 
      anidf  = cpl_mct_aVect_indexRA(bun_c2a%data,"Sx_anidf") 
      cmlwup = cpl_mct_aVect_indexRA(bun_c2a%data,"Faxx_lwup") 
      cmlat  = cpl_mct_aVect_indexRA(bun_c2a%data,"Faxx_lat") 
      cmsen  = cpl_mct_aVect_indexRA(bun_c2a%data,"Faxx_sen") 
      cmevap = cpl_mct_aVect_indexRA(bun_c2a%data,"Faxx_evap") 

      swndr  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swndr") 
      swvdr  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swvdr") 
      swvdf  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swvdf") 
      swndf  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swndf") 

      arainc = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_rainc") 
      arainl = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_rainl") 
      asnowc = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_snowc") 
      asnowl = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_snowl") 

      alwdn  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_lwdn") 

      allocate( area(nloc)) 
      call cpl_mct_aVect_getRAttr(bun_a2c%dom%lGrid,"aream",area,rcode) 

   end if

   !----------------------------------------------------------------------------
   ! do calculations on local data
   !----------------------------------------------------------------------------

   !--- do global sum ---
   h => diag_datai(:,m_atm,p_heat )
   w => diag_datai(:,m_atm,p_water)
   a => diag_datai(:,m_atm,p_area )
   h = 0.0
   w = 0.0
   a = 0.0
   do n=1,nloc
      !--- area -----------------------
      da = area(n)
      a(n_area  ) = a(n_area  ) + da 
      !--- heat flux ------------------
      h(n_hfrz  ) =  0.0
      h(n_hmelt ) =  0.0
      h(n_hswnet) = h(n_hswnet) - da *                        &
      &  ( (1.0- bun_c2a%data%rAttr(avsdr,n))*bun_a2c%data%rAttr(swvdr ,n) &
      &  + (1.0- bun_c2a%data%rAttr(anidr,n))*bun_a2c%data%rAttr(swndr ,n) &
      &  + (1.0- bun_c2a%data%rAttr(avsdf,n))*bun_a2c%data%rAttr(swvdf ,n) &
      &  + (1.0- bun_c2a%data%rAttr(anidf,n))*bun_a2c%data%rAttr(swndf ,n) )
      h(n_hlwdn ) = h(n_hlwdn ) - da * bun_a2c%data%rAttr(alwdn ,n)
      h(n_hlwup ) = h(n_hlwup ) - da * bun_c2a%data%rAttr(cmlwup,n)
      h(n_hlat  ) = h(n_hlat  ) - da * bun_c2a%data%rAttr(cmlat ,n)
      h(n_hsen  ) = h(n_hsen  ) - da * bun_c2a%data%rAttr(cmsen ,n)
      !--- water flux -----------------
      w(n_wfrz  ) =  0.0
      w(n_wmelt ) =  0.0
      w(n_wrain ) = w(n_wrain ) - da * bun_a2c%data%rAttr(arainc,n)
      w(n_wrain ) = w(n_wrain ) - da * bun_a2c%data%rAttr(arainl,n)
      w(n_wsnow ) = w(n_wsnow ) - da * bun_a2c%data%rAttr(asnowc,n)
      w(n_wsnow ) = w(n_wsnow ) - da * bun_a2c%data%rAttr(asnowl,n)
      w(n_wevap ) = w(n_wevap ) - da * bun_c2a%data%rAttr(cmevap,n)
      w(n_wroff ) =  0.0
      !--------------------------------
   enddo

   !--- sum over sources/sinks ---
   h(n_hnet) = sum(h)
   w(n_wnet) = sum(w)

   swnet_atm = h(n_hswnet)
 
   first_call = .false.

END subroutine diag_atm

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_lnd - compute land diagnostics
!
! !DESCRIPTION:
!    Compute land diagnostics (instantaneous global averages)
!
! !REMARKS:
!    Area averages are relative to the entire unit sphere, area = 4*pi rad^2
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE diag_lnd(bun_l2c,bun_c2l,bun_r2c,bun_lfrac)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(cpl_bundle),intent(in) :: bun_l2c   ! lnd->cpl bundle
   type(cpl_bundle),intent(in) :: bun_c2l   ! cpl->lnd bundle
   type(cpl_bundle),intent(in) :: bun_r2c   ! rof->cpl bundle
   type(cpl_bundle),intent(in) :: bun_lfrac ! surface fractions on lnd domain

!EOP

   !----- local -----
   integer(IN),save    :: nloc  = -1          ! size of local data array
   integer(IN),save    :: nlocr = -1          ! size of local data array
   logical,save        :: first_call = .true. ! flags 1st time init
   integer(IN)         :: rcode               ! return code
   integer(IN)         :: n                   ! generic index
   real(R8)   ,pointer :: h(:),w(:),a(:)      ! heat/water/area data for this model
   real(R8)            :: dal                 ! area of one grid cell
   real(R8)   ,allocatable,save ::   area (:) ! area of all grid cells
   real(R8)   ,allocatable      ::   mask (:) ! domain mask array
   integer(IN),allocatable,save ::   imask(:) ! domain mask array
   real(R8)   ,allocatable,save ::  area_r(:) !
   real(R8)   ,allocatable      ::  mask_r(:) !
   integer(IN),allocatable,save :: imask_r(:) !

   integer(IN),save :: arainc,arainl,asnowc,asnowl ! aVect index: rain, snow
   integer(IN),save :: anidr,avsdr,anidf,avsdf     ! aVect index: albedoes
   integer(IN),save :: swvdr,swndr,swvdf,swndf     ! aVect index: shortwave down
   integer(IN),save :: alwdn                       ! aVect index: longwave down
   integer(IN),save :: llwup                       ! aVect index: longwave up
   integer(IN),save :: llat,lsen,levap             ! aVect index: lat,sen,evap
   integer(IN),save :: lfrac                       ! aVect index: land fraction
   integer(IN),save :: rroff                       ! aVect index: runoff

   !----- formats -----
   character(*),parameter :: subName = '(diag_lnd) '
   character(*),parameter :: F00   = "('(diag_lnd) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! one-time initializations
   !----------------------------------------------------------------------------
   if (first_call) then
      !--- get local size of land bundle ---
      nloc  = cpl_mct_aVect_lsize(bun_l2c%data)

      !--- get local size of river bundle ---
      nlocr = cpl_mct_aVect_lsize(bun_r2c%data)

      !--- get indicies to fields in bundles ---
      avsdr  = cpl_mct_aVect_indexRA(bun_l2c%data,"Sl_avsdr") 
      anidr  = cpl_mct_aVect_indexRA(bun_l2c%data,"Sl_anidr") 
      avsdf  = cpl_mct_aVect_indexRA(bun_l2c%data,"Sl_avsdf") 
      anidf  = cpl_mct_aVect_indexRA(bun_l2c%data,"Sl_anidf") 

      swndr  = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_swndr") 
      swvdr  = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_swvdr") 
      swvdf  = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_swvdf") 
      swndf  = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_swndf") 

      arainc = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_rainc") 
      arainl = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_rainl") 
      asnowc = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_snowc") 
      asnowl = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_snowl") 

      alwdn  = cpl_mct_aVect_indexRA(bun_c2l%data,"Faxa_lwdn") 

      llwup  = cpl_mct_aVect_indexRA(bun_l2c%data,"Fall_lwup") 
      llat   = cpl_mct_aVect_indexRA(bun_l2c%data,"Fall_lat" ) 
      lsen   = cpl_mct_aVect_indexRA(bun_l2c%data,"Fall_sen" ) 
      levap  = cpl_mct_aVect_indexRA(bun_l2c%data,"Fall_evap") 

      rroff  = cpl_mct_aVect_indexRA(bun_r2c%data,"Forr_roff") 

      lfrac  = cpl_mct_aVect_indexRA(bun_lfrac%data,"lfrac")

      !--- allocate local arrays ---
      allocate( area  (nloc ),mask  (nloc ),imask  (nloc )) 
      allocate( area_r(nlocr),mask_r(nlocr),imask_r(nlocr)) 

      !--- get the cell areas and domain masks data ---
      call cpl_mct_aVect_getRAttr(bun_l2c%dom%lGrid,"aream",area  ,rcode)
      call cpl_mct_aVect_getRAttr(bun_r2c%dom%lGrid,"aream",area_r,rcode)
      call cpl_mct_aVect_getRAttr(bun_l2c%dom%lGrid,"mask",mask  ,rcode)
      call cpl_mct_aVect_getRAttr(bun_r2c%dom%lGrid,"mask",mask_r,rcode)

      !--- create masks with integer values ---
      imask  (:) = 0 ; where (mask   /= 0.0) imask   = 1
      imask_r(:) = 0 ; where (mask_r /= 0.0) imask_r = 1
      deallocate(mask,mask_r)

   end if

   !----------------------------------------------------------------------------
   ! do calculations on local data
   !----------------------------------------------------------------------------
   !--- do global sum ---
   h => diag_datai(:,m_lnd,p_heat )
   w => diag_datai(:,m_lnd,p_water)
   a => diag_datai(:,m_lnd,p_area )
   h = 0.0
   w = 0.0
   a = 0.0
   do n=1,nloc
      if (imask(n) /= 0) then
         !--- area -----------------------
         dal = area(n)*bun_lfrac%data%rAttr(lfrac,n)
         a(n_area  ) = a(n_area  ) + dal
         !--- heat flux ------------------
         h(n_hfrz  ) =  0.0
         h(n_hmelt ) =  0.0
         h(n_hswnet) = h(n_hswnet) + dal *                       &
         &   ( (1.0-bun_l2c%data%rAttr(avsdr,n))*bun_c2l%data%rAttr(swvdr,n) &
         &   + (1.0-bun_l2c%data%rAttr(anidr,n))*bun_c2l%data%rAttr(swndr,n) &
         &   + (1.0-bun_l2c%data%rAttr(avsdf,n))*bun_c2l%data%rAttr(swvdf,n) &
         &   + (1.0-bun_l2c%data%rAttr(anidf,n))*bun_c2l%data%rAttr(swndf,n) )
         h(n_hlwdn ) = h(n_hlwdn ) + dal * bun_c2l%data%rAttr(alwdn,n)
         h(n_hlwup ) = h(n_hlwup ) + dal * bun_l2c%data%rAttr(llwup,n)
         h(n_hlat  ) = h(n_hlat  ) + dal * bun_l2c%data%rAttr(llat ,n)
         h(n_hsen  ) = h(n_hsen  ) + dal * bun_l2c%data%rAttr(lsen ,n)
         !--- water flux -----------------
         w(n_wfrz  ) =  0.0
         w(n_wmelt ) =  0.0
         w(n_wrain ) = w(n_wrain ) + dal * bun_c2l%data%rAttr(arainc,n)
         w(n_wrain ) = w(n_wrain ) + dal * bun_c2l%data%rAttr(arainl,n)
         w(n_wsnow ) = w(n_wsnow ) + dal * bun_c2l%data%rAttr(asnowc,n)
         w(n_wsnow ) = w(n_wsnow ) + dal * bun_c2l%data%rAttr(asnowl,n)
         w(n_wevap ) = w(n_wevap ) + dal * bun_l2c%data%rAttr(levap ,n)
         !--------------------------------
      endif
   enddo

   !--- do river model data ---
   do n=1,nlocr
      if (imask_r(n) /= 0) then
         w(n_wroff) = w(n_wroff) - area_r(n)*bun_r2c%data%rAttr(rroff,n)
      endif
   enddo

   !--- form sum across all heat/water sources/sinks
   h(n_hnet) = sum(h)
   w(n_wnet) = sum(w)

   !--- save swnet
   swnet_lnd = h(n_hswnet)

   first_call = .false.

END subroutine diag_lnd

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_ice - compute atmosphere diagnostics
!
! !DESCRIPTION:
!    Compute ice model diagnostics (instantaneous global averages)
!
! !REMARKS:
!    Area averages are relative to the entire unit sphere, area = 4*pi rad^2
!    This routine assumes ice and ocean are on the same grid.
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE diag_ice(bun_i2c,bun_c2i,bun_o2c,bun_ifrac)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(cpl_bundle),intent(in ) :: bun_i2c   ! ice->cpl bundle
   type(cpl_bundle),intent(in ) :: bun_c2i   ! cpl->ice bundle
   type(cpl_bundle),intent(in ) :: bun_o2c   ! ocn->cpl bundle
   type(cpl_bundle),intent(in ) :: bun_ifrac ! surface fractions on ice domain

!EOP

   !----- local -----
   integer(IN),save     :: nloc = -1           ! size of local data array
   logical,save         :: first_call = .true. ! flags 1-time initialization
   integer(IN)          :: rcode               ! return code
   integer(IN)          :: n                   ! generic index
   real(R8)   ,pointer  :: hn(:),wn(:),an(:)   ! heat/water/area  data for n-hemi
   real(R8)   ,pointer  :: hs(:),ws(:),as(:)   ! heat/water/area  data for s-hemi
   real(R8)             :: da,dai              ! area of one grid cell
   real(R8)             :: heatf               ! heat from freezing xor melt pot
   real(R8)   ,allocatable,save ::  area(:)    ! cell area data
   real(R8)   ,allocatable,save ::  lats(:)    ! cell latitude
   real(R8)   ,allocatable      ::  mask(:)    ! domain mask
   integer(IN),allocatable,save :: imask(:)    ! domain mask
   integer(IN),save :: melth                   ! aVect index: melt heat
   integer(IN),save :: aswnet                  ! aVect index: shortwave, net 
   integer(IN),save :: swvdr,swndr,swvdf,swndf ! aVect index: shortwave downward
   integer(IN),save :: avsdr,anidr,avsdf,anidf ! aVect index: albedoes
   integer(IN),save :: pswnet                  ! aVect index: penetrating shortwave
   integer(IN),save :: lwdn                    ! aVect index: longwave down
   integer(IN),save :: lwup                    ! aVect index: longwave up
   integer(IN),save :: lat,sen                 ! aVect index: lat,sen
   integer(IN),save :: meltw                   ! aVect index: meltwater
   integer(IN),save :: rain,snow,evap          ! aVect index: rain,snow,evap
   integer(IN),save :: oq                      ! aVect index: heat from melting ice
   integer(IN),save :: ifrac                   ! aVect index: ice fraction

   !----- formats -----
   character(*),parameter :: subName = '(diag_ice) '
   character(*),parameter :: F00   = "('(diag_ice) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! one-time initializations
   !----------------------------------------------------------------------------
   if (first_call) then

      !--- get local size of atmosphere bundle ---
      nloc = cpl_mct_aVect_lsize(bun_i2c%data)

      avsdr  = cpl_mct_aVect_indexRA(bun_i2c%data,"Si_avsdr") 
      anidr  = cpl_mct_aVect_indexRA(bun_i2c%data,"Si_anidr") 
      avsdf  = cpl_mct_aVect_indexRA(bun_i2c%data,"Si_avsdf") 
      anidf  = cpl_mct_aVect_indexRA(bun_i2c%data,"Si_anidf") 

      swndr  = cpl_mct_aVect_indexRA(bun_c2i%data,"Faxa_swndr") 
      swvdr  = cpl_mct_aVect_indexRA(bun_c2i%data,"Faxa_swvdr") 
      swvdf  = cpl_mct_aVect_indexRA(bun_c2i%data,"Faxa_swvdf") 
      swndf  = cpl_mct_aVect_indexRA(bun_c2i%data,"Faxa_swndf") 
      lwdn   = cpl_mct_aVect_indexRA(bun_c2i%data,"Faxa_lwdn" ) 

      melth  = cpl_mct_aVect_indexRA(bun_i2c%data,"Fioi_melth") 
      aswnet = cpl_mct_aVect_indexRA(bun_i2c%data,"Faii_swnet") 
      pswnet = cpl_mct_aVect_indexRA(bun_i2c%data,"Fioi_swpen") 
      lwup   = cpl_mct_aVect_indexRA(bun_i2c%data,"Faii_lwup" ) 
      lat    = cpl_mct_aVect_indexRA(bun_i2c%data,"Faii_lat"  ) 
      sen    = cpl_mct_aVect_indexRA(bun_i2c%data,"Faii_sen"  ) 

      meltw  = cpl_mct_aVect_indexRA(bun_i2c%data,"Fioi_meltw") 
      rain   = cpl_mct_aVect_indexRA(bun_c2i%data,"Faxc_rain") 
      snow   = cpl_mct_aVect_indexRA(bun_c2i%data,"Faxc_snow") 
      evap   = cpl_mct_aVect_indexRA(bun_i2c%data,"Faii_evap" ) 
      
      ifrac  = cpl_mct_aVect_indexRA(bun_ifrac%data,"ifrac") 

      oq     = cpl_mct_aVect_indexRA(bun_o2c%data,"Fioo_q") 

      allocate( area(nloc),mask(nloc),imask(nloc),lats(nloc))
      call cpl_mct_aVect_getRAttr(bun_i2c%dom%lGrid,"aream",area,rcode)
      call cpl_mct_aVect_getRAttr(bun_i2c%dom%lGrid,"lat" ,lats,rcode)
      call cpl_mct_aVect_getRAttr(bun_i2c%dom%lGrid,"mask",mask,rcode)
      imask(:) = 0 ; where (mask /= 0.0) imask = 1
      deallocate(mask)

   end if

   !----------------------------------------------------------------------------
   ! do calculations on local data
   !----------------------------------------------------------------------------
 
   !--- do global sum ---
   hn => diag_datai(:,m_ice_nh,p_heat )
   wn => diag_datai(:,m_ice_nh,p_water)
   an => diag_datai(:,m_ice_nh,p_area )

   hs => diag_datai(:,m_ice_sh,p_heat )
   ws => diag_datai(:,m_ice_sh,p_water)
   as => diag_datai(:,m_ice_sh,p_area )
 
   hn = 0.0
   wn = 0.0
   an = 0.0
 
   hs = 0.0
   ws = 0.0
   as = 0.0
 
   swnet_ice_nh = 0.0 ! need swnet w/o penetrating sw for solar verification
   swnet_ice_sh = 0.0 
 
   do n=1,nloc

      da  = area(n)
      dai = da*bun_ifrac%data%rAttr(ifrac,n)
   
      if (imask(n) /= 0 .and. lats(n) > 0.0) then ! NORTHERN HEMISPHERE
         an(n_area  ) = an(n_area  ) + dai
         !--- heat flux ------------------
         heatf=bun_o2c%data%rAttr(oq,n)
         hn(n_hfrz  ) = hn(n_hfrz  ) - da *max(heatf,0.0)
         hn(n_hmelt ) = hn(n_hmelt ) - dai*bun_i2c%data%rAttr(melth,n)
         swnet_ice_nh = swnet_ice_nh + dai*                           &
         & ( (1.0-bun_i2c%data%rAttr(avsdr,n))*bun_c2i%data%rAttr(swvdr ,n) &
         & + (1.0-bun_i2c%data%rAttr(anidr,n))*bun_c2i%data%rAttr(swndr ,n) &
         & + (1.0-bun_i2c%data%rAttr(avsdf,n))*bun_c2i%data%rAttr(swvdf ,n) &
         & + (1.0-bun_i2c%data%rAttr(anidf,n))*bun_c2i%data%rAttr(swndf ,n) )
         hn(n_hswnet)  = hn(n_hswnet) - dai*bun_i2c%data%rAttr(pswnet,n)
         hn(n_hlwdn )  = hn(n_hlwdn ) + dai*bun_c2i%data%rAttr(lwdn  ,n)
         hn(n_hlwup )  = hn(n_hlwup ) + dai*bun_i2c%data%rAttr(lwup  ,n)
         hn(n_hlat  )  = hn(n_hlat  ) + dai*bun_i2c%data%rAttr(lat   ,n)
         hn(n_hsen  )  = hn(n_hsen  ) + dai*bun_i2c%data%rAttr(sen   ,n)
         !--- water flux -----------------
!        wn(n_wfrz  )  = hn(n_hfrz  )*HFLXtoWFLX ! implied water flux
         wn(n_wmelt )  = wn(n_wmelt ) - dai*bun_i2c%data%rAttr(meltw,n)
         wn(n_wrain )  = wn(n_wrain ) + dai*bun_c2i%data%rAttr(rain ,n)
         wn(n_wsnow )  = wn(n_wsnow ) + dai*bun_c2i%data%rAttr(snow ,n)
         wn(n_wevap )  = wn(n_wevap ) + dai*bun_i2c%data%rAttr(evap ,n)
         wn(n_wroff )  =  0.0
   
      else if (imask(n) /= 0 .and. lats(n) < 0.0) then ! SOUTHERN HEMISPHERE
         as(n_area  ) = as(n_area  ) + dai
         !--- heat flux ------------------
         heatf=bun_o2c%data%rAttr(oq,n)
         hs(n_hfrz  ) = hs(n_hfrz  ) - da *max(heatf,0.0)
         hs(n_hmelt ) = hs(n_hmelt ) - dai*bun_i2c%data%rAttr(melth,n)
         swnet_ice_sh = swnet_ice_sh + dai*                           &
         & ( (1.0-bun_i2c%data%rAttr(avsdr,n))*bun_c2i%data%rAttr(swvdr ,n) &
         & + (1.0-bun_i2c%data%rAttr(anidr,n))*bun_c2i%data%rAttr(swndr ,n) &
         & + (1.0-bun_i2c%data%rAttr(avsdf,n))*bun_c2i%data%rAttr(swvdf ,n) &
         & + (1.0-bun_i2c%data%rAttr(anidf,n))*bun_c2i%data%rAttr(swndf ,n) )
         hs(n_hswnet) = hs(n_hswnet) - dai*bun_i2c%data%rAttr(pswnet,n)
         hs(n_hlwdn ) = hs(n_hlwdn ) + dai*bun_c2i%data%rAttr(lwdn  ,n)
         hs(n_hlwup ) = hs(n_hlwup ) + dai*bun_i2c%data%rAttr(lwup  ,n)
         hs(n_hlat  ) = hs(n_hlat  ) + dai*bun_i2c%data%rAttr(lat   ,n)
         hs(n_hsen  ) = hs(n_hsen  ) + dai*bun_i2c%data%rAttr(sen   ,n)
         !--- water flux -----------------
!        ws(n_wfrz  ) = hs(n_hfrz  )*HFLXtoWFLX ! implied water flux
         ws(n_wmelt ) = ws(n_wmelt ) - dai*bun_i2c%data%rAttr(meltw,n)
         ws(n_wrain ) = ws(n_wrain ) + dai*bun_c2i%data%rAttr(rain ,n)
         ws(n_wsnow ) = ws(n_wsnow ) + dai*bun_c2i%data%rAttr(snow ,n)
         ws(n_wevap ) = ws(n_wevap ) + dai*bun_i2c%data%rAttr(evap ,n)
         ws(n_wroff ) =  0.0
      endif
   enddo
         ws(n_wfrz  ) = hs(n_hfrz  )*HFLXtoWFLX ! implied water flux
         wn(n_wfrz  ) = hn(n_hfrz  )*HFLXtoWFLX ! implied water flux
   hn(n_hswnet) = hn(n_hswnet) + swnet_ice_nh
   hs(n_hswnet) = hs(n_hswnet) + swnet_ice_sh
 
   hn(n_hnet) = sum(hn)
   wn(n_wnet) = sum(wn)

   hs(n_hnet) = sum(hs)
   ws(n_wnet) = sum(ws)
 
   first_call = .false.

END subroutine diag_ice

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_ocn - compute ocean diagnostics
!
! !DESCRIPTION:
!    Compute ocean diagnostics (instantaneous global averages)
!
! !REMARKS:
!    Area averages are relative to the entire unit sphere, area = 4*pi rad^2
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE diag_ocn(bun_o2c,bun_c2o,bun_a2c,bun_i2c,bun_alb,bun_ofrac)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(cpl_bundle),intent(in ) :: bun_o2c    ! ocn->cpl bundle
   type(cpl_bundle),intent(in ) :: bun_c2o    ! cpl->ocn bundle
   type(cpl_bundle),intent(in ) :: bun_a2c    ! atm->cpl bundle
   type(cpl_bundle),intent(in ) :: bun_i2c    ! ice->cpl bundle
   type(cpl_bundle),intent(in ) :: bun_alb    ! albedo   bundle
   type(cpl_bundle),intent(in ) :: bun_ofrac  ! surface fractions on ocn domain

!EOP

   !----- local -----
   integer(IN),save    :: nloc = -1           ! size of local data array
   logical,save        :: first_call = .true. ! flags 1st time initialization
   integer(IN)         :: rcode               ! return code
   integer(IN)         :: n                   ! generic index
   real(R8)   ,pointer :: h(:),w(:),a(:)      ! heat/water/area data this mode only 
   real(R8)            :: da,dai,daa          ! area of one grid cell
   real(R8)            :: heatf               ! heat from freezing xor melt pot

   real(R8),allocatable,save ::  area(:)      ! cell area
   real(R8),allocatable,save :: imask(:)      ! domain mask
   real(R8),allocatable      ::  mask(:)      ! domain mask

   integer(IN) ,save :: anidr,avsdr,anidf,avsdf ! aVect index: albedoes
   integer(IN) ,save :: swvdr,swndr,swvdf,swndf ! aVect index: shortwave down
   integer(IN) ,save :: rain,snow               ! aVect index: rain, snow
   integer(IN) ,save :: melth                   ! aVect index: heat from melting
   integer(IN) ,save :: lwdn,lwup               ! aVect index: longwave down, up
   integer(IN) ,save :: lat,sen                 ! aVect index: latent, sensible
   integer(IN) ,save :: meltw,evap,rroff        ! aVect index: meltwater, evap, roff
   integer(IN) ,save :: pswnet,oq               ! aVect index: pen sw, melt heat
   integer(IN) ,save :: ifrac,afrac             ! aVect index: ice & atm fraction

   !----- formats -----
   character(*),parameter :: subName = '(diag_ocn) '
   character(*),parameter :: F00   = "('(diag_ocn) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! one-time initializations
   !----------------------------------------------------------------------------
   if (first_call) then

      !--- get local size of atmosphere bundle ---
      nloc = cpl_mct_aVect_lsize(bun_o2c%data)

      swndr  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swndr") 
      swvdr  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swvdr") 
      swvdf  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swvdf") 
      swndf  = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swndf") 

      avsdr  = cpl_mct_aVect_indexRA(bun_alb%data,"So_avsdr") 
      anidr  = cpl_mct_aVect_indexRA(bun_alb%data,"So_anidr") 
      avsdf  = cpl_mct_aVect_indexRA(bun_alb%data,"So_avsdf") 
      anidf  = cpl_mct_aVect_indexRA(bun_alb%data,"So_anidf") 

      melth  = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_melth") 
      lwdn   = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_lwdn" ) 
      lwup   = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_lwup") 
      lat    = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_lat" ) 
      sen    = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_sen" ) 
      meltw  = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_meltw") 
      evap   = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_evap") 
      rain   = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_rain") 
      snow   = cpl_mct_aVect_indexRA(bun_c2o%data,"Foxx_snow") 
      rroff  = cpl_mct_aVect_indexRA(bun_c2o%data,"Forr_roff") 
      pswnet = cpl_mct_aVect_indexRA(bun_i2c%data,"Fioi_swpen") 
      oq     = cpl_mct_aVect_indexRA(bun_o2c%data,"Fioo_q"    ) 

      ifrac  = cpl_mct_aVect_indexRA(bun_ofrac%data,"ifrac") 
      afrac  = cpl_mct_aVect_indexRA(bun_ofrac%data,"afrac") 

      !--- form global area, mask & lat data ---
      allocate( area(nloc),mask(nloc),imask(nloc)) 
      call cpl_mct_aVect_getRAttr(bun_o2c%dom%lGrid,"aream",area,rcode)
      call cpl_mct_aVect_getRAttr(bun_o2c%dom%lGrid,"mask",mask,rcode)
      imask(:) = 0 ; where (mask /= 0.0) imask = 1
      deallocate(mask)

   end if

   !----------------------------------------------------------------------------
   ! do calculations on local data
   !----------------------------------------------------------------------------
   !--- do global sum ---
   h => diag_datai(:,m_ocn,p_heat )
   w => diag_datai(:,m_ocn,p_water)
   a => diag_datai(:,m_ocn,p_area )
   h = 0.0
   w = 0.0
   a = 0.0
   swnet_ocn = 0.0 ! need swnet w/o penetrating sw for solar verificatio
   do n=1,nloc
      if (imask(n) /= 0) then
         da  = area(n)
         daa = da*bun_ofrac%data%rAttr(afrac,n)
         dai = da*bun_ofrac%data%rAttr(ifrac,n)
         !--- area -----------------------
         a(n_area  ) = a(n_area  ) + daa
         !--- heat flux ------------------
         heatf       =                   bun_o2c%data%rAttr(oq,n)
         h(n_hfrz  ) = h(n_hfrz  ) + da *max(heatf,0.0)
         h(n_hmelt ) = h(n_hmelt ) + da *bun_c2o%data%rAttr(melth ,n)
         h(n_hswnet) = h(n_hswnet) + dai*bun_i2c%data%rAttr(pswnet,n)
         swnet_ocn = swnet_ocn     + daa*                      &
         & ( (1.0-bun_alb%data%rAttr(avsdr,n))*bun_a2c%data%rAttr(swvdr,n) &
         & + (1.0-bun_alb%data%rAttr(anidr,n))*bun_a2c%data%rAttr(swndr,n) &
         & + (1.0-bun_alb%data%rAttr(avsdf,n))*bun_a2c%data%rAttr(swvdf,n) &
         & + (1.0-bun_alb%data%rAttr(anidf,n))*bun_a2c%data%rAttr(swndf,n) )
         h(n_hlwdn ) = h(n_hlwdn ) + da*bun_c2o%data%rAttr(lwdn ,n)
         h(n_hlwup ) = h(n_hlwup ) + da*bun_c2o%data%rAttr(lwup ,n)
         h(n_hlat  ) = h(n_hlat  ) + da*bun_c2o%data%rAttr(lat  ,n)
         h(n_hsen  ) = h(n_hsen  ) + da*bun_c2o%data%rAttr(sen  ,n)
         !--- water flux -----------------
!        w(n_wfrz  ) = h(n_hfrz  )*HFLXtoWFLX ! implied water flux
         w(n_wmelt ) = w(n_wmelt ) + da*bun_c2o%data%rAttr(meltw,n)
         w(n_wrain ) = w(n_wrain ) + da*bun_c2o%data%rAttr(rain ,n)
         w(n_wsnow ) = w(n_wsnow ) + da*bun_c2o%data%rAttr(snow ,n)
         w(n_wevap ) = w(n_wevap ) + da*bun_c2o%data%rAttr(evap ,n)
         w(n_wroff ) = w(n_wroff ) + da*bun_c2o%data%rAttr(rroff,n)
         !--------------------------------
      endif
   enddo
         w(n_wfrz )  = h(n_hfrz )*HFLXtoWFLX ! implied water flux
   h(n_hswnet) = h(n_hswnet) + swnet_ocn

   h(n_hnet) = sum(h)
   w(n_wnet) = sum(w)
 
   first_call = .false.

END subroutine diag_ocn

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_print - print out diagnostics
!
! !DESCRIPTION:
!   Print the diagnostics and their sum in each category
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE diag_print(date)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(shr_date),intent(in)  :: date ! current model date

!EOP
   !--- local ---
   integer(IN) :: m           ! data array indicies
   integer(IN) :: nday        ! number of days in time avg
   real(R8)    :: area (  6)  ! area data
   real(R8)    :: heat (8,6)  ! heat data 
   real(R8)    :: water(8,6)  ! water data 
   integer(IN) :: cdate,sec   ! coded date, seconds

   !----- formats -----
   character(*),parameter :: subName = '(diag_print) '
   character(*),parameter :: F00   = "('(diag_print) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !---------------------------------------------------------------
   ! instantaneous data
   !---------------------------------------------------------------
   area  = diag_Gdatai(n_area,:,p_area )
   heat  = diag_Gdatai(   :  ,:,p_heat )
   water = diag_Gdatai(   :  ,:,p_water)

   water = water*1.0e6 ! scale data for printout

   call shr_date_getCDate(date,cdate,sec)

   write(6,F10)'heat  inst       '
   write(6,F11)'heat  inst       '
   write(6,F12)'heat  inst atm   ',cdate, sec,area(m_atm   ) ,heat(1:8,m_atm   )
   write(6,F12)'heat  inst lnd   ',cdate, sec,area(m_lnd   ) ,heat(1:8,m_lnd   )
   write(6,F12)'heat  inst ice_nh',cdate, sec,area(m_ice_nh) ,heat(1:8,m_ice_nh)
   write(6,F12)'heat  inst ice_sh',cdate, sec,area(m_ice_sh) ,heat(1:8,m_ice_sh)
   write(6,F12)'heat  inst ocn   ',cdate, sec,area(m_ocn   ) ,heat(1:8,m_ocn   )
   write(6,F11)'heat  inst       '
   write(6,F12)'heat  inst sum   ',cdate, sec,area(m_sum   ) ,heat(1:8,m_sum   )
   write(6,*)

   write(6,F20)'water inst       '
   write(6,F21)'water inst       '
   write(6,F22)'water inst atm   ',cdate, sec,area(m_atm   ),water(1:7,m_atm   )
   write(6,F22)'water inst lnd   ',cdate, sec,area(m_lnd   ),water(1:7,m_lnd   )
   write(6,F22)'water inst ice_nh',cdate, sec,area(m_ice_nh),water(1:7,m_ice_nh)
   write(6,F22)'water inst ice_sh',cdate, sec,area(m_ice_sh),water(1:7,m_ice_sh)
   write(6,F22)'water inst ocn   ',cdate, sec,area(m_ocn   ),water(1:7,m_ocn   )
   write(6,F21)'water inst       '
   write(6,F22)'water inst sum   ',cdate, sec,area(m_sum   ),water(1:7,m_sum   )
   write(6,*)

END SUBROUTINE diag_print

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_printAvg - print out diagnostics for time-avg data
!
! !DESCRIPTION:
!   Print the diagnostics and their sum in each category
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE diag_printAvg(date)

! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(shr_date),intent(in)  :: date ! current model date

!EOP

   !--- local ---
   integer(IN) :: m           ! data array indicies
   integer(IN) :: nday        ! number of days in time avg
   real(R8)    :: area (  6)  ! area data
   real(R8)    :: heat (8,6)  ! heat data 
   real(R8)    :: water(8,6)  ! water data 
   integer(IN) :: cdate,sec,eday

   !----- formats -----
   character(*),parameter :: subName = '(diag_printAvg) '
   character(*),parameter :: F00   = "('(diag_printAvg) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call shr_date_getCDate(date,cdate,sec)

   !--- return if there is no data to print out ---
   if (diag_ns < 1) then
      write(6,F00) 'NO SAMPLES IN TIME AVERAGE DATA'
      write(6,F12) 'heat  tavg nodata',cdate
      write(6,F22) 'water tavg nodata',cdate
      return
   end if

   !---------------------------------------------------------------
   ! average data
   !---------------------------------------------------------------
   area  = diag_datas(n_area,:,p_area )
   heat  = diag_datas(   :  ,:,p_heat )
   water = diag_datas(   :  ,:,p_water)

   water = water*1.0e6 ! scale data for printout

   area  =  area/diag_ns
   heat  =  heat/diag_ns
   water = water/diag_ns
   nday  = nint(diag_eday1 - diag_eday0)

   write(6,F30)'heat  tavg       '
   write(6,F11)'heat  tavg       '
   write(6,F12)'heat  tavg atm   ',cdate,nday,area(m_atm   ) ,heat(1:8,m_atm   )
   write(6,F12)'heat  tavg lnd   ',cdate,nday,area(m_lnd   ) ,heat(1:8,m_lnd   )
   write(6,F12)'heat  tavg ice_nh',cdate,nday,area(m_ice_nh) ,heat(1:8,m_ice_nh)
   write(6,F12)'heat  tavg ice_sh',cdate,nday,area(m_ice_sh) ,heat(1:8,m_ice_sh)
   write(6,F12)'heat  tavg ocn   ',cdate,nday,area(m_ocn   ) ,heat(1:8,m_ocn   )
   write(6,F11)'heat  tavg       '
   write(6,F12)'heat  tavg sum   ',cdate,nday,area(m_sum   ) ,heat(1:8,m_sum   )
   write(6,*)

   write(6,F40)'water tavg       '
   write(6,F21)'water tavg       '
   write(6,F22)'water tavg atm   ',cdate,nday,area(m_atm   ),water(1:7,m_atm   )
   write(6,F22)'water tavg lnd   ',cdate,nday,area(m_lnd   ),water(1:7,m_lnd   )
   write(6,F22)'water tavg ice_nh',cdate,nday,area(m_ice_nh),water(1:7,m_ice_nh)
   write(6,F22)'water tavg ice_sh',cdate,nday,area(m_ice_sh),water(1:7,m_ice_sh)
   write(6,F22)'water tavg ocn   ',cdate,nday,area(m_ocn   ),water(1:7,m_ocn   )
   write(6,F21)'water tavg       '
   write(6,F22)'water tavg sum   ',cdate,nday,area(m_sum   ),water(1:7,m_sum   )
   write(6,*)

   call shr_sys_flush(6)

END SUBROUTINE diag_printAvg

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: diag_solar - compares expected vs. actual short-wave radiation
!
! !DESCRIPTION:
!    Compare expected vs. actual short-wave net (absorbed solar)
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE diag_solar(bun_a2c,bun_l2c,bun_i2c,bun_lfrac,bun_ifrac)

! !INPUT PARAMETERS:

   implicit none

   type(cpl_bundle),intent(in) :: bun_a2c   !  atm->cpl bundle
   type(cpl_bundle),intent(in) :: bun_l2c   !  lnd->cpl bundleand data
   type(cpl_bundle),intent(in) :: bun_i2c   !  ice->cpl bundlece data
   type(cpl_bundle),intent(in) :: bun_lfrac !  surface fractions on lnd domain
   type(cpl_bundle),intent(in) :: bun_ifrac !  surface fractions on ice domain

!EOP

   !----- local -----
   integer(IN) :: n                     ! generic loop index
   integer(IN) :: rcode                 ! return code
   integer(IN),save :: nloca            ! local atm points
   integer(IN),save :: nloci            ! local ice points
   integer(IN),save :: nlocl            ! local ocean points
   real(R8)    :: sa1,sl1,sin1,sis1,so1 ! net solar computed in cpl diagnostics
   real(R8)    :: sa2,sl2,sin2,sis2,so2 ! net solar sent to cpl from comps, local sum
   real(R8)    :: sa2g,sl2g,sin2g,sis2g ! net solar sent to cpl from comps, global

   logical,save                 :: first_call = .true.
   integer(IN),save             :: aswnet
   integer(IN),save             :: lswnet
   integer(IN),save             :: iswnet
   integer(IN),save             :: lfrac,ifrac
   real(R8)   ,allocatable,save :: areai(:),areaa(:),areal(:)
   real(R8)   ,allocatable,save :: latsi(:)
   real(R8)   ,allocatable      ::  maski(:), maskl(:)
   integer(IN),allocatable,save :: imaski(:),imaskl(:)

   !----- formats -----
   character(*),parameter :: subName = '(diag_solar) '
   character(*),parameter :: F00="('(diag_solar) ',a,2f16.6)"
   character(*),parameter :: F05="('(diag_solar) ', &
   &       'swnet        atm         lnd      ice-nh      ice-sh         ocn ')"
   character(*),parameter :: F06="('(diag_solar) ', &
   &       '      ----------- ----------- ----------- ----------- -----------')"
   character(*),parameter :: F07="('(diag_solar)   cpl',5f12.6)"
   character(*),parameter :: F08="('(diag_solar) model',5f12.6)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! one-time initializations
   !----------------------------------------------------------------------------
   if (first_call) then

      nloca = cpl_mct_aVect_lsize(bun_a2c%data)
      nlocl = cpl_mct_aVect_lsize(bun_l2c%data)
      nloci = cpl_mct_aVect_lsize(bun_i2c%data)

      aswnet = cpl_mct_aVect_indexRA(bun_a2c%data,"Faxa_swnet") 

      lswnet = cpl_mct_aVect_indexRA(bun_l2c%data,"Fall_swnet") 
      lfrac = cpl_mct_aVect_indexRA(bun_lfrac%data,"lfrac") 

      iswnet = cpl_mct_aVect_indexRA(bun_i2c%data,"Faii_swnet") 
      ifrac = cpl_mct_aVect_indexRA(bun_ifrac%data,"ifrac") 

      !--- form area, mask & lat data ---
      allocate( areai(nloci),maski(nloci),imaski(nloci),stat=rcode)
      if (rcode /= 0) write(6,*) "Can't allocate ice data",rcode

      allocate(latsi(nloci),stat=rcode)
      if (rcode /= 0) write(6,*) "Can't allocate ice lats",rcode

      allocate( areaa(nloca),stat=rcode)
      if (rcode /= 0) write(6,*) "Can't allocate atmo data",rcode

      allocate( areal(nlocl),maskl(nlocl),imaskl(nlocl),stat=rcode)
      if (rcode /= 0) write(6,*)"Cant allocate land data",rcode

      call cpl_mct_aVect_getRAttr(bun_i2c%dom%lGrid,"aream",areai,rcode)
      call cpl_mct_aVect_getRAttr(bun_a2c%dom%lGrid,"aream",areaa,rcode)
      call cpl_mct_aVect_getRAttr(bun_l2c%dom%lGrid,"aream",areal,rcode)
      call cpl_mct_aVect_getRAttr(bun_i2c%dom%lGrid,"lat" ,latsi,rcode)
      call cpl_mct_aVect_getRAttr(bun_i2c%dom%lGrid,"mask",maski,rcode)
      call cpl_mct_aVect_getRAttr(bun_l2c%dom%lGrid,"mask",maskl,rcode)

      imaski(:) = 0 ; where (maski /= 0.0) imaski = 1
      imaskl(:) = 0 ; where (maskl /= 0.0) imaskl = 1
      deallocate(maski,maskl)

   endif

   !----------------------------------------------------------------------------
   ! do calculations on local data
   !----------------------------------------------------------------------------

   !--- atm --------------------------------------------------------------
   sa2 = 0.0
   do n=1,nloca
      sa2 = sa2 - areaa(n)*bun_a2c%data%rAttr(aswnet,n)
   end do

   call shr_mpi_sum(swnet_atm,sa1,cpl_comm_comp,subName//" swnet_atm")
   call shr_mpi_sum(sa2,     sa2g,cpl_comm_comp,subName//" sa2")  

   if (cpl_comm_comp_pid == 0) then
      sa1  = sa1 /(4.0*cpl_const_pi)
      sa2g = sa2g/(4.0*cpl_const_pi)
   endif
 
   !--- lnd --------------------------------------------------------------
   sl2 = 0.0
   do n=1,nlocl
      if ( imaskl(n) /= 0 ) then
         sl2 = sl2 + areal(n)*bun_lfrac%data%rAttr(lfrac,n)*bun_l2c%data%rAttr(lswnet,n)
      end if
   end do

   call shr_mpi_sum(swnet_lnd,sl1,cpl_comm_comp,subName//" swnet_lnd")
   call shr_mpi_sum(sl2,     sl2g,cpl_comm_comp,subName//" sl2")  
   if (cpl_comm_comp_pid == 0) then
      sl1  = sl1 /(4.0*cpl_const_pi)
      sl2g = sl2g/(4.0*cpl_const_pi)
   endif
 
   !--- ice-nh and ice-sh ------------------------------------------------
   sin2 = 0.0
   sis2 = 0.0
   do n=1,nloci
      if ( imaski(n) /= 0 .and. latsi(n) > 0.0) then
         sin2 = sin2 + areai(n)*bun_ifrac%data%rAttr(ifrac,n)*bun_i2c%data%rAttr(iswnet,n)
      end if
      if ( imaski(n) /= 0 .and. latsi(n) < 0.0) then
         sis2 = sis2 + areai(n)*bun_ifrac%data%rAttr(ifrac,n)*bun_i2c%data%rAttr(iswnet,n)
      end if
   end do

   call shr_mpi_sum(swnet_ice_nh,sin1,cpl_comm_comp,subName//" swnet_icen")
   call shr_mpi_sum(swnet_ice_sh,sis1,cpl_comm_comp,subName//" swnet_ices")
   call shr_mpi_sum(sin2,       sin2g,cpl_comm_comp,subName//" sa2")  
   call shr_mpi_sum(sis2,       sis2g,cpl_comm_comp,subName//" sa2")  

   if (cpl_comm_comp_pid == 0) then
      sin1  = sin1 /(4.0*cpl_const_pi)
      sis1  = sis1 /(4.0*cpl_const_pi)
      sin2g = sin2g/(4.0*cpl_const_pi)
      sis2g = sis2g/(4.0*cpl_const_pi)
   endif
 
   !--- ocn --------------------------------------------------------------
   ! note: ocn does not provide cpl with any sw-net info
   !----------------------------------------------------------------------
   call shr_mpi_sum(swnet_ocn,so1,cpl_comm_comp,subName//" swnet_ocn")
   so2 = -999.0     ! flags invalid/missing data
 
   !--- print ------------------------------------------------------------
   if (cpl_comm_comp_pid == 0) then
      so1=so1/(4.0*cpl_const_pi)
      write(6,F05)
      write(6,F06)
      write(6,F07) sa1 ,sl1 ,sin1 ,sis1 ,so1
      write(6,F08) sa2g,sl2g,sin2g,sis2g
      call shr_sys_flush(6)
   endif
 
   first_call = .false.

END SUBROUTINE diag_solar

!===============================================================================
!===============================================================================

end module diag_mod
