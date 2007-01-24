!===============================================================================
! CVS: $Id: flux_mod.F90,v 1.1.1.1 2005/02/03 22:29:00 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/flux_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: flux_mod -- does coupler's flux calculations.
!
! !DESCRIPTION:
!
!     The coupler is required to do certain flux calculations --
!     those calculations are located in this module.
!     
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

module flux_mod

! !USES:

   use shr_sys_mod     ! shared system routines
   use shr_date_mod    ! shared date module
   use shr_mpi_mod     ! shared mpi layer
   use cpl_kind_mod    ! kinds
   use cpl_const_mod   ! constants module
   use cpl_mct_mod     ! mct library
   use cpl_comm_mod    ! communicator module
   use cpl_fields_mod  ! list of fields found in bundles
   use cpl_domain_mod  ! domain data types
   use cpl_bundle_mod  ! bundle data types
   use cpl_control_mod, dbug=>cpl_control_infoDBug

   implicit none

   private ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: flux_atmOcn  ! computes atm/ocn fluxes
   public :: flux_albo    ! computes ocn albedos
   public :: flux_albi    ! modifies ice reference albedo
   public :: flux_solar   ! computes ocn net solar
   public :: flux_epbal   ! forces evap/precip/runoff balance

! !PUBLIC DATA MEMBERS:

  ! none

!EOP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: flux_atmOcn - wrapper to atm/ocn flux calculation
!
! !DESCRIPTION:
!     wrapper to atm/ocn flux calculation
!
! !REMARKS:
!     All data must be on the ocean domain (note: a domain includes a 
!     particular decomposition).
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine flux_atmOcn(bun_ocn,bun_atm,fabricate,bun_flux)

! !USES:
   use shr_timer_mod

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in ) :: bun_ocn    ! ocn state fields on ocn domain
   type(cpl_bundle),intent(in ) :: bun_atm    ! atm state fields on ocn domain
   logical         ,intent(in ) :: fabricate  ! T <=> fabriate/clobber input data
   type(cpl_bundle),intent(out) :: bun_flux   ! flux fields on ocn grid

!EOP

   !--- local ---
   integer(IN)             :: nloc = -1  ! size of local data array
   integer(IN)             :: rcode      ! return code
   integer(IN)             :: m,n        ! generic indicies
   integer(IN)             :: nchunks    ! openmp: number of chunks
   integer(IN)             :: chunksize  ! openmp: size of chunks except last one
   integer(IN)             :: chunklast  ! openmp: size of last chunk
   integer(IN)             :: cs         ! openmp: size of a particular chunck
   integer(IN)             ::t1,t2,t3,t4 ! timers
   integer(IN)             :: k1,k2,k3,k4,k5,k6,k7,k8,k9,k10 ! field indices

   real(R8)   ,allocatable ::  rmask(:)  ! ocn domain mask
   integer(IN),allocatable ::  mask (:)  ! ocn domain mask: 0 <=> inactive cell
   real(R8)   ,allocatable ::  uocn (:)  ! ocn velocity, zonal
   real(R8)   ,allocatable ::  vocn (:)  ! ocn velocity, meridional
   real(R8)   ,allocatable ::  tocn (:)  ! ocn SST
   real(R8)   ,allocatable ::  z    (:)  ! atm level height
   real(R8)   ,allocatable ::  u    (:)  ! atm velocity, zonal     
   real(R8)   ,allocatable ::  v    (:)  ! atm velocity, meridional
   real(R8)   ,allocatable ::  ptem (:)  ! atm potential T
   real(R8)   ,allocatable ::  shum (:)  ! atm specific humidity
   real(R8)   ,allocatable ::  dens (:)  ! atm density
   real(R8)   ,allocatable ::  tbot (:)  ! atm bottom surface T
   real(R8)   ,allocatable ::  sen  (:)  ! heat flux: sensible 
   real(R8)   ,allocatable ::  lat  (:)  ! heat flux: latent   
   real(R8)   ,allocatable ::  lwup (:)  ! heat flux: longwave, upward
   real(R8)   ,allocatable ::  evap (:)  ! water flux: evaporation
   real(R8)   ,allocatable ::  taux (:)  ! wind stress, zonal
   real(R8)   ,allocatable ::  tauy (:)  ! wind stress, meridional
   real(R8)   ,allocatable ::  tref (:)  ! diagnostic:  2m ref T
   real(R8)   ,allocatable ::  qref (:)  ! diagnostic:  2m ref Q
   real(R8)   ,allocatable :: duu10n(:)  ! diagnostic: 10m wind speed squared

#ifdef  UNTESTED_OPENMP
   integer,external :: omp_get_max_threads
#endif

   save

   !--- formats ---
   character(*),parameter :: subName = '(flux_atmOcn) '
   character(*),parameter :: F00 = '(  "(flux_atmOcn) ",4a  )'
   character(*),parameter :: F01 = '(  "(flux_atmOcn) ",a,i5)'

!-------------------------------------------------------------------------------
! NOTES:
! o all data is on ocn grid
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! first call allocations & setup
   !----------------------------------------------------------------------------
   if (nloc < 1 ) then
      nloc = cpl_mct_aVect_lsize(bun_ocn%data)
      call shr_timer_get(t1,'flux_ao_t1')
      call shr_timer_get(t2,'flux_ao_t2')
      call shr_timer_get(t3,'flux_ao_t3')
      call shr_timer_get(t4,'flux_ao_t4')

      !--- input fields ---
      allocate( uocn(nloc),vocn(nloc),tocn(nloc))
      allocate(rmask(nloc),mask(nloc))
      allocate(z    (nloc),u   (nloc),v   (nloc))
      allocate(ptem (nloc),shum(nloc),dens(nloc),tbot(nloc))

      !--- output fields ---
      allocate(sen  (nloc),lat (nloc),lwup(nloc),evap(nloc))
      allocate(taux (nloc),tauy(nloc),tref(nloc),qref(nloc),duu10n(nloc))

      !--- get mask, which is time-invariant ---
      call cpl_mct_aVect_getRAttr(bun_ocn%dom%lgrid,"mask",rmask,rcode)
      mask = nint(rmask)

      !--- setup chunk size for threading ---
#ifdef  UNTESTED_OPENMP
      nchunks   = omp_get_max_threads()
      chunksize = nloc/nchunks ! size of all but last chunk
      if (nchunks*chunksize /= nloc) chunklast = mod(nloc,chunksize)
      write(6,F01) 'FYI: routine is threaded, number of threads = ',nchunks
#else
      write(6,F01) 'FYI: this routine is not threaded'
#endif

   end if

   !----------------------------------------------------------------------------
   ! unpack input state fields
   !----------------------------------------------------------------------------
   call shr_timer_start(t1)
#if (1 == 0) 
   call cpl_mct_aVect_getRAttr(bun_ocn%data,"So_u   ",uocn ,rcode)
   call cpl_mct_aVect_getRAttr(bun_ocn%data,"So_v   ",vocn ,rcode)
   call cpl_mct_aVect_getRAttr(bun_ocn%data,"So_t   ",tocn ,rcode)

   call cpl_mct_aVect_getRAttr(bun_atm%data,"Sa_z   ",z    ,rcode)
   call cpl_mct_aVect_getRAttr(bun_atm%data,"Sa_u   ",u    ,rcode)
   call cpl_mct_aVect_getRAttr(bun_atm%data,"Sa_v   ",v    ,rcode)
   call cpl_mct_aVect_getRAttr(bun_atm%data,"Sa_ptem",ptem ,rcode)
   call cpl_mct_aVect_getRAttr(bun_atm%data,"Sa_shum",shum ,rcode)
   call cpl_mct_aVect_getRAttr(bun_atm%data,"Sa_dens",dens ,rcode)
   call cpl_mct_aVect_getRAttr(bun_atm%data,"Sa_tbot",tbot ,rcode)
#else
   k1  = cpl_mct_aVect_indexRA(bun_ocn%data,"So_u",perrWith=subName)
   k2  = cpl_mct_aVect_indexRA(bun_ocn%data,"So_v",perrWith=subName)
   k3  = cpl_mct_aVect_indexRA(bun_ocn%data,"So_t",perrWith=subName)
   k4  = cpl_mct_aVect_indexRA(bun_atm%data,"Sa_z",perrWith=subName)
   k5  = cpl_mct_aVect_indexRA(bun_atm%data,"Sa_u",perrWith=subName)
   k6  = cpl_mct_aVect_indexRA(bun_atm%data,"Sa_v",perrWith=subName)
   k7  = cpl_mct_aVect_indexRA(bun_atm%data,"Sa_ptem",perrWith=subName)
   k8  = cpl_mct_aVect_indexRA(bun_atm%data,"Sa_shum",perrWith=subName)
   k9  = cpl_mct_aVect_indexRA(bun_atm%data,"Sa_dens",perrWith=subName)
   k10 = cpl_mct_aVect_indexRA(bun_atm%data,"Sa_tbot",perrWith=subName)
   do n = 1,nloc
     uocn(n)  = bun_ocn%data%rAttr(k1 ,n)
     vocn(n)  = bun_ocn%data%rAttr(k2 ,n)
     tocn(n)  = bun_ocn%data%rAttr(k3 ,n)
     z(n)     = bun_atm%data%rAttr(k4 ,n)
     u(n)     = bun_atm%data%rAttr(k5 ,n)
     v(n)     = bun_atm%data%rAttr(k6 ,n)
     ptem(n)  = bun_atm%data%rAttr(k7 ,n)
     shum(n)  = bun_atm%data%rAttr(k8 ,n)
     dens(n)  = bun_atm%data%rAttr(k9 ,n)
     tbot(n)  = bun_atm%data%rAttr(k10,n)
   enddo
#endif
   !--- must fabricate "reasonable" data (using dead components?) ---
   if (fabricate) then
      mask =   1   ! ocn domain mask            ~ 0 <=> inactive cell
      tocn = 290.0 ! ocn temperature            ~ Kelvin
      uocn =   0.0 ! ocn velocity, zonal        ~ m/s
      vocn =   0.0 ! ocn velocity, meridional   ~ m/s
      z    =  55.0 ! atm height of bottom layer ~ m
      u    =   0.0 ! atm velocity, zonal        ~ m/s
      v    =   2.0 ! atm velocity, meridional   ~ m/s
      ptem = 301.0 ! atm potential temperature  ~ Kelvin
      shum = 1.e-2 ! atm specific humidity      ~ kg/kg
      dens =   1.0 ! atm density                ~ kg/m^3
      tbot = 300.0 ! atm temperature            ~ Kelvin
   endif

   call shr_timer_stop(t1)
   call shr_timer_start(t2)

   !----------------------------------------------------------------------------
   ! call the scientist-written physics routine
   !----------------------------------------------------------------------------
#ifdef  UNTESTED_OPENMP
   DO n=1,nchunks
      if (n /= nchunks) cs = chunksize
      if (n == nchunks) cs = chunklast
      m = (n-1)*chunksize
      call srfflx_ao( cs      ,   z(m),   u(m),   v(m),ptem(m), &
      &              shum  (m),dens(m),tbot(m),uocn(m),vocn(m), &
      &              tocn  (m),mask(m), sen(m), lat(m),lwup(m), &
      &              evap  (m),taux(m),tauy(m),tref(m),qref(m), &
      &              duu10n(m)                                  )
   END DO
#else
      call srfflx_ao( nloc    ,   z   ,   u   ,   v   ,ptem   , &
      &              shum     ,dens   ,tbot   ,uocn   ,vocn   , &
      &              tocn     ,mask   , sen   , lat   ,lwup   , &
      &              evap     ,taux   , tauy  ,tref   ,qref   , &
      &              duu10n                                     )
#endif

   call shr_timer_stop(t2)
   call shr_timer_start(t3)
 
   !----------------------------------------------------------------------------
   ! pack output flux fields
   !----------------------------------------------------------------------------
#if (1 == 0)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_sen "  ,sen   ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_lat "  ,lat   ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_lwup"  ,lwup  ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_evap"  ,evap  ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_taux"  ,taux  ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_tauy"  ,tauy  ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_tref"  ,tref  ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_qref"  ,qref  ,rcode)
   call cpl_mct_aVect_putRAttr(bun_flux%data,"Faoc_duu10n",duu10n,rcode)
#else
   k1  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_sen" ,perrWith=subName)
   k2  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_lat" ,perrWith=subName)
   k3  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_lwup",perrWith=subName)
   k4  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_evap",perrWith=subName)
   k5  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_taux",perrWith=subName)
   k6  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_tauy",perrWith=subName)
   k7  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_tref",perrWith=subName)
   k8  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_qref",perrWith=subName)
   k9  = cpl_mct_aVect_indexRA(bun_flux%data,"Faoc_duu10n",perrWith=subName)
   do n = 1,nloc
      bun_flux%data%rAttr(k1,n) = sen(n)
      bun_flux%data%rAttr(k2,n) = lat(n)
      bun_flux%data%rAttr(k3,n) = lwup(n)
      bun_flux%data%rAttr(k4,n) = evap(n)
      bun_flux%data%rAttr(k5,n) = taux(n)
      bun_flux%data%rAttr(k6,n) = tauy(n)
      bun_flux%data%rAttr(k7,n) = tref(n)
      bun_flux%data%rAttr(k8,n) = qref(n)
      bun_flux%data%rAttr(k9,n) = duu10n(n)
   enddo
#endif

   bun_flux%cnt = 1
   call shr_timer_stop(t3)
 
end subroutine flux_atmOcn

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: srfflx_ao -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!     
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - brought in from cpl5.
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as per Bill Large
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE srfflx_ao(imax  ,zbot  ,ubot  ,vbot  ,thbot ,   & 
           &         qbot  ,rbot  ,tbot  ,us    ,vs    ,   &
           &         ts    ,mask  ,sen   ,lat   ,lwup  ,   &
           &         evap  ,taux  ,tauy  ,tref  ,qref  ,   &
           &         duu10n                                )

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
   integer(IN),intent(in) ::       imax  ! array dimensions
   integer(IN),intent(in) :: mask (imax) ! ocn domain mask 0
   real(R8)   ,intent(in) :: zbot (imax) ! atm level height      (m)
   real(R8)   ,intent(in) :: ubot (imax) ! atm u wind            (m/s)
   real(R8)   ,intent(in) :: vbot (imax) ! atm v wind            (m/s)
   real(R8)   ,intent(in) :: thbot(imax) ! atm potential T       (K)
   real(R8)   ,intent(in) :: qbot (imax) ! atm specific humidity (kg/kg)
   real(R8)   ,intent(in) :: rbot (imax) ! atm air density       (kg/m^3)
   real(R8)   ,intent(in) :: tbot (imax) ! atm T                 (K) 
   real(R8)   ,intent(in) :: us   (imax) ! ocn u-velocity        (m/s)
   real(R8)   ,intent(in) :: vs   (imax) ! ocn v-velocity        (m/s)
   real(R8)   ,intent(in) :: ts   (imax) ! ocn temperature       (K)

   !--- output arguments -------------------------------
   real(R8),intent(out)  ::  sen  (imax) ! heat flux: sensible    (W/m^2)
   real(R8),intent(out)  ::  lat  (imax) ! heat flux: latent      (W/m^2)
   real(R8),intent(out)  ::  lwup (imax) ! heat flux: lw upward   (W/m^2)
   real(R8),intent(out)  ::  evap (imax) ! water flux: evap  ((kg/s)/m^2)
   real(R8),intent(out)  ::  taux (imax) ! surface stress, zonal      (N)
   real(R8),intent(out)  ::  tauy (imax) ! surface stress, maridional (N)
   real(R8),intent(out)  ::  tref (imax) ! diag:  2m ref height T     (K)
   real(R8),intent(out)  ::  qref (imax) ! diag:  2m ref humidity (kg/kg)
   real(R8),intent(out)  :: duu10n(imax) ! diag: 10m wind speed squared (m/s)^2
 
!EOP

   !--- local constants --------------------------------
   real(R8),parameter :: umin  =  0.5    ! minimum wind speed       (m/s)
   real(R8),parameter :: zref  = 10.0    ! reference height           (m)
   real(R8),parameter :: ztref =  2.0    ! reference height for air T (m)

   !--- local variables --------------------------------
   integer(IN) :: i      ! vector loop index
   real(R8)    :: vmag   ! surface wind magnitude   (m/s)
   real(R8)    :: thvbot ! virtual temperature      (K)
   real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
   real(R8)    :: delt   ! potential T difference   (K)
   real(R8)    :: delq   ! humidity difference      (kg/kg)
   real(R8)    :: stable ! stability factor
   real(R8)    :: rdn    ! sqrt of neutral exchange coeff (momentum) 
   real(R8)    :: rhn    ! sqrt of neutral exchange coeff (heat)     
   real(R8)    :: ren    ! sqrt of neutral exchange coeff (water)    
   real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)         
   real(R8)    :: rh     ! sqrt of exchange coefficient (heat)             
   real(R8)    :: re     ! sqrt of exchange coefficient (water)            
   real(R8)    :: ustar  ! ustar             
   real(R8)    :: qstar  ! qstar             
   real(R8)    :: tstar  ! tstar             
   real(R8)    :: hol    ! H (at zbot) over L
   real(R8)    :: xsq    ! ?
   real(R8)    :: xqq    ! ?
   real(R8)    :: psimh  ! stability function at zbot (momentum)
   real(R8)    :: psixh  ! stability function at zbot (heat and water)
   real(R8)    :: psix2  ! stability function at ztref reference height
   real(R8)    :: alz    ! ln(zbot/zref)
   real(R8)    :: al2    ! ln(zref/ztref)
   real(R8)    :: u10n   ! 10m neutral wind 
   real(R8)    :: tau    ! stress at zbot
   real(R8)    :: cp     ! specific heat of moist air
   real(R8)    :: bn     ! exchange coef funct for interpolation
   real(R8)    :: bh     ! exchange coef funct for interpolation
   real(R8)    :: fac    ! vertical interpolation factor

   !--- local functions --------------------------------
   real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
   real(R8)    :: cdn    ! function: neutral drag coeff at 10m
   real(R8)    :: psimhu ! function: unstable part of psimh
   real(R8)    :: psixhu ! function: unstable part of psimx
   real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
   real(R8)    :: Tk     ! dummy arg ~ temperature (K)
   real(R8)    :: xd     ! dummy arg ~ ?
 
   qsat(Tk)   = 640380.0 / exp(5107.4/Tk)
   cdn(Umps)  = 0.0027 / Umps + 0.000142 + 0.0000764 * Umps
   psimhu(xd) = log((1.0+xd*(2.0+xd))*(1.0+xd*xd)/8.0) - 2.0*atan(xd) + 1.571
   psixhu(xd) = 2.0 * log((1.0 + xd*xd)/2.0)
 
!-------------------------------------------------------------------------------
! PURPOSE:
!   computes atm/ocn surface fluxes
!
! NOTES: 
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + sen + lat
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
! 
! ASSUMPTIONS:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!-------------------------------------------------------------------------------
 
   al2 = log(zref/ztref)

   DO i=1,imax
     if (mask(i) /= 0) then
    
        !--- compute some needed quantities ---
        vmag   = max(umin, sqrt( (ubot(i)-us(i))**2 + (vbot(i)-vs(i))**2) )
        thvbot = thbot(i) * (1.0 + cpl_const_zvir * qbot(i)) ! virtual temp (K)
        ssq    = 0.98 * qsat(ts(i)) / rbot(i)      ! sea surf hum (kg/kg)
        delt   = thbot(i) - ts(i)                  ! pot temp diff (K)
        delq   = qbot(i) - ssq                     ! spec hum dif (kg/kg)
        alz    = log(zbot(i)/zref) 
        cp     = cpl_const_cpdair*(1.0 + cpl_const_cpvir*ssq) 
   
        !------------------------------------------------------------
        ! first estimate of Z/L and ustar, tstar and qstar
        !------------------------------------------------------------
   
        !--- neutral coefficients, z/L = 0.0 ---
        stable = 0.5 + sign(0.5 , delt)
        rdn    = sqrt(cdn(vmag))
        rhn    = (1.0-stable) * 0.0327 + stable * 0.018 
        ren    = 0.0346 
   
        !--- ustar, tstar, qstar ---
        ustar = rdn * vmag
        tstar = rhn * delt  
        qstar = ren * delq  
   
        !--- compute stability & evaluate all stability functions ---
        hol  = cpl_const_karman*cpl_const_g*zbot(i)*  &
               (tstar/thvbot+qstar/(1.0/cpl_const_zvir+qbot(i)))/ustar**2
        hol  = sign( min(abs(hol),10.0), hol )
        stable = 0.5 + sign(0.5 , hol)
        xsq    = max(sqrt(abs(1.0 - 16.0*hol)) , 1.0)
        xqq    = sqrt(xsq)
        psimh  = -5.0*hol*stable + (1.0-stable)*psimhu(xqq)
        psixh  = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
   
        !--- shift wind speed using old coefficient ---
        rd   = rdn / (1.0 + rdn/cpl_const_karman*(alz-psimh))
        u10n = vmag * rd / rdn 
   
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.0-stable)*0.0327 + stable * 0.018 
    
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0 + rdn/cpl_const_karman*(alz-psimh)) 
        rh = rhn / (1.0 + rhn/cpl_const_karman*(alz-psixh)) 
        re = ren / (1.0 + ren/cpl_const_karman*(alz-psixh)) 
   
        !--- update ustar, tstar, qstar using updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
    
        !------------------------------------------------------------
        ! iterate to converge on Z/L, ustar, tstar and qstar
        !------------------------------------------------------------
    
        !--- compute stability & evaluate all stability functions ---
        hol  = cpl_const_karman*cpl_const_g*zbot(i)* &
               (tstar/thvbot+qstar/(1.0/cpl_const_zvir+qbot(i)))/ustar**2
        hol  = sign( min(abs(hol),10.0), hol )
        stable = 0.5 + sign(0.5 , hol)
        xsq    = max(sqrt(abs(1.0 - 16.0*hol)) , 1.0)
        xqq    = sqrt(xsq)
        psimh  = -5.0*hol*stable + (1.0-stable)*psimhu(xqq)
        psixh  = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
    
        !--- shift wind speed using old coeffs ---
        rd   = rdn / (1.0 + rdn/cpl_const_karman*(alz-psimh))
        u10n = vmag * rd/rdn 
    
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.0 - stable)*0.0327 + stable * 0.018 
   
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0 + rdn/cpl_const_karman*(alz-psimh)) 
        rh = rhn / (1.0 + rhn/cpl_const_karman*(alz-psixh)) 
        re = ren / (1.0 + ren/cpl_const_karman*(alz-psixh)) 
    
        !--- update ustar, tstar, qstar using updated, shifted coeffs ---
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
    
        !------------------------------------------------------------
        ! compute the fluxes
        !------------------------------------------------------------
    
        tau = rbot(i) * ustar * ustar 
       
        !--- momentum flux ---
        taux(i) = tau * (ubot(i)-us(i)) / vmag 
        tauy(i) = tau * (vbot(i)-vs(i)) / vmag 
        
        !--- heat flux ---
        sen (i) =                cp * tau * tstar / ustar 
        lat (i) =  cpl_const_latvap * tau * qstar / ustar
        lwup(i) = -cpl_const_stebol * ts(i)**4 
      
        !--- water flux ---
        evap(i) = lat(i)/cpl_const_latvap 
    
        !------------------------------------------------------------
        ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
        !------------------------------------------------------------
        hol = hol*ztref/zbot(i)
        xsq = max( 1.0, sqrt(abs(1.0-16.0*hol)) )
        xqq = sqrt(xsq)
        psix2   = -5.0*hol*stable + (1.0-stable)*psixhu(xqq)
        fac     = (rh/cpl_const_karman) * (alz + al2 - psixh + psix2 )
        tref(i) = thbot(i) - delt*fac 
        tref(i) = tref(i) - 0.01*ztref   ! pot temp to temp correction
        fac     = (re/cpl_const_karman) * (alz + al2 - psixh + psix2 )
        qref(i) =  qbot(i) - delq*fac
    
        duu10n(i) = u10n*u10n ! 10m wind speed squared

     else
        !------------------------------------------------------------
        ! no valid data here -- out of domain
        !------------------------------------------------------------
        sen   (i) = cpl_const_spval  ! sensible         heat flux  (W/m^2)
        lat   (i) = cpl_const_spval  ! latent           heat flux  (W/m^2)
        lwup  (i) = cpl_const_spval  ! long-wave upward heat flux  (W/m^2)
        evap  (i) = cpl_const_spval  ! evaporative water flux ((kg/s)/m^2)
        taux  (i) = cpl_const_spval  ! x surface stress (N)
        tauy  (i) = cpl_const_spval  ! y surface stress (N)
        tref  (i) = cpl_const_spval  !  2m reference height temperature (K)
        qref  (i) = cpl_const_spval  !  2m reference height humidity (kg/kg)
        duu10n(i) = cpl_const_spval  ! 10m wind speed squared (m/s)^2
     endif
   ENDDO 

END subroutine srfflx_ao

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: flux_albo - ocean albedo calculation
!
! !DESCRIPTION:
!   if flux\_albav/=0 (ie. "on" or "true") \\
!      Compute four effective daily avg surface albedos for all
!      combinations of visible/near-infrared and direct/diffuse radiation
!      without accounting for zenith angle (ie. a "daily average" albedo) \\
!   else \\
!      Compute four surface albedos for all combinations of visible/
!      near-infrared and direct/diffuse radiation, accounting for
!      instantaneous zenith angle
!
! !REMARKS:
!   o upon input, albedos are assumed to be a 60 degree reference albedo
!   o albedos are computed by taking the 60 deg reference albedo
!     and then adjusting this value based on zenith angle
!   o Albedos are independent of spectral interval and other physical
!     factors such as surface wind speed.
!
!   For more details see Briegleb, Bruce P., 1992: Delta-Eddington
!   Approximation for Solar Radiation in the NCAR Community Climate
!   Model, Journal of Geophysical Research, Vol 97, D7, pp7603-7612.
!
! !REVISION HISTORY:
!    198x        - CCM1, original version
!    1992-Jun    - J. Rosinski -- standardized
!    1994-May    - J. Rosinski -- rewritten for land only
!    1994-Jul    - B. Kauffman -- rewritten for ocean only
!    2002-Oct-26 - R. Jacob -- Rewritten for cpl6
!
! !INTERFACE: ------------------------------------------------------------------

subroutine flux_albo(date,bun_ocn)

! !USES:

   use shr_orb_mod ! orbital constants and methods

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout ) :: bun_ocn ! ocn albedo fields
   type(shr_date)  ,intent(in )    :: date    ! current date

!EOP
   !--- local ---
   real(R8),save,allocatable :: anidr (:)   ! albedo: near infrared, direct
   real(R8),save,allocatable :: avsdr (:)   ! albedo: visible      , direct
   real(R8),save,allocatable :: anidf (:)   ! albedo: near infrared, diffuse
   real(R8),save,allocatable :: avsdf (:)   ! albedo: visible      , diffuse
   real(R8),save,allocatable :: lats  (:)   ! latitudes
   real(R8),save,allocatable :: lons  (:)   ! longitudes

   integer(IN),save          :: nloc = -1   ! size of local data array

   integer(IN) :: n                   ! loop index
   real(R8)    :: rlat                ! gridcell latitude in radians
   real(R8)    :: rlon                ! gridcell longitude in radians
   real(R8)    :: eccen               ! orbital eccentricity
   real(R8)    :: mvelpp              ! moving vernal equinox long
   real(R8)    :: lambm0              ! Mean long of perihelion (rad)
   real(R8)    :: obliqr              ! Earths obliquity (rad)
   real(R8)    :: delta               ! Solar declination angle in radians
   real(R8)    :: eccf                ! Earth-sun distance factor
   real(R8)    :: cosz                ! Cosine of solar zenith angle
   real(R8)    :: calday              ! calendar day including fraction, at 0e
   logical     :: first_call = .true. ! flags 1st envocation of this routine
   integer(IN) :: n_o                 ! number of points
   integer(IN) :: rcode               ! error code

   real(R8),parameter :: albdif = 0.06 ! 60 deg reference albedo, diffuse
   real(R8),parameter :: albdir = 0.07 ! 60 deg reference albedo, direct 

   SAVE

   !--- formats ---
   character(*),parameter :: subName = '(flux_albo) '
   character(*),parameter :: F00 = '(  "(flux_albo) ",4a  )'
   character(*),parameter :: F01 = '(  "(flux_albo) ",a,i5)'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !--- allocate albedo, lat & lon arrays ---
   if (nloc < 1 ) then
      nloc = cpl_mct_aVect_lsize(bun_ocn%data)
      allocate( anidr(nloc),avsdr(nloc),anidf(nloc),avsdf(nloc))
      allocate( lats(nloc),lons(nloc))
      nloc =1
   endif

   n_o = cpl_mct_aVect_lsize(bun_ocn%data)

   IF (cpl_control_fluxAlbav) THEN
      !--------------------------------------------------------------------
      if (first_call) write(6,F00) 'ocn albedo is zenith angle independant'
      !--------------------------------------------------------------------
      do n=1,n_o
         anidr(n) = albdir
         avsdr(n) = albdir
         anidf(n) = albdif
         avsdf(n) = albdif
      end do
   ELSE
      !--------------------------------------------------------------------
      if (first_call) write(6,F00) 'ocn albedo is zenith angle dependant'
      !--------------------------------------------------------------------

      !--- julian day ---
      calday = shr_date_getJulian(date,cpl_control_fluxAshift)

      !--- solar declination ---
      eccen  = cpl_control_orbEccen
      mvelpp = cpl_control_orbMvelpp
      lambm0 = cpl_control_orbLambm0
      obliqr = cpl_control_orbObliqr
      call shr_orb_decl(calday,eccen,mvelpp,lambm0,obliqr,delta,eccf)
      call cpl_mct_aVect_getRAttr(bun_ocn%dom%lGrid,"lon",lons    ,rcode)
      call cpl_mct_aVect_getRAttr(bun_ocn%dom%lGrid,"lat",lats    ,rcode)

      do n=1,n_o
         rlat = cpl_const_deg2rad * lats(n)
         rlon = cpl_const_deg2rad * lons(n)
         cosz = shr_orb_cosz( calday, rlat, rlon, delta )
         if (cosz  >  0.0) then !--- sun hit --
            anidr(n) = (.026/(cosz**1.7 + 0.065)) +   &
                       (.150*(cosz      - 0.10 )  *   &
                             (cosz      - 0.50 )  *   &
                             (cosz      - 1.00 )  )
            avsdr(n) = anidr(n)
            anidf(n) = albdif
            avsdf(n) = albdif
         else !--- dark side of earth ---
            anidr(n) = 1.0
            avsdr(n) = 1.0
            anidf(n) = 1.0
            avsdf(n) = 1.0
         end if
      end do
   END IF 

   call cpl_mct_aVect_putRAttr(bun_ocn%data,"So_anidr",anidr,rcode)
   call cpl_mct_aVect_putRAttr(bun_ocn%data,"So_avsdr",avsdr,rcode)
   call cpl_mct_aVect_putRAttr(bun_ocn%data,"So_anidf",anidf,rcode)
   call cpl_mct_aVect_putRAttr(bun_ocn%data,"So_avsdf",avsdf,rcode)

   bun_ocn%cnt = 1
 
   first_call=.false.

END subroutine flux_albo

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: flux_albi - ice albedo modification
!
! !DESCRIPTION:
!   if flux\_albav/=0 (ie. "on" or "true") \\
!     Impose a zenith angle dependance on the ice model "reference albedo".
!     Currently this only involves setting albedos to zero
!     on the dark side of the earth. \\
!   else \\
!     do not alter ice albedos
!
! !REMARKS:
!   o upon input, albedos are assumed to be a 60 degree reference albedo
!
! !REVISION HISTORY:
!   199x-       - B. Kauffman -- original cpl5 version
!   2002-Oct-26 - R. Jacob -- rewritten for cpl6
!
! !INTERFACE: ------------------------------------------------------------------

subroutine flux_albi(date,bun_ice)

! !USES:

   use shr_orb_mod

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout)  :: bun_ice ! contains ice albedo fields
   type(shr_date)  ,intent(in )    :: date    ! current date

!EOP

   !--- local ---
   real(R8),save,allocatable :: lats  (:)  ! latitudes
   real(R8),save,allocatable :: lons  (:)  ! longitudes
   integer(IN),save          :: nloc = -1  ! size of local data array

   integer(IN) :: n                   ! generic loop index
   integer(IN) :: anidf               ! aVect index for anidf
   integer(IN) :: avsdf               ! aVect index for avsdf
   integer(IN) :: anidr               ! aVect index for anidr
   integer(IN) :: avsdr               ! aVect index for avsdr
   real(R8)    :: rlat                ! gridcell latitude in radians
   real(R8)    :: rlon                ! gridcell longitude in radians
   real(R8)    :: eccen               ! orbital eccentricity
   real(R8)    :: mvelpp              ! moving vernal equinox long
   real(R8)    :: lambm0              ! Mean long of perihelion (rad)
   real(R8)    :: obliqr              ! Earths obliquity (rad)
   real(R8)    :: delta               ! Solar declination angle in radians
   real(R8)    :: eccf                ! Earth-sun distance factor
   real(R8)    :: cosz                ! Cosine of solar zenith angle
   real(R8)    :: calday              ! calendar day including fraction, at 0e
   integer(IN) :: n_o                 ! number of points
   integer(IN) :: rcode               ! error code
   logical     :: first_call = .true. ! flags 1st envocation of this routine

   real(R8),parameter :: albdif = 0.06 ! 60 deg ref albedo, diffuse
   real(R8),parameter :: albdir = 0.07 ! 60 deg ref albedo, direct 

   SAVE

   !--- formats ---
   character(*),parameter :: subName = '(flux_albi) '
   character(*),parameter :: F00 = '(  "(flux_albi) ",4a  )'
   character(*),parameter :: F01 = '(  "(flux_albi) ",a,i5)'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! allocate lat & lon arrays 
   !--------------------------------------------------------------------
   if (nloc < 1 ) then
      nloc = cpl_mct_aVect_lsize(bun_ice%data)
      allocate( lats(nloc),lons(nloc))
      nloc = 1
   endif

   IF (cpl_control_fluxAlbav) THEN
      !--------------------------------------------------------------------
      if (first_call) write(6,F00) 'ice albedo is zenith angle independant'
      !--------------------------------------------------------------------
   ELSE
      !--------------------------------------------------------------------
      if (first_call) write(6,F00) 'ice albedo is zenith angle dependant'
      !--------------------------------------------------------------------

      n_o = cpl_mct_aVect_lsize(bun_ice%data)

      anidr = cpl_mct_aVect_indexRA(bun_ice%data,"Si_anidr") 
      avsdr = cpl_mct_aVect_indexRA(bun_ice%data,"Si_avsdr") 
      anidf = cpl_mct_aVect_indexRA(bun_ice%data,"Si_anidf") 
      avsdf = cpl_mct_aVect_indexRA(bun_ice%data,"Si_avsdf") 

      !--- julian day ---
      calday = shr_date_getJulian(date,cpl_control_fluxAshift)

      !--- solar declination ---
      eccen  = cpl_control_orbEccen
      mvelpp = cpl_control_orbMvelpp
      lambm0 = cpl_control_orbLambm0
      obliqr = cpl_control_orbObliqr
      call shr_orb_decl(calday,eccen,mvelpp,lambm0,obliqr,delta,eccf)
      call cpl_mct_aVect_getRAttr(bun_ice%dom%lGrid,"lon",lons    ,rcode)
      call cpl_mct_aVect_getRAttr(bun_ice%dom%lGrid,"lat",lats    ,rcode)

      do n=1,n_o
         rlat = cpl_const_deg2rad * lats(n)
         rlon = cpl_const_deg2rad * lons(n)
         cosz = shr_orb_cosz( calday, rlat, rlon, delta )
         if (cosz < 0.0) then !--- dark side of earth ---
            bun_ice%data%rAttr(anidr,n) = 1.0
            bun_ice%data%rAttr(avsdr,n) = 1.0
            bun_ice%data%rAttr(anidf,n) = 1.0
            bun_ice%data%rAttr(avsdf,n) = 1.0
         end if
      end do

   END IF 
 
   first_call=.false.

END subroutine flux_albi

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: flux_solar - compute atm/ocn absorbed short-wave (net sw)
!
! !DESCRIPTION:
!   compute atm/ocn absorbed short-wave (net sw)
!
! !REVISION HISTORY:
!   2000-Jan-03 - B. Kauffman -- original cpl5 version
!   2002-Oct-26 - R. Jacob -- rewritten for cpl6
!
! !INTERFACE: ------------------------------------------------------------------

subroutine flux_solar(bun_atm,bun_ocn,bun_aoflux)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in  ) :: bun_atm    ! contains atm sw down fields
   type(cpl_bundle),intent(in  ) :: bun_ocn    ! contains ocn albedo fields
   type(cpl_bundle),intent(out ) :: bun_aoflux ! contains a/o net-sw fields

!EOP

   !--- local ---
   integer(IN) :: swnet       ! aV index ~ short wave: net
   integer(IN) :: anidr       ! aV index ~ albedo: near-infra, direct
   integer(IN) :: avsdr       ! aV index ~ albedo: visible   , direct
   integer(IN) :: anidf       ! aV index ~ albedo: near-infra, diffuse
   integer(IN) :: avsdf       ! aV index ~ albedo: visible   , diffuse
   integer(IN) :: swndr       ! aV index ~ short wave: near-infra, direct
   integer(IN) :: swvdr       ! aV index ~ short wave: visible   , direct
   integer(IN) :: swndf       ! aV index ~ short wave: near-infra, diffuse
   integer(IN) :: swvdf       ! aV index ~ short wave: visible   , diffuse

   integer(IN) :: n           ! loop index
   integer(IN) :: n_o         ! number of ocean points

   SAVE

   !--- formats ---
   character(*),parameter :: subName = '(flux_solar) '
   character(*),parameter :: F00 = '(  "(flux_solar) ",4a  )'
   character(*),parameter :: F01 = '(  "(flux_solar) ",a,i5)'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   swnet = cpl_mct_aVect_indexRA(bun_aoflux%data,"Faoc_swnet") 

   anidr = cpl_mct_aVect_indexRA(bun_ocn%data,"So_anidr") 
   avsdr = cpl_mct_aVect_indexRA(bun_ocn%data,"So_avsdr") 
   anidf = cpl_mct_aVect_indexRA(bun_ocn%data,"So_anidf") 
   avsdf = cpl_mct_aVect_indexRA(bun_ocn%data,"So_avsdf") 
   
   swndr = cpl_mct_aVect_indexRA(bun_atm%data,"Faxa_swndr") 
   swvdr = cpl_mct_aVect_indexRA(bun_atm%data,"Faxa_swvdr") 
   swndf = cpl_mct_aVect_indexRA(bun_atm%data,"Faxa_swndf") 
   swvdf = cpl_mct_aVect_indexRA(bun_atm%data,"Faxa_swvdf") 

   n_o = cpl_mct_aVect_lsize(bun_aoflux%data)

   do n=1,n_o
      bun_aoflux%data%rAttr(swnet,n) =                                      &
      & ((1.0 - bun_ocn%data%rAttr(anidr,n))*bun_atm%data%rAttr(swndr,n)) + &
      & ((1.0 - bun_ocn%data%rAttr(avsdr,n))*bun_atm%data%rAttr(swvdr,n)) + &
      & ((1.0 - bun_ocn%data%rAttr(anidf,n))*bun_atm%data%rAttr(swndf,n)) + &
      & ((1.0 - bun_ocn%data%rAttr(avsdf,n))*bun_atm%data%rAttr(swvdf,n))
   enddo

end subroutine flux_solar

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: flux_epbal -  Calculate precip/runoff adjustment factor
!
! !DESCRIPTION:
!   adjust precip (an atm output flux) and runoff (a lnd output) sent to
!   ice \& ocn by a scalar factor f, st  P'+R'+E = f(P+R)+E=0.
!   Why?  This will insure a net zero fresh-water flux into ocn+ice.
!   This could be used to compensate for fresh-water flux imbalances,
!   eg. due to the lack of river runoff from the lnd model.
!
! !REVISION HISTORY:
!   199x        - B. Kauffman -- Original cpl5 version
!   2003-Feb-17 - R. Jacob -- rewritten for cpl6
!
! !INTERFACE: ------------------------------------------------------------------

subroutine flux_epbal(date,bun_aoflux,bun_i2c,bun_prec,bun_r2c,bun_frac)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date)  ,intent(in   ) :: date       ! current date
   type(cpl_bundle),intent(in   ) :: bun_i2c    ! ice to cpl bundle: ice evap
   type(cpl_bundle),intent(in   ) :: bun_aoflux ! a/o flux bundle  : ocn evap
   type(cpl_bundle),intent(inout) :: bun_prec   ! a/x precip bundle: i+o prec
   type(cpl_bundle),intent(inout) :: bun_r2c    ! runoff bundle    : ocn roff
   type(cpl_bundle),intent(in   ) :: bun_frac   ! fractions on ocn domain

!EOP

   !--- local ---
   integer(IN),parameter :: k_area = 1 ! index: area of ocn domain
   integer(IN),parameter :: k_prec = 2 ! index: water ~ precipitation
   integer(IN),parameter :: k_evap = 3 ! index: water ~ evaporation
   integer(IN),parameter :: k_roff = 4 ! index: water ~ runoff
   real(R8)              :: psum(4)    ! partial/local sum of area,prec,evap,roff
   real(R8)              :: gsum(4)    !   full/global sum of area,prec,evap,roff
   real(R8)              :: tprec      ! total precip
   real(R8)              :: tevap      ! total evap
   real(R8)              :: troff      ! total runoff
   real(R8)              :: tarea      ! total area
   real(R8)              :: da         ! area of one ocn grid cell = dth*dph
   real(R8)              :: dai        ! area of ocn grid covered by ice
   real(R8)              :: dao        ! area of ocn grid covered by atm
   real(R8)              :: factor     ! prec adjustment factor: evap/prec
   integer(IN),save      :: nloc       ! size of local data array
   integer(IN),save      :: k_xrain    ! i+o-rain aVect indicies
   integer(IN),save      :: k_xsnow    ! i+o-snow aVect indicies
   integer(IN),save      :: k_ievap    ! ice-evap aVect indicies
   integer(IN),save      :: k_oevap    ! ocn-evap aVect indicies
   integer(IN),save      :: k_oroff    ! ocn-roff aVect indicies
   integer(IN),save      :: k_ifrac    ! ice-frac aVect indicies
   integer(IN),save      :: k_afrac    ! atm-frac aVect indicies
   real(R8),allocatable,save ::  area(:)  ! cell area
   real(R8),allocatable,save :: imask(:)  ! domain mask
   real(R8),allocatable      ::  mask(:)  ! domain mask
   logical    ,save :: first_call =.true. ! flags 1st invocation of routine
   integer(IN)      :: n                  ! generic loop index
   integer(IN)      :: rcode              ! return code
   integer(IN)      :: year,month,day,sec ! date & time info

   !----- formats -----
   character(*),parameter :: subName = '(flux_epbal) '
   character(*),parameter :: F00 = "('(flux_epbal) ',4a)"
   character(*),parameter :: F01 = "('(flux_epbal) ',a,3e11.3,a,f9.6)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (cpl_control_fluxEPbal(1:3) == 'off') return

   if (first_call) then
      write(6,F00) "apply adjustment factor to prec & runoff over ice & ocn"
      nloc  =  cpl_mct_aVect_lsize(bun_aoflux%data)

      k_xrain = cpl_mct_aVect_indexRA(  bun_prec%data,"Faxc_rain")
      k_xsnow = cpl_mct_aVect_indexRA(  bun_prec%data,"Faxc_snow")
      k_ievap = cpl_mct_aVect_indexRA(   bun_i2c%data,"Faii_evap")
      k_oevap = cpl_mct_aVect_indexRA(bun_aoflux%data,"Faoc_evap")
      k_oroff = cpl_mct_aVect_indexRA(   bun_r2c%data,"Forr_roff")
      k_ifrac = cpl_mct_aVect_indexRA(  bun_frac%data,"ifrac")
      k_afrac = cpl_mct_aVect_indexRA(  bun_frac%data,"afrac")

      allocate( area(nloc),mask(nloc),imask(nloc))
      call cpl_mct_aVect_getRAttr(bun_aoflux%dom%lGrid,"aream",area,rcode)
      call cpl_mct_aVect_getRAttr(bun_aoflux%dom%lGrid,"mask" ,mask,rcode)

      imask(:) = 0  
      where (mask /= 0.0) imask = 1

      deallocate(mask)
   end if


   !----------------------------------------------------------------------------
   ! compute local integrals
   !----------------------------------------------------------------------------
   psum(:) = 0 ! local/parial sum
   do n=1,nloc
      if (imask(n) /= 0) then
         da  = area(n)
         dai = da*bun_frac%data%rAttr(k_ifrac,n) 
         dao = da*bun_frac%data%rAttr(k_afrac,n)
         psum(k_area) = psum(k_area) + dai + dao
         psum(k_prec) = psum(k_prec) + dai*( bun_prec%data%rAttr(k_xrain,n) + &
         &                                   bun_prec%data%rAttr(k_xsnow,n)   )
         psum(k_prec) = psum(k_prec) + dao*( bun_prec%data%rAttr(k_xrain,n) + &
         &                                   bun_prec%data%rAttr(k_xsnow,n)   )
         psum(k_evap) = psum(k_evap) + dai*   bun_i2c%data%rAttr(k_ievap,n)
         psum(k_evap) = psum(k_evap) + dao*bun_aoflux%data%rAttr(k_oevap,n)
         psum(k_roff) = psum(k_roff) + da *   bun_r2c%data%rAttr(k_oroff,n)
      endif
   enddo 

   !----------------------------------------------------------------------------
   ! compute factor (on master process only)
   !----------------------------------------------------------------------------
   gsum(:) = 0 ! global sum
   call shr_mpi_sum(psum,gsum,cpl_comm_comp,subName//" psum")
       
   if (cpl_comm_comp_pid == 0) then

      tarea = gsum(k_area)
      tprec = gsum(k_prec)/tarea
      tevap = gsum(k_evap)/tarea
      troff = gsum(k_roff)/tarea

      if (cpl_control_fluxEPbal(1:3) == 'ocn') then
         !-------------------------------------------------------------
         ! use factor supplied by ocn, NOTE: may not cause E+f(P+R)=0 
         !-------------------------------------------------------------
         if (first_call) write(6,F01) 'use adjustment factor provided by ocn'

         factor = cpl_control_fluxEPfac
         if (factor .le. 0.0) then
            write(6,F01) 'WARNING: factor from ocn = ',factor
            write(6,F01) 'WARNING: resetting factor to 1.0'
            factor = 1.0
         end if

      else if (cpl_control_fluxEPbal(1:4) == 'inst') then
         !-------------------------------------------------------------
         ! compute factor st f(P+R)+E=0 at every timestep 
         !-------------------------------------------------------------
         if (first_call) write(6,F00) "choosing factor st P'+E+R'=0"
 
         if ( (tprec+troff) > 0.0) then
            factor = -tevap/(tprec+troff)
         else
            factor=1.0
            write(6,F01) 'WARNING: avg  P,R,(P+R) = ',tprec,troff,tprec+troff
            write(6,F01) 'WARNING: setting factor = 1.0'
         end if
      else
         !-------------------------------------------------------------
         ! invalid option
         !-------------------------------------------------------------
         write(6,F00) 'ERROR: unknown epbal option: ',cpl_control_fluxEPbal
         call shr_sys_abort(subName)
      end if

     !--- document factor ---
     call shr_date_getYMD(date,year,month,day,sec)

     if ( factor < .75  .or.  1.25 < factor ) then
        write(6,F01) '<E> <P> <R> =',tevap,tprec,troff,'  f =',factor
        write(6,F01) 'WARNING: erroneous adjustment factor?'
     else if ( day+sec .eq. 1 ) then
        write(6,F01) '<E> <P> <R> =',tevap,tprec,troff,'  f =',factor
     else if ( dbug > 1 .and. sec == 0 ) then
        write(6,F01) '<E> <P> <R> =',tevap,tprec,troff,'  f =',factor
     end if

   endif

   !----------------------------------------------------------------------------
   ! apply factor locally on all processes
   !----------------------------------------------------------------------------
   call shr_mpi_bcast(factor,cpl_comm_comp,subName//" factor")

   do n=1,nloc
      if (imask(n) /= 0) then
         bun_prec%data%rAttr(k_xrain,n) = bun_prec%data%rAttr(k_xrain,n)*factor
         bun_prec%data%rAttr(k_xsnow,n) = bun_prec%data%rAttr(k_xsnow,n)*factor
          bun_r2c%data%rAttr(k_oroff,n) =  bun_r2c%data%rAttr(k_oroff,n)*factor
      end if
   end do

   first_call = .false.

end subroutine flux_epbal

!===============================================================================

end module flux_mod
