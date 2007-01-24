!===============================================================================
! CVS: $Id: cpl_fields_mod.F90,v 1.1.1.1 2005/02/03 22:29:02 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/dead6build/cpl_fields_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_fields_mod -- coupler/component  list of exchanged fields
!
! !DESCRIPTION:
!     Contains \& coordinates fields names and buffer locations for use
!     in coupler/component interface.
!
! !REVISION HISTORY:
!     2002 Feb 11 - full, realistic list of fields
!     2001 Apr 13 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_fields_mod

! !USES:

   use cpl_mct_mod     ! mct
   use cpl_kind_mod    ! kinds

   private  ! except

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_fields_getField    ! returns string for nth aVect attribute 
   public :: cpl_fields_getLongName ! returns netCDF longname and unit strings

   !----------------------------------------------------------------------------
   ! component names
   !----------------------------------------------------------------------------

   character(32),parameter,public :: cpl_fields_atmname='atm'
   character(32),parameter,public :: cpl_fields_ocnname='ocn'
   character(32),parameter,public :: cpl_fields_icename='ice'
   character(32),parameter,public :: cpl_fields_lndname='lnd'
   character(32),parameter,public :: cpl_fields_rtmname='roff'
   character(32),parameter,public :: cpl_fields_cplname='cpl'

   !----------------------------------------------------------------------------
   ! "info-buffer" index of entries 
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: cpl_fields_ibuf_total    = 100 ! size of info-buffer

   integer(IN),parameter,public :: cpl_fields_ibuf_rcode    =   1 ! error code
   integer(IN),parameter,public :: cpl_fields_ibuf_cdate    =   2 ! current date: yymmdd
   integer(IN),parameter,public :: cpl_fields_ibuf_sec      =   3 ! elapsed sec on date
   integer(IN),parameter,public :: cpl_fields_ibuf_ncpl     =   4 ! cpl comm's per day
   integer(IN),parameter,public :: cpl_fields_ibuf_nfields  =  10
   integer(IN),parameter,public :: cpl_fields_ibuf_gsize    =  11
   integer(IN),parameter,public :: cpl_fields_ibuf_lsize    =  12
   integer(IN),parameter,public :: cpl_fields_ibuf_gisize   =  13
   integer(IN),parameter,public :: cpl_fields_ibuf_gjsize   =  14
   integer(IN),parameter,public :: cpl_fields_ibuf_lisize   =  15
   integer(IN),parameter,public :: cpl_fields_ibuf_ljsize   =  16
   integer(IN),parameter,public :: cpl_fields_ibuf_stopeod  =  19
   integer(IN),parameter,public :: cpl_fields_ibuf_stopnow  =  20
   integer(IN),parameter,public :: cpl_fields_ibuf_resteod  =  21
   integer(IN),parameter,public :: cpl_fields_ibuf_restnow  =  22
   integer(IN),parameter,public :: cpl_fields_ibuf_histeod  =  23
   integer(IN),parameter,public :: cpl_fields_ibuf_histnow  =  24
   integer(IN),parameter,public :: cpl_fields_ibuf_histavg  =  25
   integer(IN),parameter,public :: cpl_fields_ibuf_diageod  =  26
   integer(IN),parameter,public :: cpl_fields_ibuf_diagnow  =  27
   integer(IN),parameter,public :: cpl_fields_ibuf_infotim  =  28
   integer(IN),parameter,public :: cpl_fields_ibuf_infobug  =  29
   integer(IN),parameter,public :: cpl_fields_ibuf_precadj  =  31 ! precip adjustment factor (* 1.0e+6)
   integer(IN),parameter,public :: cpl_fields_ibuf_ashift   =  32 ! albedo calculation time shift
   integer(IN),parameter,public :: cpl_fields_ibuf_nbasins  =  33 ! number of active runoff basins
   integer(IN),parameter,public :: cpl_fields_ibuf_xalbic   =  34 ! request extra albedo solar init msg
   integer(IN),parameter,public :: cpl_fields_ibuf_inimask  =  36 ! flag cpl to send back domain spec for lnd
   integer(IN),parameter,public :: cpl_fields_ibuf_dead     =  37 ! non-0 <=> dead model
   integer(IN),parameter,public :: cpl_fields_ibuf_domain   =  40
   integer(IN),parameter,public :: cpl_fields_ibuf_userest  =  41 ! non-0 <=> use restart data sent to cpl

   integer(IN),parameter,public :: cpl_fields_rbuf_total    =  50 ! size of real info-buffer

   integer(IN),parameter,public :: cpl_fields_rbuf_spval    =   1 ! the special value
   integer(IN),parameter,public :: cpl_fields_rbuf_eccen    =  10 ! Earth's eccentricity
   integer(IN),parameter,public :: cpl_fields_rbuf_obliqr   =  11 ! Earth's Obliquity
   integer(IN),parameter,public :: cpl_fields_rbuf_lambm0   =  12 ! longitude of perihelion at v-equinox
   integer(IN),parameter,public :: cpl_fields_rbuf_mvelpp   =  13 ! Earth's Moving vernal equinox of orbit +pi

   !----------------------------------------------------------------------------
   ! initial fields, generally a domain description 
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: cpl_fields_grid_total  = 7

   character(*),parameter,public :: cpl_fields_grid_fields = &
      &'lat&
      &:lon&
      &:area&
      &:aream&
      &:index&
      &:mask&
      &:pid'

   integer(IN),parameter,public :: cpl_fields_grid_lat    = 1    ! lat from component
   integer(IN),parameter,public :: cpl_fields_grid_lon    = 2    ! lon from component
   integer(IN),parameter,public :: cpl_fields_grid_area   = 3    ! area from component
   integer(IN),parameter,public :: cpl_fields_grid_aream  = 4    ! area from mapping file
   integer(IN),parameter,public :: cpl_fields_grid_index  = 5    ! global index
   integer(IN),parameter,public :: cpl_fields_grid_mask   = 6    ! mask, 0 = inactive cell
   integer(IN),parameter,public :: cpl_fields_grid_pid    = 7    ! proc id number

   !----------------------------------------------------------------------------
   ! atm fields 
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: cpl_fields_a2c_total  = 19

   character(*), parameter,public :: cpl_fields_a2c_states = &
      &'Sa_z&
      &:Sa_u&
      &:Sa_v&
      &:Sa_tbot&
      &:Sa_ptem&
      &:Sa_shum&
      &:Sa_dens&
      &:Sa_pbot&
      &:Sa_pslv'
   character(*), parameter,public :: cpl_fields_a2c_fluxes = &
      &'Faxa_lwdn&
      &:Faxa_rainc&
      &:Faxa_rainl&
      &:Faxa_snowc&
      &:Faxa_snowl&
      &:Faxa_swndr&
      &:Faxa_swvdr&
      &:Faxa_swndf&
      &:Faxa_swvdf&
      &:Faxa_swnet'
   character(*), parameter,public :: cpl_fields_a2c_fields = &
      trim(cpl_fields_a2c_states)//":"//trim(cpl_fields_a2c_fluxes)

   !----- state fields -----
   integer(IN),parameter,public :: cpl_fields_a2c_z      =  1 ! bottom atm level height
   integer(IN),parameter,public :: cpl_fields_a2c_u      =  2 ! bottom atm level zon wind
   integer(IN),parameter,public :: cpl_fields_a2c_v      =  3 ! bottom atm level mer wind
   integer(IN),parameter,public :: cpl_fields_a2c_tbot   =  4 ! bottom atm level temp
   integer(IN),parameter,public :: cpl_fields_a2c_ptem   =  5 ! bottom atm level pot temp
   integer(IN),parameter,public :: cpl_fields_a2c_shum   =  6 ! bottom atm level spec hum
   integer(IN),parameter,public :: cpl_fields_a2c_dens   =  7 ! bottom atm level air den
   integer(IN),parameter,public :: cpl_fields_a2c_pbot   =  8 ! bottom atm level pressure
   integer(IN),parameter,public :: cpl_fields_a2c_pslv   =  9 ! sea level atm pressure
   !----- fluxes computed by atm ----
   integer(IN),parameter,public :: cpl_fields_a2c_lwdn   = 10 ! downward lw heat flux
   integer(IN),parameter,public :: cpl_fields_a2c_rainc  = 11 ! prec: liquid "convective"
   integer(IN),parameter,public :: cpl_fields_a2c_rainl  = 12 ! prec: liquid "large scale"
   integer(IN),parameter,public :: cpl_fields_a2c_snowc  = 13 ! prec: frozen "convective"
   integer(IN),parameter,public :: cpl_fields_a2c_snowl  = 14 ! prec: frozen "large scale"
   integer(IN),parameter,public :: cpl_fields_a2c_swndr  = 15 ! sw: nir direct  downward
   integer(IN),parameter,public :: cpl_fields_a2c_swvdr  = 16 ! sw: vis direct  downward
   integer(IN),parameter,public :: cpl_fields_a2c_swndf  = 17 ! sw: nir diffuse downward
   integer(IN),parameter,public :: cpl_fields_a2c_swvdf  = 18 ! sw: vis diffuse downward
   integer(IN),parameter,public :: cpl_fields_a2c_swnet  = 19 ! sw: net


   integer(IN),parameter,public :: cpl_fields_c2a_total  = 17

   character(*), parameter,public :: cpl_fields_c2a_states = &
      &'Sx_tref&
      &:Sx_qref&
      &:Sx_avsdr&
      &:Sx_anidr&
      &:Sx_avsdf&
      &:Sx_anidf&
      &:Sx_t&
      &:So_t&
      &:Sx_snowh&
      &:Sx_ifrac&
      &:Sx_ofrac'
   character(*), parameter,public :: cpl_fields_c2a_fluxes = &
      &'Faxx_taux&
      &:Faxx_tauy&
      &:Faxx_lat&
      &:Faxx_sen&
      &:Faxx_lwup&
      &:Faxx_evap'
   character(*), parameter,public :: cpl_fields_c2a_fields = &
      trim(cpl_fields_c2a_states)//":"//trim(cpl_fields_c2a_fluxes)

   !----- states given to atm ----
   integer(IN),parameter,public :: cpl_fields_c2a_tref  =  1 ! 2m reference temperature
   integer(IN),parameter,public :: cpl_fields_c2a_qref  =  2 ! 2m reference specific humidity
   integer(IN),parameter,public :: cpl_fields_c2a_avsdr =  3 ! albedo, visible, direct
   integer(IN),parameter,public :: cpl_fields_c2a_anidr =  4 ! albedo, near-ir, direct
   integer(IN),parameter,public :: cpl_fields_c2a_avsdf =  5 ! albedo, visible, diffuse
   integer(IN),parameter,public :: cpl_fields_c2a_anidf =  6 ! albedo, near-ir, diffuse
   integer(IN),parameter,public :: cpl_fields_c2a_t     =  7 ! surface temperature
   integer(IN),parameter,public :: cpl_fields_c2a_sst   =  8 ! sea surface temperature
   integer(IN),parameter,public :: cpl_fields_c2a_snowh =  9 ! surface snow depth
   integer(IN),parameter,public :: cpl_fields_c2a_ifrac = 10 ! surface ice fraction
   integer(IN),parameter,public :: cpl_fields_c2a_ofrac = 11 ! surface ocn fraction
   !----- fluxes given to atm ----
   integer(IN),parameter,public :: cpl_fields_c2a_taux  = 12 ! wind stress, zonal
   integer(IN),parameter,public :: cpl_fields_c2a_tauy  = 13 ! wind stress, meridional
   integer(IN),parameter,public :: cpl_fields_c2a_lat   = 14 ! latent          heat flux
   integer(IN),parameter,public :: cpl_fields_c2a_sen   = 15 ! sensible        heat flux
   integer(IN),parameter,public :: cpl_fields_c2a_lwup  = 16 ! upward longwave heat flux
   integer(IN),parameter,public :: cpl_fields_c2a_evap  = 17 ! evaporation    water flux

   !----------------------------------------------------------------------------
   ! ice fields
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: cpl_fields_i2c_total  = 22

   character(*), parameter,public :: cpl_fields_i2c_states = &
      &'Si_t&
      &:Si_tref&
      &:Si_qref&
      &:Si_ifrac&
      &:Si_avsdr&
      &:Si_anidr&
      &:Si_avsdf&
      &:Si_anidf&
      &:index'
   character(*), parameter,public :: cpl_fields_i2c_fluxes = &
      &'Faii_taux&
      &:Faii_tauy&
      &:Faii_lat&
      &:Faii_sen&
      &:Faii_lwup&
      &:Faii_evap&
      &:Faii_swnet&
      &:Fioi_swpen&
      &:Fioi_melth&
      &:Fioi_meltw&
      &:Fioi_salt&
      &:Fioi_taux&
      &:Fioi_tauy'
   character(*), parameter,public :: cpl_fields_i2c_fields = &
      trim(cpl_fields_i2c_states)//":"//trim(cpl_fields_i2c_fluxes)

   !----- ice states -----
   integer(IN),parameter,public :: cpl_fields_i2c_t     =  1 ! temperature
   integer(IN),parameter,public :: cpl_fields_i2c_tref  =  2 ! 2m reference temperature
   integer(IN),parameter,public :: cpl_fields_i2c_qref  =  3 ! 2m reference specific humidity
   integer(IN),parameter,public :: cpl_fields_i2c_ifrac =  4 ! fractional ice coverage
   integer(IN),parameter,public :: cpl_fields_i2c_avsdr =  5 ! albedo: visible, direct
   integer(IN),parameter,public :: cpl_fields_i2c_anidr =  6 ! albedo: near ir, direct
   integer(IN),parameter,public :: cpl_fields_i2c_avsdf =  7 ! albedo: visible, diffuse
   integer(IN),parameter,public :: cpl_fields_i2c_anidf =  8 ! albedo: near ir, diffuse
   !----- compression index -----
   integer(IN),parameter,public :: cpl_fields_i2c_index =  9 ! global data compr index
   !----- a/i fluxes computed by ice -----
   integer(IN),parameter,public :: cpl_fields_i2c_taux  = 10 ! wind stress, zonal
   integer(IN),parameter,public :: cpl_fields_i2c_tauy  = 11 ! wind stress, meridional
   integer(IN),parameter,public :: cpl_fields_i2c_lat   = 12 ! latent          heat flux
   integer(IN),parameter,public :: cpl_fields_i2c_sen   = 13 ! sensible        heat flux
   integer(IN),parameter,public :: cpl_fields_i2c_lwup  = 14 ! upward longwave heat flux
   integer(IN),parameter,public :: cpl_fields_i2c_evap  = 15 ! evaporation    water flux
   integer(IN),parameter,public :: cpl_fields_i2c_swnet = 16 ! shortwave: net absorbed
   !----- i/o fluxes computed by ice -----
   integer(IN),parameter,public :: cpl_fields_i2c_swpen = 17 ! net SW penetrating ice
   integer(IN),parameter,public :: cpl_fields_i2c_melth = 18 ! heat  flux from melting ice
   integer(IN),parameter,public :: cpl_fields_i2c_meltw = 19 ! water flux from melting ice
   integer(IN),parameter,public :: cpl_fields_i2c_salt  = 20 ! salt  flux from melting ice
   integer(IN),parameter,public :: cpl_fields_i2c_otaux = 21 ! ice/ocn stress, zonal
   integer(IN),parameter,public :: cpl_fields_i2c_otauy = 22 ! ice/ocn stress, meridional


   integer(IN),parameter,public :: cpl_fields_c2i_total  = 21

   character(*), parameter,public :: cpl_fields_c2i_states = &
      &'So_t&
      &:So_s&
      &:So_u&
      &:So_v&
      &:Sa_z&
      &:Sa_u&
      &:Sa_v&
      &:Sa_ptem&
      &:Sa_tbot&
      &:Sa_shum&
      &:Sa_dens&
      &:So_dhdx&
      &:So_dhdy'
   character(*), parameter,public :: cpl_fields_c2i_fluxes = &
      &'Fioo_q&
      &:Faxa_swndr&
      &:Faxa_swvdr&
      &:Faxa_swndf&
      &:Faxa_swvdf&
      &:Faxa_lwdn&
      &:Faxc_rain&
      &:Faxc_snow'
   character(*), parameter,public :: cpl_fields_c2i_fields = &
      trim(cpl_fields_c2i_states)//":"//trim(cpl_fields_c2i_fluxes)

   !----- ocn states -----
   integer(IN),parameter,public :: cpl_fields_c2i_ot    =  1 ! ocn temp
   integer(IN),parameter,public :: cpl_fields_c2i_os    =  2 ! ocn salinity
   integer(IN),parameter,public :: cpl_fields_c2i_ou    =  3 ! ocn u velocity
   integer(IN),parameter,public :: cpl_fields_c2i_ov    =  4 ! ocn v velocity
   integer(IN),parameter,public :: cpl_fields_c2i_dhdx  = 12 ! ocn surface slope, zonal
   integer(IN),parameter,public :: cpl_fields_c2i_dhdy  = 13 ! ocn surface slope, merid
   !----- atm states ----- 
   integer(IN),parameter,public :: cpl_fields_c2i_z     =  5 ! atm bottom layer height
   integer(IN),parameter,public :: cpl_fields_c2i_u     =  6 ! atm u velocity
   integer(IN),parameter,public :: cpl_fields_c2i_v     =  7 ! atm v velocity
   integer(IN),parameter,public :: cpl_fields_c2i_ptem  =  8 ! atm potential temp
   integer(IN),parameter,public :: cpl_fields_c2i_tbot  =  9 ! atm bottom temp
   integer(IN),parameter,public :: cpl_fields_c2i_shum  = 10 ! atm specfic humidity
   integer(IN),parameter,public :: cpl_fields_c2i_dens  = 11 ! atm air density
   !----- ocn fluxes -----
   integer(IN),parameter,public :: cpl_fields_c2i_q     = 14 ! ocn freeze or melt heat
   !----- atm fluxes -----
   integer(IN),parameter,public :: cpl_fields_c2i_swndr = 15 ! atm sw near-ir, direct
   integer(IN),parameter,public :: cpl_fields_c2i_swvdr = 16 ! atm sw visable, direct
   integer(IN),parameter,public :: cpl_fields_c2i_swndf = 17 ! atm sw near-ir, diffuse
   integer(IN),parameter,public :: cpl_fields_c2i_swvdf = 18 ! atm sw visable, diffuse
   integer(IN),parameter,public :: cpl_fields_c2i_lwdn  = 19 ! long-wave down
   integer(IN),parameter,public :: cpl_fields_c2i_rain  = 20 ! rain
   integer(IN),parameter,public :: cpl_fields_c2i_snow  = 21 ! snow

   !----------------------------------------------------------------------------
   ! lnd fields 
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: cpl_fields_l2c_total  = 15

   character(*), parameter,public :: cpl_fields_l2c_states = &
      &'Sl_t&
      &:Sl_tref&
      &:Sl_qref&
      &:Sl_avsdr&
      &:Sl_anidr&
      &:Sl_avsdf&
      &:Sl_anidf&
      &:Sl_snowh'
   character(*), parameter,public :: cpl_fields_l2c_fluxes = &
      &'Fall_taux&
      &:Fall_tauy&
      &:Fall_lat&
      &:Fall_sen&
      &:Fall_lwup&
      &:Fall_evap&
      &:Fall_swnet'
   character(*), parameter,public :: cpl_fields_l2c_fields = &
      trim(cpl_fields_l2c_states)//":"//trim(cpl_fields_l2c_fluxes)

   !----- lnd states -----
   integer(IN),parameter,public :: cpl_fields_l2c_t     =  1 ! temperature
   integer(IN),parameter,public :: cpl_fields_l2c_tref  =  2 ! 2m reference temperature
   integer(IN),parameter,public :: cpl_fields_l2c_qref  =  3 ! 2m reference specific humidity
   integer(IN),parameter,public :: cpl_fields_l2c_avsdr =  4 ! albedo: direct , visible
   integer(IN),parameter,public :: cpl_fields_l2c_anidr =  5 ! albedo: direct , near-ir
   integer(IN),parameter,public :: cpl_fields_l2c_avsdf =  6 ! albedo: diffuse, visible
   integer(IN),parameter,public :: cpl_fields_l2c_anidf =  7 ! albedo: diffuse, near-ir
   integer(IN),parameter,public :: cpl_fields_l2c_snowh =  8 ! snow height
   !----- computed by lnd -----
   integer(IN),parameter,public :: cpl_fields_l2c_taux  =  9 ! wind stress, zonal
   integer(IN),parameter,public :: cpl_fields_l2c_tauy  = 10 ! wind stress, meridional
   integer(IN),parameter,public :: cpl_fields_l2c_lat   = 11 ! latent          heat flux
   integer(IN),parameter,public :: cpl_fields_l2c_sen   = 12 ! sensible        heat flux
   integer(IN),parameter,public :: cpl_fields_l2c_lwup  = 13 ! upward longwave heat flux
   integer(IN),parameter,public :: cpl_fields_l2c_evap  = 14 ! evaporation    water flux
   integer(IN),parameter,public :: cpl_fields_l2c_swnet = 15 ! 2m reference temperature


   integer(IN),parameter,public :: cpl_fields_c2l_total  = 18

   character(*), parameter,public :: cpl_fields_c2l_states = &
      &'Sa_z&
      &:Sa_u&
      &:Sa_v&
      &:Sa_tbot&
      &:Sa_ptem&
      &:Sa_shum&
      &:Sa_dens&
      &:Sa_pbot&
      &:Sa_pslv'
   character(*), parameter,public :: cpl_fields_c2l_fluxes = &
      &'Faxa_lwdn&
      &:Faxa_rainc&
      &:Faxa_rainl&
      &:Faxa_snowc&
      &:Faxa_snowl&
      &:Faxa_swndr&
      &:Faxa_swvdr&
      &:Faxa_swndf&
      &:Faxa_swvdf'
   character(*), parameter,public :: cpl_fields_c2l_fields = &
      trim(cpl_fields_c2l_states)//":"//trim(cpl_fields_c2l_fluxes)

   !----- atm states -----
   integer(IN),parameter,public :: cpl_fields_c2l_z     =  1 ! bottom atm level height
   integer(IN),parameter,public :: cpl_fields_c2l_u     =  2 ! bottom atm level zon wind
   integer(IN),parameter,public :: cpl_fields_c2l_v     =  3 ! bottom atm level mer wind
   integer(IN),parameter,public :: cpl_fields_c2l_tbot  =  4 ! bottom atm level temp
   integer(IN),parameter,public :: cpl_fields_c2l_ptem  =  5 ! bottom atm level pot temp
   integer(IN),parameter,public :: cpl_fields_c2l_shum  =  6 ! bottom atm level spec hum
   integer(IN),parameter,public :: cpl_fields_c2l_dens  =  7 ! bottom atm level air dens
   integer(IN),parameter,public :: cpl_fields_c2l_pbot  =  8 ! bottom atm level pressure
   integer(IN),parameter,public :: cpl_fields_c2l_pslv  =  9 ! sea level atm pressure
   !----- computed by atm -----
   integer(IN),parameter,public :: cpl_fields_c2l_lwdn  = 10 ! downward longwave heat flux
   integer(IN),parameter,public :: cpl_fields_c2l_rainc = 11 ! precip: liquid, convective
   integer(IN),parameter,public :: cpl_fields_c2l_rainl = 12 ! precip: liquid, large-scale
   integer(IN),parameter,public :: cpl_fields_c2l_snowc = 13 ! precip: frozen, convective
   integer(IN),parameter,public :: cpl_fields_c2l_snowl = 14 ! precip: frozen, large-scale
   integer(IN),parameter,public :: cpl_fields_c2l_swndr = 15 ! shortwave: nir direct  down
   integer(IN),parameter,public :: cpl_fields_c2l_swvdr = 16 ! shortwave: vis direct  down
   integer(IN),parameter,public :: cpl_fields_c2l_swndf = 17 ! shortwave: nir diffuse down
   integer(IN),parameter,public :: cpl_fields_c2l_swvdf = 18 ! shortwave: vis diffuse down

   !----- special lnd grid initialization -----

   integer(IN),parameter,public :: cpl_fields_c2lg_total  = 6

   character(*), parameter,public :: cpl_fields_c2lg_fields = &
      &'lon&
      &:lat&
      &:area&
      &:lfrac&
      &:maskl&
      &:maska'

   integer(IN),parameter,public :: cpl_fields_c2lg_alon  =  1 ! longitude
   integer(IN),parameter,public :: cpl_fields_c2lg_alat  =  2 ! latitude
   integer(IN),parameter,public :: cpl_fields_c2lg_aarea =  3 ! cell area
   integer(IN),parameter,public :: cpl_fields_c2lg_lfrac =  4 ! lnd fraction
   integer(IN),parameter,public :: cpl_fields_c2lg_lmask =  5 ! lnd mask
   integer(IN),parameter,public :: cpl_fields_c2lg_amask =  6 ! atm mask

   !----------------------------------------------------------------------------
   ! ocn fields
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: cpl_fields_o2c_total  = 7

   character(*), parameter,public :: cpl_fields_o2c_states = &
      &'So_t&
      &:So_u&
      &:So_v&
      &:So_s&
      &:So_dhdx&
      &:So_dhdy'
   character(*), parameter,public :: cpl_fields_o2c_fluxes = &
      &'Fioo_q'
   character(*), parameter,public :: cpl_fields_o2c_fields = &
      trim(cpl_fields_o2c_states)//":"//trim(cpl_fields_o2c_fluxes)


   !----- ocn states -----
   integer(IN),parameter,public :: cpl_fields_o2c_t    =  1 ! temperature
   integer(IN),parameter,public :: cpl_fields_o2c_u    =  2 ! velocity, zonal
   integer(IN),parameter,public :: cpl_fields_o2c_v    =  3 ! velocity, meridional
   integer(IN),parameter,public :: cpl_fields_o2c_s    =  4 ! salinity
   integer(IN),parameter,public :: cpl_fields_o2c_dhdx =  5 ! surface slope, zonal
   integer(IN),parameter,public :: cpl_fields_o2c_dhdy =  6 ! surface slope, meridional
   integer(IN),parameter,public :: cpl_fields_o2c_q    =  7 ! heat of fusion (q>0) melt pot (q<0)


   integer(IN),parameter,public :: cpl_fields_c2o_total  = 18

   character(*), parameter,public :: cpl_fields_c2o_states = &
      &'Si_ifrac&
      &:Sa_pslv&
      &:Faoc_duu10n'
   character(*), parameter,public :: cpl_fields_c2o_fluxes = &
      &'Foxx_taux&
      &:Foxx_tauy&
      &:Foxx_swnet&
      &:Foxx_lat&
      &:Foxx_sen&
      &:Foxx_lwup&
      &:Foxx_lwdn&
      &:Foxx_melth&
      &:Foxx_salt&
      &:Foxx_prec&
      &:Foxx_snow&
      &:Foxx_rain&
      &:Foxx_evap&
      &:Foxx_meltw&
      &:Forr_roff'
   character(*), parameter,public :: cpl_fields_c2o_fields = &
      trim(cpl_fields_c2o_states)//":"//trim(cpl_fields_c2o_fluxes)

   !----- ocn model input -----
   integer(IN),parameter,public :: cpl_fields_c2o_ifrac =  1 ! state: ice fraction
   integer(IN),parameter,public :: cpl_fields_c2o_press =  2 ! state: sea level pressure
   integer(IN),parameter,public :: cpl_fields_c2o_duu10 =  3 ! state: 10m wind speed squared
   integer(IN),parameter,public :: cpl_fields_c2o_taux  =  4 ! wind stress: zonal
   integer(IN),parameter,public :: cpl_fields_c2o_tauy  =  5 ! wind stress: meridional
   integer(IN),parameter,public :: cpl_fields_c2o_swnet =  6 ! heat flux: shortwave net
   integer(IN),parameter,public :: cpl_fields_c2o_lat   =  7 ! heat flux: latent
   integer(IN),parameter,public :: cpl_fields_c2o_sen   =  8 ! heat flux: sensible
   integer(IN),parameter,public :: cpl_fields_c2o_lwup  =  9 ! heat flux: long-wave up
   integer(IN),parameter,public :: cpl_fields_c2o_lwdn  = 10 ! heat flux: long-wave down
   integer(IN),parameter,public :: cpl_fields_c2o_melth = 11 ! heat flux: melt
   integer(IN),parameter,public :: cpl_fields_c2o_salt  = 12 ! salt flux
   integer(IN),parameter,public :: cpl_fields_c2o_prec  = 13 ! water flux: rain+snow
   integer(IN),parameter,public :: cpl_fields_c2o_snow  = 14 ! water flux: snow
   integer(IN),parameter,public :: cpl_fields_c2o_rain  = 15 ! water flux: rain
   integer(IN),parameter,public :: cpl_fields_c2o_evap  = 16 ! water flux: evap
   integer(IN),parameter,public :: cpl_fields_c2o_meltw = 17 ! water flux: melt
   integer(IN),parameter,public :: cpl_fields_c2o_roff  = 18 ! water flux: runoff

   !----------------------------------------------------------------------------
   ! run-off field 
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: cpl_fields_r2c_total  = 1

   character(*), parameter,public :: cpl_fields_r2c_states = &
      &''
   character(*), parameter,public :: cpl_fields_r2c_fluxes = &
      &'Forr_roff'
!  character(*), parameter,public :: cpl_fields_r2c_fields = &
!    trim(cpl_fields_r2c_states)//":"//trim(cpl_fields_r2c_fluxes)
   character(*), parameter,public :: cpl_fields_r2c_fields = &
                                        trim(cpl_fields_r2c_fluxes)

   integer(IN),parameter,public :: cpl_fields_r2c_runoff = 1

!EOP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_fields_getField
!
! !DESCRIPTION:
!     Returns the field character string based on colon delimited
!     string and field number
!
! !REVISION HISTORY:
!     2003-Jan-24  - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine cpl_fields_getField(outfield,nfld,cstring)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(out) :: outfield   ! output field name
   integer     ,intent(in ) :: nfld       ! field number
   character(*),intent(in ) :: cstring    ! colon delimited field string

!EOP

  type(cpl_mct_list)   :: mctIstr  ! mct list from input cstring
  type(cpl_mct_string) :: mctOStr  ! mct string for output outfield

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  outfield = ' '

  call cpl_mct_list_init(mctIstr,cstring)
  call cpl_mct_list_get(mctOStr,nfld,mctIstr)
  outfield = cpl_mct_string_toChar(mctOStr)
  call cpl_mct_list_clean(mctIstr)
  call cpl_mct_string_clean(mctOStr)

end subroutine cpl_fields_getField

!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_fields_getLongName -- get netCDF attributes for a field
!
! !DESCRIPTION:
!    By parsing a field name and using a lookup table, get the netCDF attribute 
!    character strings corresponding to the given field name.
!
! !REVISION HISTORY:
!     2003-may-12 - B. Kauffman - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_fields_getLongName(fldstr,longname,units)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)   :: fldstr     ! field name
   character(*),intent(out)  :: longname   ! corresponding longname
   character(*),intent(out)  :: units      ! corresponding units

!EOP

   !--- local ---
   integer,parameter    :: nl = 60             ! max number of longnames
   integer              :: n                   ! generic loop index
   character(32),save   :: lookup(nl,3)        ! longname lookup table
   character(32)        :: shortname           ! short name
   logical, save        :: firstCall = .true.  ! flags initializtion

   !----- formats -----
   character(*),parameter :: subname = "(cpl_fields_getLongName)"
   character(*),parameter :: F00 = "('(cpl_fields_getLongName) ',4a)"

!-------------------------------------------------------------------------------

   if (FirstCall) then
      FirstCall = .false.
      lookup( :,:) = "undefined"
      lookup( 1,1) = "area"  ; lookup( 1,2) = "cell area                     " ; lookup( 1,3) = "rad^2   "
      lookup( 2,1) = "aream" ; lookup( 2,2) = "cell area from map file       " ; lookup( 2,3) = "rad^2   "
      lookup( 3,1) = "index" ; lookup( 3,2) = "global index                  " ; lookup( 3,3) = "unitless"
      lookup( 4,1) = "lat"   ; lookup( 4,2) = "latitude                      " ; lookup( 4,3) = "degrees north "
      lookup( 5,1) = "lon"   ; lookup( 5,2) = "longitude                     " ; lookup( 5,3) = "degrees east  "
      lookup( 6,1) = "mask"  ; lookup( 6,2) = "domain mask                   " ; lookup( 6,3) = "unitless"
      lookup( 7,1) = "maska" ; lookup( 7,2) = "domain mask, atm              " ; lookup( 7,3) = "unitless"
      lookup( 8,1) = "maski" ; lookup( 8,2) = "domain mask, ice              " ; lookup( 8,3) = "unitless"
      lookup( 9,1) = "maskl" ; lookup( 9,2) = "domain mask, lnd              " ; lookup( 9,3) = "unitless"
      lookup(10,1) = "masko" ; lookup(10,2) = "domain mask, ocn              " ; lookup(10,3) = "unitless"
      lookup(11,1) = "pid"   ; lookup(11,2) = "process ID                    " ; lookup(11,3) = "unitless"

      lookup(12,1) = "anidf" ; lookup(12,2) = "albedo, near-infrared, diffuse" ; lookup(12,3) = "unitless"
      lookup(13,1) = "anidr" ; lookup(13,2) = "albedo, near-infrared, direct " ; lookup(13,3) = "unitless"
      lookup(14,1) = "avsdf" ; lookup(14,2) = "albedo, visible, diffuse      " ; lookup(14,3) = "unitless"
      lookup(15,1) = "avsdr" ; lookup(15,2) = "albedo, visible, direct       " ; lookup(15,3) = "unitless"
      lookup(16,1) = "dens"  ; lookup(16,2) = "density                       " ; lookup(16,3) = "kg/m^3  "
      lookup(17,1) = "dhdx"  ; lookup(17,2) = "surface slope, zonal          " ; lookup(17,3) = "m/m     "
      lookup(18,1) = "dhdy"  ; lookup(18,2) = "surface slope, meridional     " ; lookup(18,3) = "m/m     "
      lookup(19,1) = "duu10n"; lookup(19,2) = "10m neutral wind speed squared" ; lookup(19,3) = "m^2/s^2 "
      lookup(20,1) = "evap"  ; lookup(20,2) = "evaporation                   " ; lookup(20,3) = "kg/s/m^2"
      lookup(21,1) = "frac"  ; lookup(21,2) = "fraction                      " ; lookup(21,3) = "unitless"
      lookup(22,1) = "afrac" ; lookup(22,2) = "fraction atm                  " ; lookup(22,3) = "unitless"
      lookup(23,1) = "ifrac" ; lookup(23,2) = "fraction ice                  " ; lookup(23,3) = "unitless"
      lookup(24,1) = "lfrac" ; lookup(24,2) = "fraction lnd                  " ; lookup(24,3) = "unitless"
      lookup(25,1) = "oraco" ; lookup(25,2) = "fraction ocn                  " ; lookup(25,3) = "unitless"
      lookup(26,1) = "lwdn"  ; lookup(26,2) = "longwave radiation, upward    " ; lookup(26,3) = "W/m^2   "
      lookup(27,1) = "lwup"  ; lookup(27,2) = "longwave radiation, downward  " ; lookup(27,3) = "W/m^2   "
      lookup(28,1) = "melth" ; lookup(28,2) = "melt heat                     " ; lookup(28,3) = "W/m^2   "
      lookup(29,1) = "meltw" ; lookup(29,2) = "melt water                    " ; lookup(29,3) = "kg/s/m^2"
      lookup(30,1) = "pbot"  ; lookup(30,2) = "pressure, bottom              " ; lookup(30,3) = "Pa      "
      lookup(31,1) = "prec"  ; lookup(31,2) = "precipitation                 " ; lookup(31,3) = "kg/s/m^2"
      lookup(32,1) = "pslv"  ; lookup(32,2) = "pressure, sea level           " ; lookup(32,3) = "Pa      "
      lookup(33,1) = "ptem"  ; lookup(33,2) = "potential temperature         " ; lookup(33,3) = "kelvin  "
      lookup(34,1) = "q"     ; lookup(34,2) = "q<0 = heat of fusion, &
                                              &q>0 = melting potential"        ; lookup(34,3) = "W/m^2   "
      lookup(35,1) = "qref"  ; lookup(35,2) = "humidity, reference           " ; lookup(35,3) = "kg/kg   "
      lookup(36,1) = "rain"  ; lookup(36,2) = "precipitation, liquid         " ; lookup(36,3) = "kg/s/m^2"
      lookup(37,1) = "rainc" ; lookup(37,2) = "precip, liquid, convective    " ; lookup(37,3) = "kg/s/m^2"
      lookup(38,1) = "rainl" ; lookup(38,2) = "precip, liquid, large-scale   " ; lookup(38,3) = "kg/s/m^2"
      lookup(39,1) = "roff"  ; lookup(39,2) = "water flux: runoff            " ; lookup(39,3) = "kg/s/m^2"
      lookup(40,1) = "salt"  ; lookup(40,2) = "salt flux                     " ; lookup(40,3) = "kg/s/m^2"
      lookup(41,1) = "sen"   ; lookup(41,2) = "sensible heat flux            " ; lookup(41,3) = "W/m^2   "
      lookup(42,1) = "shum"  ; lookup(42,2) = "humidity, specific            " ; lookup(42,3) = "kg/kg   "
      lookup(43,1) = "snow"  ; lookup(43,2) = "precipitation, frozen         " ; lookup(43,3) = "kg/s/m^2"
      lookup(44,1) = "snowc" ; lookup(44,2) = "precip, frozen, convective    " ; lookup(44,3) = "kg/s/m^2"
      lookup(45,1) = "snowl" ; lookup(45,2) = "precip, frozen, large-scale   " ; lookup(45,3) = "kg/s/m^2"
      lookup(46,1) = "swndf" ; lookup(46,2) = "rad, sw, near-infr, diffuse   " ; lookup(46,3) = "W/m^2   "
      lookup(47,1) = "swndr" ; lookup(47,2) = "rad, sw, near-infr, direct    " ; lookup(47,3) = "W/m^2   "
      lookup(48,1) = "swnet" ; lookup(48,2) = "rad, sw, net                  " ; lookup(48,3) = "W/m^2   "
      lookup(49,1) = "swpen" ; lookup(49,2) = "rad, sw, penetrating          " ; lookup(49,3) = "W/m^2   "
      lookup(50,1) = "swvdf" ; lookup(50,2) = "rad, sw, visible, diffuse     " ; lookup(50,3) = "W/m^2   "
      lookup(51,1) = "swvdr" ; lookup(51,2) = "rad, sw, visible, direct      " ; lookup(51,3) = "W/m^2   "
      lookup(52,1) = "t"     ; lookup(52,2) = "temperature                   " ; lookup(52,3) = "kelvin  "
      lookup(53,1) = "taux"  ; lookup(53,2) = "stress, zonal                 " ; lookup(53,3) = "N/m^2   "
      lookup(54,1) = "tauy"  ; lookup(54,2) = "stress, meridional            " ; lookup(54,3) = "N/m^2   "
      lookup(55,1) = "tbot"  ; lookup(55,2) = "temperature, bottom           " ; lookup(55,3) = "kelvin  "
      lookup(56,1) = "tref"  ; lookup(56,2) = "temperature, reference        " ; lookup(56,3) = "kelvin  "
      lookup(57,1) = "u"     ; lookup(57,2) = "velocity, zonal               " ; lookup(57,3) = "N/m^2   "
      lookup(58,1) = "v"     ; lookup(58,2) = "velocity, meridional          " ; lookup(58,3) = "N/m^2   "
      lookup(59,1) = "z"     ; lookup(59,2) = "height                        " ; lookup(59,3) = "m       "
      lookup(60,1) = "snowh" ; lookup(60,2) = "snow depth                    " ; lookup(60,3) = "m       "
   end if

   !--- find shortname in suffix of field name string ---
   n = len_trim(fldstr)
   do while ( n>0 )
      if ( fldstr(n:n) == "_") exit
      n = n - 1
   end do
   shortname = ""
   if (n < len_trim(fldstr)) shortname = fldstr(n+1:len_trim(fldstr))
   if (len_trim(shortname) < 1) write(6,F00) "WARNING: short name has length 0"

   !--- find/set longname & units corresponding to shortname ---
   longname = "undefined"
   units    = "undefined"
   do n=1,nl
      if ( trim(shortname) == trim(lookup(n,1)) ) then
         longname = trim(lookup(n,2))
         units    = trim(lookup(n,3))
         exit
      end if
   end do

end subroutine cpl_fields_getLongName

!===============================================================================

end module cpl_fields_mod

