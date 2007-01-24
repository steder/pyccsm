!===============================================================================
! CVS $Id: merge_mod.F90,v 1.1.1.1 2005/02/03 22:29:00 steder Exp $
! CVS $Source: /home/cvsroot/steder/pyCPL/merge_mod.F90,v $
! CVS $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: merge_mod -- field merging module.
!
! !DESCRIPTION:
!    Merges fields to be sent to a component.  "Merging" means combining one or
!    more fields to create a new field.  Typically this is two or more fields
!    of the same type, e.g. combining atm/ice, atm/lnd, and atm/ocn sensible heat
!    flux fields, weighted by surface fraction, to create an atm/surface heat 
!    flux.  But it could also involve somewhat differing fields, e.g. 
!    precipitation plus snow melt.  Merging is normally not an automatic process, 
!    some hand-tuning is generally necessary to achieve the results necessary for 
!    valid science.
!
! !REVISION HISTORY:
!     2002-Sep-27 - B. Kauffman - created initial version
!
! !INTERFACE: ------------------------------------------------------------------

module merge_mod

! !USES:

   use cpl_kind_mod        ! kinds
   use cpl_control_mod     ! control flags
   use cpl_bundle_mod      ! bundle data type and methods
   use frac_mod            ! surface fractions
   use data_mod            ! lengthly data declarations/inits for main program
   use shr_sys_mod         ! wrappers to system calls
   use shr_timer_mod       ! timing utilities

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: merge_atm        ! merges atm input data
   public :: merge_ocn        ! merges ocn input data

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: merge_atm -- merge bundles to form atm input bundle
!
! !DESCRIPTION:
!    merge bundles to form atm input bundle
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2002-Jun-09 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine merge_atm(fix_So2c_a)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in) :: fix_So2c_a ! KLUDGE: use alt bun_So2a_a, temp fix

!EOP

   !----- local -----
   integer(IN),save :: t01 = -1 ! timer id
   integer(IN),save :: t02 = -1 ! timer id
   integer(IN),save :: t03 = -1 ! timer id
   integer(IN)      :: j,k      ! So_t index, wrt KLUDGE fix wrt mapping bug
   
   !----- formats -----
   character(*),parameter :: F00 = "('(merge_atm) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   if (t01 == -1) call shr_timer_get(t01,"merge_atm")
   if (t02 == -1) call shr_timer_get(t02,"merge_atm 0gather")
   if (t03 == -1) call shr_timer_get(t03,"merge_atm add")

   call shr_timer_start(t01)
   call shr_timer_start(t02)

   !--- merge atm states & fluxes ---
   call cpl_bundle_zero(con_Xc2a%bundle) ! zero before adding

   call cpl_bundle_gather(con_Xc2a%bundle, bun_Sa2c_a, bun_Fa2c_a, &
   &                                       bun_Sl2c_l, bun_Fl2c_l, &
   &                                       bun_So2c_a, bun_Fo2c_a, &
   &                                       bun_Si2c_a, bun_Fi2c_a )

   call shr_timer_stop(t02)
   call shr_timer_start(t03)

   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_taux',bun_aoflux_a,'Faoc_taux' &
                                          ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_taux',bun_Fl2c_l  ,'Fall_taux' &
                                          ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_taux',bun_Fi2c_a  ,'Faii_taux' &
                                          ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_tauy',bun_aoflux_a,'Faoc_tauy' &
                                          ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_tauy',bun_Fl2c_l  ,'Fall_tauy' &
                                          ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_tauy',bun_Fi2c_a  ,'Faii_tauy' &
                                          ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_lat',bun_aoflux_a,'Faoc_lat' &
                                         ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_lat',bun_Fl2c_l  ,'Fall_lat' &
                                         ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_lat',bun_Fi2c_a  ,'Faii_lat' &
                                         ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_sen',bun_aoflux_a,'Faoc_sen' &
                                         ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_sen',bun_Fl2c_l  ,'Fall_sen' &
                                         ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_sen',bun_Fi2c_a  ,'Faii_sen' &
                                         ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_lwup',bun_aoflux_a,'Faoc_lwup' &
                                          ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_lwup',bun_Fl2c_l  ,'Fall_lwup' &
                                          ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_lwup',bun_Fi2c_a  ,'Faii_lwup' &
                                          ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_evap',bun_aoflux_a,'Faoc_evap' &
                                          ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_evap',bun_Fl2c_l  ,'Fall_evap' &
                                          ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Faxx_evap',bun_Fi2c_a  ,'Faii_evap' &
                                          ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_tref',bun_aoflux_a,'Faoc_tref' &
                                          ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_tref',bun_Sl2c_l    ,'Sl_tref' &
                                          ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_tref',bun_Si2c_a    ,'Si_tref' &
                                          ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_qref',bun_aoflux_a,'Faoc_qref' &
                                          ,bun_frac_a  ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_qref',bun_Sl2c_l    ,'Sl_qref' &
                                          ,bun_frac_a  ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_qref',bun_Si2c_a    ,'Si_qref' &
                                          ,bun_frac_a  ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_avsdr',bun_oalbedo_a,'So_avsdr' &
                                           ,bun_frac_a   ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_avsdr',bun_Sl2c_l   ,'Sl_avsdr' &
                                           ,bun_frac_a   ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_avsdr',bun_Si2c_a   ,'Si_avsdr' &
                                           ,bun_frac_a   ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_anidr',bun_oalbedo_a,'So_anidr' &
                                           ,bun_frac_a   ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_anidr',bun_Sl2c_l   ,'Sl_anidr' &
                                           ,bun_frac_a   ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_anidr',bun_Si2c_a   ,'Si_anidr' &
                                           ,bun_frac_a   ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_avsdf',bun_oalbedo_a,'So_avsdf' &
                                           ,bun_frac_a   ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_avsdf',bun_Sl2c_l   ,'Sl_avsdf' &
                                           ,bun_frac_a   ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_avsdf',bun_Si2c_a   ,'Si_avsdf' &
                                           ,bun_frac_a   ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_anidf',bun_oalbedo_a,'So_anidf' &
                                           ,bun_frac_a   ,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_anidf',bun_Sl2c_l   ,'Sl_anidf' &
                                           ,bun_frac_a   ,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_anidf',bun_Si2c_a   ,'Si_anidf' &
                                           ,bun_frac_a   ,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_t',bun_So2c_a,'So_t' &
                                       ,bun_frac_a,'ofrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_t',bun_Sl2c_l,'Sl_t' &
                                       ,bun_frac_a,'lfrac')
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_t',bun_Si2c_a,'Si_t' &
                                       ,bun_frac_a,'ifrac')

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_snowh',bun_Sl2c_l,'Sl_snowh' &
                                           ,bun_frac_a,'lfrac', zero=.true.)

   call cpl_bundle_add(con_Xc2a%bundle,'Sx_ifrac',bun_frac_a,'ifrac', zero=.true.)
   call cpl_bundle_add(con_Xc2a%bundle,'Sx_ofrac',bun_frac_a,'ofrac', zero=.true.)

   !*** KLUDGE - BUG WORKAROUND ************************************************
   ! purpose: So_t values have errors due to mapping, can be less than freezing
   !****************************************************************************
!  --- FIX #1 : apply a floor on the So_t values where-ever ocn-fraction > 0
!  j  = cpl_mct_aVect_indexRA(con_Xc2a%bundle%data,'Sx_ofrac')      
!  k  = cpl_mct_aVect_indexRA(con_Xc2a%bundle%data,'So_t')
!   where (con_Xc2a%bundle%data%rAttr(j,:)  > 0.0_R8) &             
!      con_Xc2a%bundle%data%rAttr(k,:) = max(270.,con_Xc2a%bundle%data%rAttr(k,:)) 
!
!  --- FIX #2 : apply a floor on the So_t values everywhere
!  k  = cpl_mct_aVect_indexRA(con_Xc2a%bundle%data,'So_t')
!  con_Xc2a%bundle%data%rAttr(k,:) = max(270.,con_Xc2a%bundle%data%rAttr(k,:)) 
!
!  --- FIX #3 : use an alt calc of So_t from a bundle passed in as an arg
   k  = cpl_mct_aVect_indexRA(  con_Xc2a%bundle%data,'So_t')
   j  = cpl_mct_aVect_indexRA(       fix_So2c_a%data,'So_t')
   con_Xc2a%bundle%data%rAttr(k,:) = fix_So2c_a%data%rAttr(j,:)
   !****************************************************************************

   call shr_timer_stop(t03)
   call shr_timer_stop(t01)

end subroutine merge_atm

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: merge_ocn -- merge bundles to form ocn input bundle
!
! !DESCRIPTION:
!    merge bundles to form ocn input bundle
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2002-Jun-06 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine merge_ocn()

   implicit none

! !INPUT/OUTPUT PARAMETERS:

!EOP

   !----- local -----
   integer(IN),save :: t01 = -1 ! timer id
   integer(IN),save :: t02 = -1 ! timer id
   integer(IN),save :: t03 = -1 ! timer id
   
   integer(IN) :: n             ! loop index
   integer(IN) :: npts          ! size of AV
   real(R8)    :: afrac,ifrac   ! atm & ice fraction
   integer(IN) :: k_Foxx_taux   ! index into aVect
   integer(IN) :: k_Foxx_tauy
   integer(IN) :: k_Foxx_swnet
   integer(IN) :: k_Foxx_lat
   integer(IN) :: k_Foxx_sen
   integer(IN) :: k_Foxx_lwup
   integer(IN) :: k_Foxx_evap
   integer(IN) :: k_Foxx_lwdn
   integer(IN) :: k_Foxx_rain
   integer(IN) :: k_Foxx_snow
   integer(IN) :: k_Foxx_prec
   integer(IN) :: k_Foxx_melth
   integer(IN) :: k_Foxx_meltw
   integer(IN) :: k_Foxx_salt
   integer(IN) :: k_Si_ifrac
   integer(IN) :: k_Faoc_taux
   integer(IN) :: k_Faoc_tauy
   integer(IN) :: k_Faoc_swnet
   integer(IN) :: k_Faoc_lat
   integer(IN) :: k_Faoc_sen
   integer(IN) :: k_Faoc_lwup
   integer(IN) :: k_Faoc_evap
   integer(IN) :: k_Faxc_rain
   integer(IN) :: k_Faxc_snow
   integer(IN) :: k_Fioi_taux
   integer(IN) :: k_Fioi_tauy
   integer(IN) :: k_Fioi_swpen
   integer(IN) :: k_Faxa_lwdn
   integer(IN) :: k_Fioi_melth
   integer(IN) :: k_Fioi_meltw
   integer(IN) :: k_Fioi_salt
   integer(IN) :: k_afrac
   integer(IN) :: k_ifrac

   !----- formats -----
   character(*),parameter :: F00 = "('(merge_ocn) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   if (t01 == -1) call shr_timer_get(t01,"merge_ocn")
   if (t02 == -1) call shr_timer_get(t02,"merge_ocn 0gather")
   if (t03 == -1) call shr_timer_get(t03,"merge_ocn add")

   call shr_timer_start(t01)
   call shr_timer_start(t02)

   call cpl_bundle_zero(bun_Xc2oSNAP_o) ! zero before adding

   call cpl_bundle_gather(bun_Xc2oSNAP_o,bun_Sa2c_o, bun_Fa2c_o, &
!--- not allowed now ------------------- bun_Sl2c_o, bun_Fl2c_o, & --------
   &                                     bun_So2c_o, bun_Fo2c_o, &
   &                                     bun_Si2c_i, bun_Fi2c_i, &
   &                                     bun_Xr2c_o, bun_aoflux_o )

   call shr_timer_stop(t02)
   call shr_timer_start(t03)

   npts         = cpl_mct_aVect_lsize(bun_Xc2oSNAP_o%data)

   k_Foxx_taux  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_taux')
   k_Foxx_tauy  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_tauy')
   k_Foxx_swnet = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_swnet')
   k_Foxx_lat   = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_lat')
   k_Foxx_sen   = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_sen')
   k_Foxx_lwup  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_lwup')
   k_Foxx_evap  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_evap')
   k_Foxx_lwdn  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_lwdn')
   k_Foxx_rain  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_rain')
   k_Foxx_snow  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_snow')
   k_Foxx_prec  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_prec')
   k_Foxx_melth = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_melth')
   k_Foxx_meltw = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_meltw')
   k_Foxx_salt  = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Foxx_salt')
   k_Si_ifrac   = cpl_mct_aVect_indexRA(bun_Xc2oSNAP_o%data,'Si_ifrac')

   k_Faoc_taux  = cpl_mct_aVect_indexRA(bun_aoflux_o%data,'Faoc_taux')
   k_Faoc_tauy  = cpl_mct_aVect_indexRA(bun_aoflux_o%data,'Faoc_tauy')
   k_Faoc_swnet = cpl_mct_aVect_indexRA(bun_aoflux_o%data,'Faoc_swnet')
   k_Faoc_lat   = cpl_mct_aVect_indexRA(bun_aoflux_o%data,'Faoc_lat')
   k_Faoc_sen   = cpl_mct_aVect_indexRA(bun_aoflux_o%data,'Faoc_sen')
   k_Faoc_lwup  = cpl_mct_aVect_indexRA(bun_aoflux_o%data,'Faoc_lwup')
   k_Faoc_evap  = cpl_mct_aVect_indexRA(bun_aoflux_o%data,'Faoc_evap')
   k_Faxc_rain  = cpl_mct_aVect_indexRA(bun_precip_o%data,'Faxc_rain')
   k_Faxc_snow  = cpl_mct_aVect_indexRA(bun_precip_o%data,'Faxc_snow')

   k_Fioi_taux  = cpl_mct_aVect_indexRA(bun_Fi2c_i%data,'Fioi_taux')
   k_Fioi_tauy  = cpl_mct_aVect_indexRA(bun_Fi2c_i%data,'Fioi_tauy')
   k_Fioi_swpen = cpl_mct_aVect_indexRA(bun_Fi2c_i%data,'Fioi_swpen')
   k_Faxa_lwdn  = cpl_mct_aVect_indexRA(bun_Fa2c_o%data,'Faxa_lwdn')
   k_Fioi_melth = cpl_mct_aVect_indexRA(bun_Fi2c_i%data,'Fioi_melth')
   k_Fioi_meltw = cpl_mct_aVect_indexRA(bun_Fi2c_i%data,'Fioi_meltw')
   k_Fioi_salt  = cpl_mct_aVect_indexRA(bun_Fi2c_i%data,'Fioi_salt')

   k_afrac      = cpl_mct_aVect_indexRA(bun_frac_o%data,'afrac')
   k_ifrac      = cpl_mct_aVect_indexRA(bun_frac_o%data,'ifrac')

   do n = 1,npts
      afrac = bun_frac_o%data%rAttr(k_afrac,n)
      ifrac = bun_frac_o%data%rAttr(k_ifrac,n)

      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_taux,n) = &
        bun_aoflux_o%data%rAttr(k_Faoc_taux,n)*afrac + &
          bun_Fi2c_i%data%rAttr(k_Fioi_taux,n)*ifrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_tauy,n) = &
        bun_aoflux_o%data%rAttr(k_Faoc_tauy,n)*afrac + &
          bun_Fi2c_i%data%rAttr(k_Fioi_tauy,n)*ifrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_swnet,n) = &
        bun_aoflux_o%data%rAttr(k_Faoc_swnet,n)*afrac + &
          bun_Fi2c_i%data%rAttr(k_Fioi_swpen,n)*ifrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_lat,n) = &
        bun_aoflux_o%data%rAttr(k_Faoc_lat,n)*afrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_sen,n) = &
        bun_aoflux_o%data%rAttr(k_Faoc_sen,n)*afrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_lwup,n) = &
        bun_aoflux_o%data%rAttr(k_Faoc_lwup,n)*afrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_evap,n) = &
        bun_aoflux_o%data%rAttr(k_Faoc_evap,n)*afrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_lwdn,n) = &
          bun_Fa2c_o%data%rAttr(k_Faxa_lwdn,n)*afrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_rain,n) = &
        bun_precip_o%data%rAttr(k_Faxc_rain,n)*afrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_snow,n) = &
        bun_precip_o%data%rAttr(k_Faxc_snow,n)*afrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_prec,n) = &
         bun_Xc2oSNAP_o%data%rAttr(k_Foxx_rain,n) + &
         bun_Xc2oSNAP_o%data%rAttr(k_Foxx_snow,n) 
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_melth,n) = &
          bun_Fi2c_i%data%rAttr(k_Fioi_melth,n)*ifrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_meltw,n) = &
          bun_Fi2c_i%data%rAttr(k_Fioi_meltw,n)*ifrac
      bun_Xc2oSNAP_o%data%rAttr(k_Foxx_salt,n) = &
          bun_Fi2c_i%data%rAttr(k_Fioi_salt,n)*ifrac
      bun_Xc2oSNAP_o%data%rAttr(k_Si_ifrac,n) = ifrac

      bun_Xc2oSNAP_o%cnt = 1

   enddo

   call shr_timer_stop(t03)
   call shr_timer_stop(t01)

end subroutine merge_ocn
!===============================================================================

end module merge_mod
