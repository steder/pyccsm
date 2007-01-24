!===============================================================================
! CVS: $Id: data_mod.F90,v 1.1.1.1 2005/02/03 22:29:00 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/data_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: data_mod -- data declaration and initialazion for coupler main.
!
! !DESCRIPTION:
!  Does data declarations and initializations that might otherwise have
!  been located in the coupler main program.
! 
! !REMARKS:
!  o Bundles include an aVect (field data) and a pointer to an associated 
!    domain, thus several bundles can point to (share) the same domain.
!  o Domains include "local grid" info and domain decomposition info (gsMap), 
!    and also a pointer to an associated "global grid" (an un-decomposed grid),
!    thus several domains can point to (share) the same global grid (gGrid).
!  o Routers are specific to two domains (domains include/have a specific 
!    decomposition/gsMap): a "local" domain and a "remote" domain.  The router 
!    datatype does not keep this information, but it is implied and the cpl code
!    writer must remember which domain decompositions (gsMap's) are associated.
!  o ReArrangers are specific to two local domains: a source domain and a
!    destination domain.  The ReArranger datatype does not keep this 
!    information, but it is implied and the cpl code must keep track of this.
!    
!
! !REVISION HISTORY:
!     2002-May-xx - B. Kauffman - added bundleInit & mapInit routines
!     2002-Apr-28 - B. Kauffman - full set of declarations for CCSM cpl6.0
!     2001-Jun-08 - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

MODULE data_mod

! !USES:

   use cpl_mct_mod       ! access to  mct data types
   use cpl_domain_mod    ! defines domain data types
   use cpl_bundle_mod    ! defines bundle data types
   use cpl_map_mod       ! defines map    data types
   use cpl_fields_mod    ! indicies into bundles & ibuf
   use cpl_control_mod   ! control variables (eg. mapping file names)
   use cpl_kind_mod      ! kinds
   use shr_sys_mod       ! system call wrappers
   use cpl_contract_mod  ! contract

   implicit none

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: data_bundleInit  ! initialize the bundles declared in this module
   public :: data_mapInit     ! initialize the maps    declared in this module

! !PUBLIC DATA MEMBERS:

   !----------------------------------------------------------------------------
   ! datatypes (bundles & routers) for cpl/model communication, both ways
   !----------------------------------------------------------------------------

   !--- domains --- includes global grid + local grid and decomp info ---
   type(cpl_domain)     :: dom_a   ! atm  domain 
   type(cpl_domain)     :: dom_i   ! ice  domain
   type(cpl_domain)     :: dom_l   ! lnd  domain
   type(cpl_domain)     :: dom_r   ! roff domain
   type(cpl_domain)     :: dom_o   ! ocn  domain

   !--- bundles to/from component models ---

   type(cpl_contract) :: con_Xa2c  ! everything recv'd from atm
   type(cpl_contract) :: con_Xl2c  ! everything recv'd from lnd
   type(cpl_contract) :: con_Xr2c  ! everything recv'd from runoff
   type(cpl_contract) :: con_Xo2c  ! everything recv'd from ocn
   type(cpl_contract) :: con_Xi2c  ! everything recv'd from ice

   type(cpl_contract) :: con_Xc2a  ! everything sent   to   atm
   type(cpl_contract) :: con_Xc2l  ! everything sent   to   lnd
   type(cpl_contract) :: con_Dc2l  ! special grid info to   lnd
   type(cpl_contract) :: con_Xc2o  ! everything sent   to   ocn
   type(cpl_contract) :: con_Xc2i  ! everything sent   to   ice

   type(cpl_bundle) :: bun_Xc2oSNAP_o  ! everything sent to ocn, snapshot
   type(cpl_bundle) :: bun_Xc2oPSUM_o  ! everything sent to ocn, partial sum

   type(cpl_bundle) :: bun_aoflux_o  ! ao fluxes ocn grid
   type(cpl_bundle) :: bun_aoflux_a  ! ao fluxes atm grid
   type(cpl_bundle) :: bun_oalbedo_o ! ocean albedos on ocn grid
   type(cpl_bundle) :: bun_oalbedo_a ! ocean albedos on atm grid
   type(cpl_bundle) :: bun_precip_o  ! total snow and rain on ocn grid
   type(cpl_bundle) :: bun_precip_a  ! total snow and rain on atm grid

   type(cpl_bundle) :: bun_Sa2c_a  ! a2c states
   type(cpl_bundle) :: bun_Fa2c_a  ! a2c fluxes
   type(cpl_bundle) :: bun_Sa2c_o  ! a2c states mapped to o
   type(cpl_bundle) :: bun_Fa2c_o  ! a2c fluxes mapped to o
   type(cpl_bundle) :: bun_Sl2c_l  ! l2c states
   type(cpl_bundle) :: bun_Fl2c_l  ! l2c fluxes
!  type(cpl_bundle) :: bun_Sl2c_o  ! l2c states mapped to o
!  type(cpl_bundle) :: bun_Fl2c_o  ! l2c fluxes mapped to o
   type(cpl_bundle) :: bun_Xr2c_o  ! r2c fields mapped to o
   type(cpl_bundle) :: bun_So2c_o  ! o2c states
   type(cpl_bundle) :: bun_Fo2c_o  ! o2c fluxes
   type(cpl_bundle) :: bun_So2c_a  ! o2c states mapped to a
   type(cpl_bundle) :: bun_Fo2c_a  ! o2c fluxes mapped to a
   type(cpl_bundle) :: bun_Si2c_i  ! i2c states
   type(cpl_bundle) :: bun_Fi2c_i  ! i2c fluxes
   type(cpl_bundle) :: bun_Si2c_a  ! i2c states mapped to a
   type(cpl_bundle) :: bun_Fi2c_a  ! i2c fluxes mapped to a

   !--- fundamental maps ---
   type(cpl_map),target :: map_Sa2o  ! maps states a->o grids
   type(cpl_map),target :: map_Fa2o  ! maps fluxes a->o grids
   type(cpl_map),target :: map_So2a  ! maps states o->a grids 
   type(cpl_map),target :: map_Fo2a  ! maps fluxes o->a grids 
   type(cpl_map),target :: map_Xr2o  ! maps fluxes r->o grids 
   type(cpl_map),target :: map_ID    ! identity map

   !--- redundant maps ---
   type(cpl_map),pointer :: map_Fa2i  ! maps fluxes a->i grids
   type(cpl_map),pointer :: map_Fa2l  ! maps fluxes a->l grids
   type(cpl_map),pointer :: map_Sa2i  ! maps states a->i grids
   type(cpl_map),pointer :: map_Sa2l  ! maps states a->l grids

   type(cpl_map),pointer :: map_Fi2a  ! maps fluxes i->a grids
   type(cpl_map),pointer :: map_Fi2l  ! maps fluxes i->l grids
   type(cpl_map),pointer :: map_Fi2o  ! maps fluxes i->o grids
   type(cpl_map),pointer :: map_Si2a  ! maps states i->a grids
   type(cpl_map),pointer :: map_Si2l  ! maps states i->l grids
   type(cpl_map),pointer :: map_Si2o  ! maps states i->o grids

   type(cpl_map),pointer :: map_Fl2a  ! maps fluxes l->a grids
   type(cpl_map),pointer :: map_Fl2i  ! maps fluxes l->i grids
   type(cpl_map),pointer :: map_Fl2o  ! maps fluxes l->o grids
   type(cpl_map),pointer :: map_Sl2a  ! maps states l->a grids
   type(cpl_map),pointer :: map_Sl2i  ! maps states l->i grids
   type(cpl_map),pointer :: map_Sl2o  ! maps states l->o grids

   type(cpl_map),pointer :: map_Fo2i  ! maps fluxes o->i grids
   type(cpl_map),pointer :: map_Fo2l  ! maps fluxes o->l grids
   type(cpl_map),pointer :: map_So2i  ! maps states o->i grids
   type(cpl_map),pointer :: map_So2l  ! maps states o->l grids

   save

!EOP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: data_bundleInit - initialize all bundles
!
! !DESCRIPTION:
!     initialize all the bundle's needed
!
! !REVISION HISTORY:
!     2002-Mar-06 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine data_bundleInit()

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   ! input/output are all variables declared in this module

!EOP

   !--- formats ---
   character(len=*),parameter :: F00 = '("(data_bundleInit) ",60a)'

!-------------------------------------------------------------------------------
! METHOD:
!  1) assign the domain (an input) to the bundle
!  2) set size of bundle aVect to be the same size as that of the domain
!  3) use input str to define real aVect fields (must be a valid aVect list str)
!  4) hard-coded st the aVect has *no* integer data
! NOTE:
!  o memory is allocated for bun%data, but data values are left undefined
!-------------------------------------------------------------------------------

   write(6,F00) "initializing bundles"
   call shr_sys_flush(6)

   call cpl_bundle_init(bun_Xc2oSNAP_o,"Xc2oSNAP_o",cpl_fields_c2o_fields,dom_o)
   call cpl_bundle_init(bun_Xc2oPSUM_o,"Xc2oPSUM_o",cpl_fields_c2o_fields,dom_o)

   call cpl_bundle_init(bun_aoflux_o,"aoflux_o",&
      &"Faoc_sen&
      &:Faoc_lat&
      &:Faoc_lwup&
      &:Faoc_evap&
      &:Faoc_taux&
      &:Faoc_tauy&
      &:Faoc_tref&
      &:Faoc_qref&
      &:Faoc_duu10n&
      &:Faoc_swnet"   ,dom_o)
   call cpl_bundle_init(bun_aoflux_a,"aoflux_a",&
      &"Faoc_sen&
      &:Faoc_lat&
      &:Faoc_lwup&
      &:Faoc_evap&
      &:Faoc_taux&
      &:Faoc_tauy&
      &:Faoc_tref&
      &:Faoc_qref&
      &:Faoc_duu10n&
      &:Faoc_swnet"   ,dom_a)
   call cpl_bundle_init(bun_oalbedo_o,"oalbedo_o",&
      &"So_avsdr&
      &:So_anidr&
      &:So_avsdf&
      &:So_anidf"     ,dom_o)
   call cpl_bundle_init(bun_oalbedo_a,"oalbedo_a",&
      &"So_avsdr&
      &:So_anidr&
      &:So_avsdf&
      &:So_anidf"     ,dom_a)

   call cpl_bundle_init(bun_precip_o,"precip_o",&
      &"Faxc_rain&
      &:Faxc_snow"        ,dom_o)
   call cpl_bundle_init(bun_precip_a,"precip_a",&
      &"Faxc_rain&
      &:Faxc_snow"       ,dom_a)

   call cpl_bundle_init(bun_Sa2c_a,"Sa2c_a",cpl_fields_a2c_states,dom_a)
   call cpl_bundle_init(bun_Fa2c_a,"Fa2c_a",cpl_fields_a2c_fluxes,dom_a)
   call cpl_bundle_init(bun_Sa2c_o,"Sa2c_o",cpl_fields_a2c_states,dom_o)
   call cpl_bundle_init(bun_Fa2c_o,"Fa2c_o",cpl_fields_a2c_fluxes,dom_o)
   call cpl_bundle_init(bun_Sl2c_l,"Sl2c_l",cpl_fields_l2c_states,dom_l)
   call cpl_bundle_init(bun_Fl2c_l,"Fl2c_l",cpl_fields_l2c_fluxes,dom_l)
!  call cpl_bundle_init(bun_Sl2c_o,"Sl2c_o",cpl_fields_l2c_states,dom_o)
!  call cpl_bundle_init(bun_Fl2c_o,"Fl2c_o",cpl_fields_l2c_fluxes,dom_o)
   call cpl_bundle_init(bun_Xr2c_o,"Xr2c_o",cpl_fields_r2c_fields,dom_o)
   call cpl_bundle_init(bun_So2c_o,"So2c_o",cpl_fields_o2c_states,dom_o)
   call cpl_bundle_init(bun_Fo2c_o,"Fo2c_o",cpl_fields_o2c_fluxes,dom_o)
   call cpl_bundle_init(bun_So2c_a,"So2c_a",cpl_fields_o2c_states,dom_a)
   call cpl_bundle_init(bun_Fo2c_a,"Fo2c_a",cpl_fields_o2c_fluxes,dom_a)
   call cpl_bundle_init(bun_Si2c_i,"Si2c_i",cpl_fields_i2c_states,dom_i)
   call cpl_bundle_init(bun_Fi2c_i,"Fi2c_i",cpl_fields_i2c_fluxes,dom_i)
   call cpl_bundle_init(bun_Si2c_a,"Si2c_a",cpl_fields_i2c_states,dom_a)
   call cpl_bundle_init(bun_Fi2c_a,"Fi2c_a",cpl_fields_i2c_fluxes,dom_a)

end subroutine data_bundleInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: data_mapInit - initialize all mapping data 
!
! !DESCRIPTION:
!     initialize all the mapping data.
!
! !REVISION HISTORY:
!     2002-May-21 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine data_mapInit()

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   ! input could be map files names, but these are currently hard-coded
   ! output are the mapping data types declared in this module

!EOP

   !--- locatl ---
   character(256)    :: mapName      ! ID string given to a map
   character(256)    :: fileName     ! mapping data input file name


   !--- formats ---
   character(len=*),parameter :: F00 = '("(data_mapInit) ",60a)'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "initializing mapping data..."
   call shr_sys_flush(6)

   !------------------------------------------------------------
   ! initialize identity map
   !------------------------------------------------------------

   map_ID%IDtype = 1

   !------------------------------------------------------------
   ! initialize maps between atm & ocn 
   !------------------------------------------------------------

   fileName = cpl_control_mapFn_a2oS
   mapName  = "map_Sa2o"
   call cpl_map_init(map_Sa2o,dom_a,dom_o,mapName,fileName,"src")

   fileName = cpl_control_mapFn_a2oF
   mapName  = "map_Fa2o"
   call cpl_map_init(map_Fa2o,dom_a,dom_o,mapName,fileName,"src")

   fileName = cpl_control_mapFn_o2aF
   mapName  = "map_Fo2a"
   call cpl_map_init(map_Fo2a,dom_o,dom_a,mapName,fileName,"dst")

   fileName = cpl_control_mapFn_o2aS
   mapName  = "map_So2a"
   call cpl_map_init(map_So2a,dom_o,dom_a,mapName,fileName,"dst")

   fileName = cpl_control_mapFn_r2o
   mapName  = "map_Xr2o"
   call cpl_map_init(map_Xr2o,dom_r,dom_o,mapName,fileName,"dst")

   !------------------------------------------------------------
   ! initialize rest of redundant maps
   !------------------------------------------------------------
   map_Sa2i => map_Sa2o
   map_Fa2i => map_Fa2o
   map_Si2a => map_So2a
   map_Fi2a => map_Fo2a

   write(6,F00) "... done initializing mapping."
   call shr_sys_flush(6)

end subroutine data_mapInit

!===============================================================================

END MODULE data_mod

