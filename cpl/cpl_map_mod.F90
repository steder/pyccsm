!===============================================================================
! CVS: $Id: cpl_map_mod.F90,v 1.1.1.1 2005/02/03 22:29:01 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/cpl/cpl_map_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_map_mod -- mapping subsystem module
!
! !DESCRIPTION:
!    This is the module represents a major subsystem of cpl6. "Mapping" refers 
!    to the transfer of 2d field data from one domain/grid to another.  It is 
!    often desirable that maps have the properties of being {\it smooth} and 
!    {\it conservative}.  Common mapping techniques are bilinear interplation
!    and area-averaging.  Mapping is implemented by a sparse matrix multiply.
!    This module defines the sparse matrix data type used for mapping and 
!    handles the actual mapping of data (matrix multiply).  This module also 
!    handles the initialization and error checking of sparse matrix data.
!
! !REVISION HISTORY:
!    2001-Aug-14 - B. Kauffman -- gathered all mapping routines into this module
!    2001-May-20 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

module cpl_map_mod

! !USES:

   use cpl_mct_mod        ! mct interface
   use cpl_domain_mod     ! data type & methods
   use cpl_bundle_mod     ! data type & methods
   use cpl_comm_mod       ! global data
   use cpl_kind_mod       ! kinds
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use cpl_control_mod, only: bfbflag=>cpl_control_bfbflag
   use shr_sys_mod        ! flush
   use shr_mpi_mod        ! mpi layer

   implicit none

   private   ! except

! !PUBLIC TYPES:

   public :: cpl_map
   
   type cpl_map
     character(CL)            :: name    ! text ID of mapping data
     type(cpl_mct_sMat)       :: sMat    ! the mct sparse matrix data type
     type(cpl_domain),pointer :: src     ! the associated source domain
     type(cpl_domain),pointer :: dst     ! the associated destination domain
     type(cpl_domain)         :: new     ! new/intermediate domain required by mct
     type(cpl_mct_rearr)      :: rearr   ! rearranger to/from new
     character(3)             :: newtype ! intermediate domain type: src or dst ?
     integer(IN)              :: IDtype  ! 0=normal, 1=identity(ID)
     type(cpl_mct_Avect)      :: areasrc ! area of src grid from mapping file
     type(cpl_mct_Avect)      :: areadst ! area of dst grid from mapping file
   end type cpl_map

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_map_init      ! initialize a map
   public :: cpl_map_clean     ! clean/dealloc a map
   public :: cpl_map_info      ! obtain information about a map
   public :: cpl_map_bun       ! map from one bundle to another
   public :: cpl_map_npFix     ! fix NP values wrt mapping vector fields

   interface cpl_map_npFix; module procedure cpl_map_npFixNew3; end interface
!  interface cpl_map_npFix; module procedure cpl_map_npFixNew2; end interface
!  interface cpl_map_npFix; module procedure cpl_map_npFixNew; end interface
!  interface cpl_map_npFix; module procedure cpl_map_npFixOld; end interface
!  interface cpl_map_npFix; module procedure cpl_map_npFixNone; end interface 


! !PUBLIC DATA MEMBERS:

   character(*),parameter,public :: cpl_map_areaAV_field = 'aream'

!EOP

   !--- module variables ---
   character(*),parameter :: modName = "cpl_map_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_init - create a map between two domains
!
! !DESCRIPTION:
!    Given two domains, map\_init initializes/creates a map between them.
! 
! !REMARKS:
!    MCT currently does not support the creation of a map between two arbitrary 
!    domains, rather one of the two domains must have a decomposition of MCT's 
!    choosing, thus the calling routine can specify which of the two domains 
!    should be altered and this new/required domain can be returned.
!
! !REVISION HISTORY:
!    2001-Jun-14 - T. Craig - first functioning version
!    2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_map_init(map_X,dom_src,dom_dst,mapName,fileName,newdom,adj_areas)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map)   ,intent(out)        :: map_X    ! map_X data
   type(cpl_domain),intent( in),target :: dom_src  ! map's source domain
   type(cpl_domain),intent( in),target :: dom_dst  ! map's destination domain
   character(*),intent( in)            :: mapName  ! map's ID string 
   character(*),intent( in)            :: fileName ! file containing map data
   character(*),intent( in)            :: newdom   ! which domain to alter (for mct)
   logical,optional                    :: adj_areas! flag to adjust areas

!EOP

   !--- local ---
   type(cpl_map) :: map_0            ! temporary sMat on root only
   type(cpl_mct_rearr) :: lrearr     ! local rearranger to/from new
   integer(IN)   :: rCode            ! rCode flag
   integer(IN)   :: lSize            ! size of Attribute Vector
   integer(IN)   :: ilrow            ! index of row field in Smat
   integer(IN)   :: ilcol            ! index of col field in Smat
   integer(IN)   :: iwgt             ! index of wgt field in Smat
   integer(IN)   :: n                ! generic index
   integer(IN)   :: nflds            ! field number of AV
   integer(IN)   :: nfldd            ! field number of AV
   integer(IN)   :: nas              ! field number of area in src AV
   integer(IN)   :: nad              ! field number of area in dst AV
   integer(IN)   :: msks             ! field number of mask in src AV
   integer(IN)   :: mskd             ! field number of mask in dst AV
   integer(IN)   :: isrc             ! location of src index for wgt
   integer(IN)   :: idst             ! location of dst index for wgt

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_init) '
   character(*),parameter :: F00 = "('(cpl_map_init) ',8a)"
   character(*),parameter :: F01 = "('(cpl_map_init) ',2(a,i8))"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

    write(6,F00) "initialize map: ",trim(mapName)

   !----------------------------------------------------------------------------
   ! set map's src & dest domains, name, "new mct domain" choice 
   !----------------------------------------------------------------------------
   map_0%name    =  trim(mapName) // ", source file = " // trim(fileName)
   map_0%newtype =  newdom
   map_0%src     => dom_src
   map_0%dst     => dom_dst
   map_0%IDtype  =  0

   if (bfbflag) then
     map_0%newtype = 'src'
     write(6,*) subName,': bfbflag = ',bfbflag, &
       ': overwriting newtype from ',newdom,' to ',map_0%newtype
   endif

   !----------------------------------------------------------------------------
   ! read & test the map data on root processor only 
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid == 0) then
      call cpl_map_read(map_0,fileName)
      if (dbug > 1) then
         call cpl_map_test(map_0)
      else
         write(6,F01) "skipping map test, dbug level = ",dbug
      end if
   endif
   call shr_mpi_barrier(cpl_comm_comp,subName//" barrier")

   !----------------------------------------------------------------------------
   ! scatter map_X data and create new/intermediate domain as required by MCT
   !----------------------------------------------------------------------------
   map_X%name    =  map_0%name    
   map_X%newtype =  map_0%newtype 
   map_X%src     => dom_src
   map_X%dst     => dom_dst
   map_X%IDtype  =  map_0%IDtype

   call cpl_mct_aVect_clean(map_X%areasrc)
   call cpl_mct_aVect_scatter(map_0%areasrc,map_X%areasrc, map_X%src%gsMap,0,cpl_comm_comp,rCode)
   call cpl_mct_aVect_clean(map_X%areadst)
   call cpl_mct_aVect_scatter(map_0%areadst,map_X%areadst, map_X%dst%gsMap,0,cpl_comm_comp,rCode)

   if (dbug > 2) then
      if (cpl_comm_comp_pid == 0) then
         write(6,*) subName,'lSize of src & dest',&
                    cpl_mct_aVect_lSize(map_0%areasrc), &
                    cpl_mct_aVect_lSize(map_0%areadst)
         write(6,*) subName,'min/max src ',minval(map_0%areasrc%rAttr(1,:)), &
                                           maxval(map_0%areasrc%rAttr(1,:))
         write(6,*) subName,'min/max dst ',minval(map_0%areadst%rAttr(1,:)), &
                                           maxval(map_0%areadst%rAttr(1,:))
      endif
      call cpl_mct_aVect_info(4,map_X%areasrc,cpl_comm_comp,cpl_comm_comp_pid,istr='map areasrc')
      call cpl_mct_aVect_info(4,map_X%areadst,cpl_comm_comp,cpl_comm_comp_pid,istr='map areadst')
   endif

   if (map_0%newtype == "dst") then !--- create new/intermediate destination domain ---
      write(6,F00) 'scatter matrix by column...'  ; call shr_sys_flush(6)

      call cpl_mct_sMat_scatterByCol(dom_src%gsMap, map_0%sMat, map_X%sMat,0, cpl_comm_comp, rCode)

      if (dbug > 1 .or. rCode /= 0) then
         write(6,F01) 'scatter   sMat return code =',rCode 
         write(6,F01) 'scattered sMat rows x cols =',map_X%sMat%nrows,' x',map_X%sMat%ncols
         write(6,F01) 'scattered sMat lSize  ',cpl_mct_sMat_lSize(map_X%sMat)
         write(6,F01) 'scattered sMat GnumEl ',cpl_mct_sMat_GNumEl(map_X%sMat,cpl_comm_comp)
      end if

      map_X%new%name  =  'dst intermediate map'
      map_X%new%suffix=  dom_dst%suffix
      map_X%new%n     =  dom_dst%n
      map_X%new%ni    =  dom_dst%ni
      map_X%new%nj    =  dom_dst%nj

      call cpl_mct_sMat_2YgsMap(map_X%sMat,map_X%new%gsMap,0,cpl_comm_comp,cpl_comm_mph_cid)

      if (dbug > 1) then
         write(6,F01)'new gsMap gSize =', cpl_mct_gsMap_gSize(map_X%new%gsMap)
         write(6,F01)'new gsMap lSize =', cpl_mct_gsMap_lSize(map_X%new%gsMap,cpl_comm_comp)
         call shr_sys_flush(6)
      end if

      call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%new%gsMap,"ROW"   ,cpl_comm_comp)
      call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%src%gsMap,"COLUMN",cpl_comm_comp)
      call cpl_mct_rearr_init (map_X%new%gsMap,map_X%dst%gsMap,cpl_comm_comp,map_X%rearr)
      call cpl_mct_rearr_init (map_X%dst%gsMap,map_X%new%gsMap,cpl_comm_comp,lrearr)
      call cpl_mct_aVect_init (map_X%new%lGrid,map_X%dst%lGrid,cpl_mct_gsMap_lSize(map_X%new%gsMap,cpl_comm_comp))
      call cpl_mct_rearr_rearrange(map_X%dst%lGrid,map_X%new%lGrid,lrearr)
      call cpl_mct_rearr_clean(lrearr)

   elseif (map_0%newtype == "src") then !--- create new/intermediate source domain ---
      write(6,F00) 'scatter matrix by row...' ; call shr_sys_flush(6)

      call cpl_mct_sMat_scatterByRow(dom_dst%gsMap,map_0%sMat,map_X%sMat,0, cpl_comm_comp,rCode)
      write(6,F01) 'SM_ScatterByRow, rCode=',rCode ; call shr_sys_flush(6)

      if (dbug > 1 .or. rCode /= 0) then
         write(6,F01) 'scatter   sMat return code =',rCode 
         write(6,F01) 'scattered sMat rows x cols =',map_X%sMat%nrows,' x',map_X%sMat%ncols
         write(6,F01) 'scattered sMat lSize  ',cpl_mct_sMat_lSize(map_X%sMat)
         write(6,F01) 'scattered sMat gNumEl ',cpl_mct_sMat_GNumEl(map_X%sMat,cpl_comm_comp)
         call shr_sys_flush(6)
      end if

      map_X%new%name  =  'src intermediate map'
      map_X%new%suffix=  dom_src%suffix
      map_X%new%n     =  dom_src%n
      map_X%new%ni    =  dom_src%ni
      map_X%new%nj    =  dom_src%nj

      call cpl_mct_sMat_2XgsMap(map_X%sMat,map_X%new%gsMap,0,cpl_comm_comp,cpl_comm_mph_cid)

      if (dbug > 1) then
         write(6,F01)'new gsMap gSize =', cpl_mct_gsMap_gSize(map_X%new%gsMap)
         write(6,F01)'new gsMap lSize =', cpl_mct_gsMap_lSize(map_X%new%gsMap,cpl_comm_comp)
         call shr_sys_flush(6)
      end if

      call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%dst%gsMap,"ROW"   ,cpl_comm_comp)
      call cpl_mct_sMat_g2lMat(map_X%sMat,map_X%new%gsMap,"COLUMN",cpl_comm_comp)
      call cpl_mct_rearr_init (map_X%src%gsMap,map_X%new%gsMap,cpl_comm_comp,map_X%rearr)
      call cpl_mct_aVect_init (map_X%new%lGrid,map_X%dst%lGrid,cpl_mct_gsMap_lSize(map_X%new%gsMap,cpl_comm_comp))
      call cpl_mct_rearr_rearrange(map_X%src%lGrid,map_X%new%lGrid,map_X%rearr)

   else !--- no other valid choices ---
      write(6,F00) 'ERROR: invalid newdom value = ',map_0%newtype
      call shr_sys_abort(trim(subName)//" invalid newdom value")
   endif

#ifdef CPP_VECTOR
   ! initialize the vector parts of the sMat
   call cpl_mct_sMat_Vecinit(map_X%sMat)
#endif

   !----------------------------------------------------------------------------
   ! adjust mapping weights based on areas of domains vs areas in mappings
   !----------------------------------------------------------------------------
   if (present(adj_areas)) then
   if (adj_areas) then

     write(6,F00) 'ERROR: do not use adj_areas option right now'
     call shr_sys_abort(trim(subName)//" adj_areas used")

     !--- adjust mapping weights based on areas of domains vs areas in mappings
     !--- S * src_area/map_area -> S' -> S2D -> D' *map_area/dst_area -> D
     !--- S  = fluxes associated with src grid and src areas
     !--- S' = fluxes associated with src grid and mapping areas
     !--- D  = fluxes associated with dst grid and dst areas
     !--- D' = fluxes associated with dst grid and mapping areas
     ilrow = cpl_mct_sMat_indexIA(map_X%sMat,'lrow')
     ilcol = cpl_mct_sMat_indexIA(map_X%sMat,'lcol')
     iwgt  = cpl_mct_sMat_indexRA(map_X%sMat,'weight')
     lSize = cpl_mct_sMat_lSize(map_X%sMat)
     nflds = cpl_mct_aVect_indexRA(map_X%areasrc,cpl_map_areaAV_field)
     nfldd = cpl_mct_aVect_indexRA(map_X%areadst,cpl_map_areaAV_field)
     nas   = cpl_mct_aVect_indexRA(map_X%src%lGrid,"area",perrWith=subName)
     nad   = cpl_mct_aVect_indexRA(map_X%dst%lGrid,"area",perrWith=subName)
     msks  = cpl_mct_aVect_indexRA(map_X%src%lGrid,"mask",perrWith=subName)
     mskd  = cpl_mct_aVect_indexRA(map_X%dst%lGrid,"mask",perrWith=subName)

     write(6,*) trim(subName),' adjusting areas for ',trim(mapName),' ', &
        ilrow,ilcol,iwgt,lSize,nflds,nfldd,nas,nad,msks,mskd

     do n=1,lSize
       isrc = map_X%sMat%data%iAttr(ilcol,n)
       idst = map_X%sMat%data%iAttr(ilrow,n)
       if (abs(map_X%src%lGrid%rAttr(msks,isrc)) < 1.0e-06) &
          write(6,*) trim(subName),n,isrc,map_X%src%lGrid%rAttr(msks,isrc)
       if (abs(map_X%dst%lGrid%rAttr(mskd,idst)) < 1.0e-06) &
          write(6,*) trim(subName),n,idst,map_X%dst%lGrid%rAttr(mskd,idst)
       map_X%sMat%data%rAttr(iwgt,n) = map_X%sMat%data%rAttr(iwgt,n) * &
          map_X%areadst%rAttr(nfldd,idst) /   &
          map_X%dst%lGrid%rAttr(nas,idst) *   &
          map_X%src%lGrid%rAttr(nad,isrc) /   &
          map_X%areasrc%rAttr(nflds,isrc) 
     enddo
   endif
   endif

   map_0%name    = '<clean>'
   map_0%newtype = '<clean>'
   nullify(map_0%src)
   nullify(map_0%dst)
   if (cpl_comm_comp_pid == 0) then
     call cpl_mct_sMat_clean(map_0%sMat)
     call cpl_mct_aVect_clean(map_0%areasrc)
     call cpl_mct_aVect_clean(map_0%areadst)
   endif
   map_0%IDtype = 0

   write(6,F00) "done initializing map: ",trim(mapName) ; call shr_sys_flush(6)
   call shr_sys_flush(6)

end subroutine cpl_map_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_clean - deallocate a map data type
!
! !DESCRIPTION:
!    Clean a map
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-Jan-20 - T. Craig - first functioning version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_map_clean(mapping)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map),intent(inout)    :: mapping ! mapping data

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_clean) '
   character(*),parameter :: F00 = "('(cpl_map_clean) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (dbug>2) write(6,F00) 'cleaning map, name=',trim(mapping%name)

   mapping%name    = '<clean>'
   mapping%newtype = '<clean>'
   nullify(mapping%src)
   nullify(mapping%dst)
   call cpl_mct_sMat_clean(mapping%sMat)
   call cpl_mct_rearr_clean(mapping%rearr)
   call cpl_domain_clean(mapping%new)
   call cpl_mct_aVect_clean(mapping%areasrc)
   call cpl_mct_aVect_clean(mapping%areadst)
   mapping%IDtype = 0

end subroutine cpl_map_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_info - print some info about the map
!
! !DESCRIPTION:
!    print some info about the map datatype
! 
! !REMARKS:
!
! !REVISION HISTORY:
!    2002-Jan-14 - T. Craig - first functioning version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_map_info(mapping)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map)   ,intent(in)                  :: mapping ! mapping data

!EOP

   !--- local ---
   integer(IN)         :: i,j,k,n      ! generic indicies
   real(R8)            :: mn,mx        ! min,max

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_info) '
   character(*),parameter :: F00 = '(a,a,2f12.3)'

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   write(6,*) ' '
   write(6,*) subName,' mapping%name:',trim(mapping%name)
   write(6,*) subName,' mapping%newtype:',trim(mapping%newtype)
   write(6,*) subName,' mapping%IDtype:',mapping%IDtype
   write(6,*) subName,' mapping%src: '
   call cpl_domain_info(mapping%src)
   write(6,*) subName,' mapping%dst: '
   call cpl_domain_info(mapping%dst)
   write(6,*) subName,' newtype: ',mapping%newtype
   write(6,*) subName,' mapping%new: '
   call cpl_domain_info(mapping%new)

   write(6,*) subName,' mapping%sMat: '
   write(6,*) subName,' mapping%sMat%nrows: ',mapping%sMat%nrows
   write(6,*) subName,' mapping%sMat%ncols: ',mapping%sMat%ncols

   call cpl_mct_aVect_info(2,mapping%sMat%data,cpl_comm_comp,cpl_comm_comp_pid)
   write(6,*) ' '
   call shr_sys_flush(6)

end subroutine cpl_map_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_bun - map a bundle from one domain to domain
!
! !DESCRIPTION:
!    Map a bundle of fields from one domain to another.
!    It is assumed that both bundles have an identical list of fields.
!
! !REVISION HISTORY:
!    20May01 - T. Craig -- first prototype
!    15Jun01 - E.T. Ong -- Removed zeroing of bunn%data and buno%data -
!                          this is done in mct calls.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_bun(buni,buno,mapx,bunfs,fsname,bunfd,fdname,mvector)

! !USES:

   use shr_timer_mod       ! share timer routines

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout)       :: buni    ! input bundle
   type(cpl_bundle),intent(out)         :: buno    ! output bundle
   type(cpl_map)   ,intent(inout)       :: mapx    ! mapping between two domains
   type(cpl_bundle),intent(in),optional :: bunfs   ! src fraction input bundle
   character(*)    ,intent(in),optional :: fsname  ! name of field in bunfs
   type(cpl_bundle),intent(in),optional :: bunfd   ! dst fraction input bundle
   character(*)    ,intent(in),optional :: fdname  ! name of field in bunfd
   logical         ,intent(in),optional :: mvector  ! enable vector-friendly mapping

!EOP

   !--- local ---
   type(cpl_bundle):: buni_local  ! local copy of buni for optional arguments
   type(cpl_bundle):: bunn        ! temp bundle with decomp that mct chooses
   logical         :: Sum         ! Should rearranger do a sum?
   character(CL)   :: name        ! name string for temp bundle
   logical         :: normalize   ! true if optional arguments present
   logical         :: usevector   ! true if we want to use the vector-friendly code
   integer(IN)     :: n,m         ! generic indicies
   integer(IN)     :: npts        ! number of points (local) in an aVect field
   integer(IN)     :: nfld        ! number of fields (local) in an aVect field
   integer(IN)     :: kfld        ! field number of fld1 in bun1
   real(R8)        :: fdr         ! 1/bunfd

   logical,save :: first_call = .true.     ! first time in subroutine
   integer(IN),save :: tmap1,tmap2,tmap3,tmap4,tmap5,tmap6
   integer(IN),save :: tmap7,tmap8,tmap9,tmap10

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_bun) '
   character(*),parameter :: F00 = '("(cpl_map_bun) ",8a)'
   character(*),parameter :: F01 = '("(cpl_map_bun) ",3a,i5)'

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     first_call = .false.
     call shr_timer_get(tmap1,'tmap1')
     call shr_timer_get(tmap2,'tmap2')
     call shr_timer_get(tmap3,'tmap3')
     call shr_timer_get(tmap4,'tmap4')
     call shr_timer_get(tmap5,'tmap5')
     call shr_timer_get(tmap6,'tmap6')
     call shr_timer_get(tmap7,'tmap7')
     call shr_timer_get(tmap8,'tmap8')
     call shr_timer_get(tmap9,'tmap9')
     call shr_timer_get(tmap10,'tmap10')
   endif

   !----------------------------------------------------------------------------
   ! By default, turn on vector parts if CPP_VECTOR is defined
   ! But behaviour can be changed by input argument.
   !----------------------------------------------------------------------------
#ifdef CPP_VECTOR
   usevector = .true.
#else
   usevector = .false.
#endif
   ! override the defaults with the input argument
   if(present(mvector)) then
     usevector = mvector
   endif

   !----------------------------------------------------------------------------
   ! identity map => "mapping is a simply copy the bundle
   !----------------------------------------------------------------------------
   if (mapx%IDtype == 1) then
     call cpl_bundle_copy(buni,outbun=buno)
     return
   endif

   !----------------------------------------------------------------------------
   ! if normalization is requested, create a local copy of bunlde
   !----------------------------------------------------------------------------
   call shr_timer_start(tmap1)
   normalize = .false.
   if (present(bunfs) .and.present(bunfd).and.    &
       present(fsname).and.present(fdname)) then
     normalize = .true.
   endif
   if ((present(bunfs).and..not.present(bunfd)) .or.  &
       (present(bunfd).and..not.present(bunfs)) .or.  &
       (present(bunfs).and..not.present(fsname)) .or. &
       (present(bunfd).and..not.present(fdname))) then
     write(6,*) subName,' ERROR: optional arguments inconsistent '
     call shr_sys_abort(subName)
   endif
   call shr_timer_stop(tmap1)

   call shr_timer_start(tmap2)
   if (normalize) then
     name = trim(buni%name) // "_local"
     call cpl_bundle_initv(buni_local,name,buni,buni%dom)
!     call cpl_bundle_copy(buni,outbun=buni_local)
!     call cpl_bundle_mult(buni_local,bunfs,fsname)
      call cpl_bundle_mult(buni_local,bunfs,fsname,initbun=buni)
   endif
   call shr_timer_stop(tmap2)

   !----------------------------------------------------------------------------
   ! do the mapping 
   !----------------------------------------------------------------------------
   if (mapx%newtype == "src") then !--- use intermediate src domain ---

     !--- initialize temporary bundle ---
     name = trim(buno%name) // "_map" // "src "
     call shr_timer_start(tmap3)
     if (normalize) then
       call cpl_bundle_initv(bunn,name,buni_local,mapx%new)
     else
       call cpl_bundle_initv(bunn,name,buni,      mapx%new)
     endif
     call shr_timer_stop(tmap3)

     !--- redistribute the data ---
     call shr_timer_start(tmap4)
     if (normalize) then
       call cpl_mct_rearr_rearrange(buni_local%data,bunn%data,mapx%rearr,vector=usevector)
     else
       call cpl_mct_rearr_rearrange(buni%data      ,bunn%data,mapx%rearr,vector=usevector)
     endif
     call shr_timer_stop(tmap4)

     !--- map the data ---
     call shr_timer_start(tmap5)
     call cpl_mct_Smat_AvMult(bunn%data,mapx%sMat,buno%data,vector=usevector)
     call shr_timer_stop(tmap5)

   else if (mapx%newtype == "dst") then ! use intermediate dest domain ---

     !--- initialize temporary bundle ---
     name = trim(buno%name) // "_map" // "dst "
     call shr_timer_start(tmap6)
     call cpl_bundle_initv(bunn,name,buno,mapx%new)
     call shr_timer_stop(tmap6)

     !--- map the data ---
     call shr_timer_start(tmap7)
     if (normalize) then
       call cpl_mct_Smat_AvMult(buni_local%data,mapx%sMat,bunn%data,vector=usevector)
     else
       call cpl_mct_Smat_AvMult(buni%data      ,mapx%sMat,bunn%data,vector=usevector)
     endif
     call shr_timer_stop(tmap7)

     !--- redistribute the data ---
     call shr_timer_start(tmap8)
     Sum=.TRUE.  ! this may not be necessary -RLJ
     call cpl_mct_rearr_rearrange(bunn%data,buno%data,mapx%rearr,Sum=Sum,vector=usevector)
     call shr_timer_stop(tmap8)

   else
      write(6,*) subName,' error invalid newtype ',mapx%newtype
   endif

   !----------------------------------------------------------------------------
   ! finish normalization
   !----------------------------------------------------------------------------
   call shr_timer_start(tmap9)
   if (normalize) then
     call cpl_bundle_divide(buno,bunfd,fdname)
     call cpl_bundle_clean(buni_local)
   endif
   call shr_timer_stop(tmap9)

   !----------------------------------------------------------------------------
   ! clean up the temporary bundle
   !----------------------------------------------------------------------------
   call shr_timer_start(tmap10)
   if (buni%cnt /= 1) write(6,F01) &
       "WARNING:  bundle ",trim(buni%name)," has accum count =",buni%cnt
   buno%cnt = buni%cnt
   call cpl_bundle_clean(bunn)
   call shr_timer_stop(tmap10)

end subroutine cpl_map_bun

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE:  cpl_map_read - Read in a SparseMatrix
!
! !DESCRIPTION: 
!     Read in mapping matrix data from a SCRIP netCDF data file.
!
! !REMARKS:
!   This routine must be generalized to
!   o deal with mapping data between arbitrary domains
!   o read in data files with arbitrary file names
!   o read in netCDF data files
!
! !SEE ALSO: 
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY: 
!    2002-Mar-19 - B. Kauffman -- first version
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_read(mapping,fileName)

! !USES:

#include <netcdf.inc>

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map)   ,intent(out) :: mapping   ! mapping data
   character(*),intent(in)  :: filename  ! netCDF file to read

!EOP

   !--- local ---
   integer(IN)           :: n       ! generic loop indicies
   integer(IN)           :: na      ! size of source domain
   integer(IN)           :: nb      ! size of destination domain
   integer(IN)           :: ns      ! number of non-zero elements in matrix
   integer(IN)           :: ni,nj   ! number of row and col in the matrix
   integer(IN)           :: igrow   ! aVect index for matrix row
   integer(IN)           :: igcol   ! aVect index for matrix column
   integer(IN)           :: iwgt    ! aVect index for matrix element
   integer(IN)           :: iarea   ! aVect index for area

   real(R8)   ,allocatable :: S  (:)  ! matrix elements
   integer(IN),allocatable :: row(:)  ! correpsonding row    index
   integer(IN),allocatable :: col(:)  ! correpsonding column index
   real(R8)   ,allocatable :: area_a(:)  ! area of a grid in map
   real(R8)   ,allocatable :: area_b(:)  ! area of b grid in map

   character,allocatable :: str(:)  ! variable length char string
   character(CL)         :: attstr  ! netCDF attribute name string
   integer(IN)           :: rcode   ! netCDF routine return code
   integer(IN)           :: fid     ! netCDF file      ID
   integer(IN)           :: vid     ! netCDF variable  ID
   integer(IN)           :: did     ! netCDF dimension ID

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_read) '
   character(*),parameter :: F00 = '("(cpl_map_read) ",4a)'
   character(*),parameter :: F01 = '("(cpl_map_read) ",2(a,i7))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "reading mapping matrix data..."

   !----------------------------------------------------------------------------
   ! open & read the file
   !----------------------------------------------------------------------------
   write(6,F00) "* file name                  : ",trim(fileName)
   rcode = nf_open(filename,NF_NOWRITE,fid)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   !--- allocate memory & get matrix data ----------
   rcode = nf_inq_dimid (fid, 'n_s', did)  ! size of sparse matrix
   rcode = nf_inq_dimlen(fid, did  , ns)
   rcode = nf_inq_dimid (fid, 'n_a', did)  ! size of  input vector
   rcode = nf_inq_dimlen(fid, did  , na)
   rcode = nf_inq_dimid (fid, 'n_b', did)  ! size of output vector
   rcode = nf_inq_dimlen(fid, did  , nb)

   write(6,F01) "* matrix dimensions rows x cols :",na,' x',nb
   write(6,F01) "* number of non-zero elements: ",ns

   allocate(row(ns),col(ns),S(ns),stat=rcode)
   allocate(area_a(ns),area_b(nb),stat=rcode)

   rcode = nf_inq_varid     (fid,'S'  ,vid)
   rcode = nf_get_var_double(fid,vid  ,S  )
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
   rcode = nf_inq_varid     (fid,'row',vid)
   rcode = nf_get_var_int   (fid,vid  ,row)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
   rcode = nf_inq_varid     (fid,'col',vid)
   rcode = nf_get_var_int   (fid,vid  ,col)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   rcode = nf_inq_varid     (fid,'area_a',vid)
   rcode = nf_get_var_double(fid,vid     ,area_a)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
   rcode = nf_inq_varid     (fid,'area_b',vid)
   rcode = nf_get_var_double(fid,vid     ,area_b)
   if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)

   rcode = nf_close(fid)

   !----------------------------------------------------------------------------
   ! init & load the mct sMat data type
   !----------------------------------------------------------------------------
   ! cpl_mct_sMat_init must be given the number of rows and columns that
   ! would be in the full matrix.  Nrows= size of output vector=nb.
   ! Ncols = size of input vector = na.
   call cpl_mct_sMat_init(mapping%sMat, nb, na, ns)

   igrow = cpl_mct_sMat_indexIA(mapping%sMat,'grow')
   igcol = cpl_mct_sMat_indexIA(mapping%sMat,'gcol')
   iwgt  = cpl_mct_sMat_indexRA(mapping%sMat,'weight')

   mapping%sMat%data%iAttr(igrow,:) = row(:)
   mapping%sMat%data%iAttr(igcol,:) = col(:)
   mapping%sMat%data%rAttr(iwgt ,:) =   S(:)

   call cpl_mct_aVect_init(mapping%areasrc,' ',cpl_map_areaAV_field,na)
   call cpl_mct_aVect_init(mapping%areadst,' ',cpl_map_areaAV_field,nb)

   iarea = cpl_mct_aVect_indexRA(mapping%areasrc,cpl_map_areaAV_field)
   mapping%areasrc%rAttr(iarea,1:na) = area_a(1:na)
   iarea = cpl_mct_aVect_indexRA(mapping%areadst,cpl_map_areaAV_field)
   mapping%areadst%rAttr(iarea,1:nb) = area_b(1:nb)

   deallocate(row, col, S, stat=rcode)
   if (rcode /= 0) call cpl_mct_perr_die(subName,':: deallocate(row...',rcode)

   deallocate(area_a, area_b, stat=rcode)
   if (rcode /= 0) call cpl_mct_perr_die(subName,':: deallocate(area_a...',rcode)

   write(6,F00) "... done reading file"
   call shr_sys_flush(6)

end subroutine cpl_map_read

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE:  cpl_map_test - consistancy test for a uniprocessor SparseMatrix
!
! !DESCRIPTION: 
!     Perform a variety of tests on a *uniprocessor* mct sparse-matrix data type
!
! !SEE ALSO: 
!    mct/m_SparseMatrix.F90 (MCT source code)
!
! !REVISION HISTORY: 
!    2002-Mar-19 - B. Kauffman -- first version
! 
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_test(map_X)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_map),intent(in) :: map_X  ! mapping data

!EOP

   !--- local ---
   integer(IN)    :: ns           ! size of aVect
   integer(IN)    :: igrow        ! aVect index for matrix row
   integer(IN)    :: igcol        ! aVect index for matrix column
   integer(IN)    :: iwgt         ! aVect index for matrix element

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_test) '
   character(*),parameter :: F00 = '("(cpl_map_test) ",4a)'
   character(*),parameter :: F01 = '("(cpl_map_test) ",2(a,i9))'
   character(*),parameter :: F02 = '("(cpl_map_test) ",2(a,es11.3))'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "consistancy check on mapping matrix data..."

   write(6,F00) "* map name              : ",trim(map_X%name)
   write(6,F00) "* map src domain name   : ",trim(map_X%src%name)
   write(6,F00) "* map dst domain name   : ",trim(map_X%dst%name)
   write(6,F00) "* map newtype selection : ",trim(map_X%newtype)
   write(6,F01) "* map IDtype            : ",map_X%IDtype

   ns    = cpl_mct_sMat_lSize  (map_X%sMat)
   igRow = cpl_mct_sMat_indexIA(map_X%sMat,'grow')
   igCol = cpl_mct_sMat_indexIA(map_X%sMat,'gcol')
   iwgt  = cpl_mct_sMat_indexRA(map_X%sMat,'weight')

   write(6,F01) '* number of rows x cols = ',map_X%sMat%nRows, &
                                       ' x ',map_X%sMat%nCols

   write(6,F01) "* number non-zero elements =",ns

   write(6,F02) "* min/max   element value = ",                  &
                   minval(map_X%sMat%data%rAttr(iwgt ,:)),", ", &
                   maxval(map_X%sMat%data%rAttr(iwgt ,:))

   write(6,F01) "* first/last row    value = ",          &
                   map_X%sMat%data%iAttr(igRow, 1),", ", &
                   map_X%sMat%data%iAttr(igRow,ns)
   write(6,F01) "* min/max    row    value = ",                  &
                   minval(map_X%sMat%data%iAttr(igRow,:)),", ", &
                   maxval(map_X%sMat%data%iAttr(igRow,:))

   write(6,F01) "* first/last column value = ",          &
                   map_X%sMat%data%iAttr(igCol, 1),", ", &
                   map_X%sMat%data%iAttr(igCol,ns)
   write(6,F01) "* min/max    column value = ",                  &
                   minval(map_X%sMat%data%iAttr(igCol,:)),", ", &
                   maxval(map_X%sMat%data%iAttr(igCol,:))

   call shr_sys_flush(6)

!  NOTES: 
!  o stuff below is removed because it alters the mapping data.
!  o any testing done should alter the data being tested.
!  o if any sorting is required, do it on file input
!
!  type(cpl_mct_list) :: sortKey       ! SparseMatrix sorting key list:
!  logical        :: descending(2) ! Descending order flags for sort test 2a
!
!  !----------------------------------------------------------------------------
!  ! sort and re-test
!  !----------------------------------------------------------------------------
!
!  write(6,F00) "Descending sort test..."
!
!  call cpl_mct_list_init(sortKey,"grow:gcol")
!  descending = .true.
!  call cpl_mct_sMat_SortPermute(sMat, sortKey, descending)
!
!  write(6,*) subName, "sMat%data%iAttr(igrow, 1) = ",sMat%data%iAttr(igrow,1)
!  write(6,*) subName, "sMat%data%iAttr(igcol, 1) = ",sMat%data%iAttr(igcol,1)
!  write(6,*) subName, "sMat%data%iAttr(igrow,ns) = ",sMat%data%iAttr(igrow,ns)
!  write(6,*) subName, "sMat%data%iAttr(igcol,ns) = ",sMat%data%iAttr(igcol,ns)
!
!  !----------------------------------------------------------------------------
!  ! sort and re-test
!  !----------------------------------------------------------------------------
!
!  write(6,F00) "Ascending sort test..."
!
!  call cpl_mct_sMat_SortPermute(sMat,sortKey)
!
!  write(6,*) subName, "sMat%data%iAttr(igrow, 1) = ",sMat%data%iAttr(igrow,1)
!  write(6,*) subName, "sMat%data%iAttr(igcol, 1) = ",sMat%data%iAttr(igcol,1)
!  write(6,*) subName, "sMat%data%iAttr(igrow,ns) = ",sMat%data%iAttr(igrow,ns)
!  write(6,*) subName, "sMat%data%iAttr(igcol,ns) = ",sMat%data%iAttr(igcol,ns)

end subroutine cpl_map_test

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNone - just return, do nothing.
!
! !DESCRIPTION:
!    Stub for north pole correction of vector fields.  This one just returns
!    and does no adjustments.
!
! !REVISION HISTORY:
!    1Feb03 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNone(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP

   !--- local ---
   logical,save :: first_call = .true. ! flags 1st invocation of routine

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNone) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNone) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) ' WARNING, no npFix adjustments '
     first_call = .false.
   endif

end subroutine cpl_map_npFixNone

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNew - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered and equally spaced.
!
! !REVISION HISTORY:
!    20Sep02 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNew(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod
   use shr_timer_mod       ! share timer routines

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP

   !--- local ---
   integer(IN)  :: n,m                   ! generic indices
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kin                   ! index index
   integer(IN)  :: nmin,nmax             ! indices of highest latitude in input
   integer(IN)  :: npts                  ! local number of points in an aV
   integer(IN)  :: num                   ! number of points at highest latitude
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   real(R8)     :: latmax                ! value of highest latitude
   real(R8)     :: olon,olat             ! output bundle lon/lat
   real(R8)     :: ilon,ilat             ! input bundle lon/lat
   real(R8)     :: npu,npv               ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2         ! angles for trig functions
   integer(IN),allocatable :: nn1(:)     ! index of highest latitude
   integer(IN),allocatable :: nn2(:)     ! index of highest latitude
   real(R8),allocatable :: ilon1(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat1(:)      ! lat of input grid at highest latitude
   real(R8),allocatable :: ilon2(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat2(:)      ! lat of input grid at highest latitude
   real(R8)     :: w1,w2,w3,w4           ! weights
   real(R8)     :: f1,f2,f3,f4           ! function values
   real(R8)     :: alpha,beta            ! used to generate weights
   real(R8)     :: rtmp                  ! real temporary
   type(cpl_mct_aVect) :: gData          ! global/gathered input bundle data
   type(cpl_mct_aVect) :: gGrid          ! global/gathered input bundle data
   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   logical      :: found                 ! search for new interpolation
   integer(IN)  :: rcode                 ! error code
   integer(IN)  :: nn1x                  ! tmp
   real(R8)     :: ilon1x                ! tmp
   logical,save :: first_call = .true.   ! flags 1st invocation of routine
   integer(IN),save :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNew) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNew) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) " compute bilinear weights & indicies for NP region."
     call shr_timer_get(tnpf1,'tnpf1')
     call shr_timer_get(tnpf2,'tnpf2')
     call shr_timer_get(tnpf3,'tnpf3')
     call shr_timer_get(tnpf4,'tnpf4')
     call shr_timer_get(tnpf5,'tnpf5')
     call shr_timer_get(tnpf6,'tnpf6')
     call shr_timer_get(tnpf7,'tnpf7')
     call shr_timer_get(tnpf8,'tnpf8')
     call shr_timer_get(tnpf9,'tnpf9')
     first_call = .false.
   endif

   call shr_timer_start(tnpf1)

   kui = cpl_mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = cpl_mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = cpl_mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = cpl_mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   kin   = cpl_mct_aVect_indexRA(buni%dom%lGrid,"index",perrWith=subName)
   klati = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   nmin = (buni%dom%ni)*(buni%dom%nj-1) + 1
   nmax = buni%dom%n
   num  = buni%dom%ni

   call shr_timer_stop(tnpf1)
   call shr_timer_start(tnpf2)

   call cpl_mct_aVect_gather(buni%data,gData,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call cpl_mct_aVect_gather(buni%dom%lGrid,gGrid,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)

   call shr_timer_stop(tnpf2)

!--artificial tests, fill gdata with "known" data here
!   if (cpl_comm_comp_pid == 0) then
!   do n = 1,num
!     nn1x = nmin + n - 1
!     rtmp = gGrid%rAttr(kloni,nn1x)
!     ilon1x = mod(rtmp,360.)
!     theta1 = ilon1x*cpl_const_deg2rad
!     gData%rAttr(kui,nn1x) = 1.0
!     gData%rAttr(kvi,nn1x) = 0.0
!   enddo
!   endif

   call shr_timer_start(tnpf3)

   if (cpl_comm_comp_pid /= pid0) call cpl_mct_aVect_clean(gData)
   if (cpl_comm_comp_pid /= pid0) call cpl_mct_aVect_clean(gGrid)
   call cpl_mct_aVect_bcast(gData,pid0,cpl_comm_comp,rcode)
   call cpl_mct_aVect_bcast(gGrid,pid0,cpl_comm_comp,rcode)

   call shr_timer_stop(tnpf3)
   call shr_timer_start(tnpf4)

   allocate(nn1(num))
   allocate(nn2(num))
   allocate(ilon1(num))
   allocate(ilon2(num))
   allocate(ilat1(num))
   allocate(ilat2(num))

   call shr_timer_stop(tnpf4)
   call shr_timer_start(tnpf5)

   latmax = gGrid%rAttr(klati,nmin)
   npu = 0.
   npv = 0.
   do n = 1,num
     nn1(n) = nmin + n - 1
     nn2(n) = nn1(n) + 1
     if (nn2(n) > nmax) nn2(n) = nn2(n) - num

     rtmp = gGrid%rAttr(kloni,nn1(n))
     ilon1(n) = mod(rtmp+360.,360.)
     rtmp = gGrid%rAttr(kloni,nn2(n))
     ilon2(n) = mod(rtmp+360.,360.)
     ilat1(n) =     gGrid%rAttr(klati,nn1(n))
     ilat2(n) =     gGrid%rAttr(klati,nn2(n))
     if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360.

     latmax = max(latmax,ilat1(n))

     theta1 = ilon1(n)*cpl_const_deg2rad
     npu = npu + cos(theta1)*gData%rAttr(kui,nn1(n)) &
               - sin(theta1)*gData%rAttr(kvi,nn1(n))
     npv = npv + sin(theta1)*gData%rAttr(kui,nn1(n)) &
               + cos(theta1)*gData%rAttr(kvi,nn1(n))
   enddo
   npu = npu / float(num)
   npv = npv / float(num)

   call shr_timer_stop(tnpf5)
   call shr_timer_start(tnpf6)

   npts = cpl_mct_aVect_lSize(buno%data)
   do m = 1,npts
     olat = buno%dom%lGrid%rAttr(klato,m)
     if (olat >= latmax) then
       rtmp = buno%dom%lGrid%rAttr(klono,m)
       olon = mod(rtmp,360.)
       n = 1
       found = .false.
       do while (n <= num .and. .not.found )
         if (    olon >= ilon1(n) .and. olon < ilon2(n) .or.   &
            olon+360. >= ilon1(n) .and. olon < ilon2(n)) then
           ilat = (ilat1(n) + ilat2(n)) * 0.5
           if (ilon2(n) == ilon1(n)) then
             alpha = 0.5
           else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
           else if (olon+360.>= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon+360. - ilon1(n)) / (ilon2(n) - ilon1(n))
           else
             write(6,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
           endif
           if (ilat >= 90.) then
             beta  = 1.0
           else
             beta  = (olat - ilat) / (90. - ilat)
           endif
           w1 = (1.0-alpha)*(1.0-beta)
           w2 = (    alpha)*(1.0-beta)
           w3 = (    alpha)*(    beta)
           w4 = (1.0-alpha)*(    beta)

           theta1 = ilon1(n)*cpl_const_deg2rad
           theta2 = ilon2(n)*cpl_const_deg2rad

           f1 = gData%rAttr(kui,nn1(n))  ! 4 input velocities
           f2 = gData%rAttr(kui,nn2(n))
           f3 =  cos(theta1)*npu + sin(theta1)*npv
           f4 =  cos(theta2)*npu + sin(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

           f1 = gData%rAttr(kvi,nn1(n))
           f2 = gData%rAttr(kvi,nn2(n))
           f3 = -sin(theta1)*npu + cos(theta1)*npv
           f4 = -sin(theta2)*npu + cos(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
           found = .true.
         endif
         n = n + 1     ! normal increment
       enddo
       if ( .not.found ) then
         write(6,*) subName,' ERROR: found = false ',found,m,olon,olat
       endif
     endif
   end do

   call shr_timer_stop(tnpf6)
   call shr_timer_start(tnpf7)

   call cpl_mct_aVect_clean(gData)
   call cpl_mct_aVect_clean(gGrid)
   deallocate(nn1)
   deallocate(nn2)
   deallocate(ilon1)
   deallocate(ilon2)
   deallocate(ilat1)
   deallocate(ilat2)

   call shr_timer_stop(tnpf7)

end subroutine cpl_map_npFixNew

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNew2 - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered and equally spaced.
!
!    This version (New2) is the same as New except it saves the gGrid data
!    type from the first call.  This assumes the gGrid used in all calls
!    to npfix is the same.  This is bfb with version (New) on 9/1/2003.
!
! !REVISION HISTORY:
!    29Aug03 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNew2(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod
   use shr_timer_mod       ! share timer routines

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP

   !--- local ---
   integer(IN)  :: n,m                   ! generic indices
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kin                   ! index index
   integer(IN)  :: nmin,nmax             ! indices of highest latitude in input
   integer(IN)  :: npts                  ! local number of points in an aV
   integer(IN)  :: num                   ! number of points at highest latitude
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   real(R8)     :: latmax                ! value of highest latitude
   real(R8)     :: olon,olat             ! output bundle lon/lat
   real(R8)     :: ilon,ilat             ! input bundle lon/lat
   real(R8)     :: npu,npv               ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2         ! angles for trig functions
   integer(IN),allocatable :: nn1(:)     ! index of highest latitude
   integer(IN),allocatable :: nn2(:)     ! index of highest latitude
   real(R8),allocatable :: ilon1(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat1(:)      ! lat of input grid at highest latitude
   real(R8),allocatable :: ilon2(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat2(:)      ! lat of input grid at highest latitude
   real(R8)     :: w1,w2,w3,w4           ! weights
   real(R8)     :: f1,f2,f3,f4           ! function values
   real(R8)     :: alpha,beta            ! used to generate weights
   real(R8)     :: rtmp                  ! real temporary
   type(cpl_mct_aVect) :: gData          ! global/gathered input bundle data
   type(cpl_mct_aVect),save :: gGrid     ! global/gathered input bundle data
   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   logical      :: found                 ! search for new interpolation
   integer(IN)  :: rcode                 ! error code
   integer(IN)  :: nn1x                  ! tmp
   real(R8)     :: ilon1x                ! tmp
   logical,save :: first_call = .true.   ! flags 1st invocation of routine
   integer(IN),save :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNew2) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNew2) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) " compute bilinear weights & indicies for NP region."
     call shr_timer_get(tnpf1,'tnpf1')
     call shr_timer_get(tnpf2,'tnpf2')
     call shr_timer_get(tnpf3,'tnpf3')
     call shr_timer_get(tnpf4,'tnpf4')
     call shr_timer_get(tnpf5,'tnpf5')
     call shr_timer_get(tnpf6,'tnpf6')
     call shr_timer_get(tnpf7,'tnpf7')
     call shr_timer_get(tnpf8,'tnpf8')
     call shr_timer_get(tnpf9,'tnpf9')
     call cpl_mct_aVect_gather(buni%dom%lGrid,gGrid,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
     if (cpl_comm_comp_pid /= pid0) call cpl_mct_aVect_clean(gGrid)
     call cpl_mct_aVect_bcast(gGrid,pid0,cpl_comm_comp,rcode)
     first_call = .false.
   endif

   call shr_timer_start(tnpf1)

   kui = cpl_mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = cpl_mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = cpl_mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = cpl_mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   kin   = cpl_mct_aVect_indexRA(buni%dom%lGrid,"index",perrWith=subName)
   klati = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   nmin = (buni%dom%ni)*(buni%dom%nj-1) + 1
   nmax = buni%dom%n
   num  = buni%dom%ni

   call shr_timer_stop(tnpf1)
   call shr_timer_start(tnpf2)

   call cpl_mct_aVect_gather(buni%data,gData,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call shr_timer_stop(tnpf2)

!--artificial tests, fill gdata with "known" data here
!   if (cpl_comm_comp_pid == 0) then
!   do n = 1,num
!     nn1x = nmin + n - 1
!     rtmp = gGrid%rAttr(kloni,nn1x)
!     ilon1x = mod(rtmp,360.)
!     theta1 = ilon1x*cpl_const_deg2rad
!     gData%rAttr(kui,nn1x) = 1.0
!     gData%rAttr(kvi,nn1x) = 0.0
!   enddo
!   endif

   call shr_timer_start(tnpf3)

   if (cpl_comm_comp_pid /= pid0) call cpl_mct_aVect_clean(gData)
   call cpl_mct_aVect_bcast(gData,pid0,cpl_comm_comp,rcode)

   call shr_timer_stop(tnpf3)
   call shr_timer_start(tnpf4)

   allocate(nn1(num))
   allocate(nn2(num))
   allocate(ilon1(num))
   allocate(ilon2(num))
   allocate(ilat1(num))
   allocate(ilat2(num))

   call shr_timer_stop(tnpf4)
   call shr_timer_start(tnpf5)

   latmax = gGrid%rAttr(klati,nmin)
   npu = 0.
   npv = 0.
   do n = 1,num
     nn1(n) = nmin + n - 1
     nn2(n) = nn1(n) + 1
     if (nn2(n) > nmax) nn2(n) = nn2(n) - num

     rtmp = gGrid%rAttr(kloni,nn1(n))
     ilon1(n) = mod(rtmp+360.,360.)
     rtmp = gGrid%rAttr(kloni,nn2(n))
     ilon2(n) = mod(rtmp+360.,360.)
     ilat1(n) =     gGrid%rAttr(klati,nn1(n))
     ilat2(n) =     gGrid%rAttr(klati,nn2(n))
     if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360.

     latmax = max(latmax,ilat1(n))

     theta1 = ilon1(n)*cpl_const_deg2rad
     npu = npu + cos(theta1)*gData%rAttr(kui,nn1(n)) &
               - sin(theta1)*gData%rAttr(kvi,nn1(n))
     npv = npv + sin(theta1)*gData%rAttr(kui,nn1(n)) &
               + cos(theta1)*gData%rAttr(kvi,nn1(n))
   enddo
   npu = npu / float(num)
   npv = npv / float(num)

   call shr_timer_stop(tnpf5)
   call shr_timer_start(tnpf6)

   npts = cpl_mct_aVect_lSize(buno%data)
   do m = 1,npts
     olat = buno%dom%lGrid%rAttr(klato,m)
     if (olat >= latmax) then
       rtmp = buno%dom%lGrid%rAttr(klono,m)
       olon = mod(rtmp,360.)
       n = 1
       found = .false.
       do while (n <= num .and. .not.found )
         if (    olon >= ilon1(n) .and. olon < ilon2(n) .or.   &
            olon+360. >= ilon1(n) .and. olon < ilon2(n)) then
           ilat = (ilat1(n) + ilat2(n)) * 0.5
           if (ilon2(n) == ilon1(n)) then
             alpha = 0.5
           else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
           else if (olon+360.>= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon+360. - ilon1(n)) / (ilon2(n) - ilon1(n))
           else
             write(6,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
           endif
           if (ilat >= 90.) then
             beta  = 1.0
           else
             beta  = (olat - ilat) / (90. - ilat)
           endif
           w1 = (1.0-alpha)*(1.0-beta)
           w2 = (    alpha)*(1.0-beta)
           w3 = (    alpha)*(    beta)
           w4 = (1.0-alpha)*(    beta)

           theta1 = ilon1(n)*cpl_const_deg2rad
           theta2 = ilon2(n)*cpl_const_deg2rad

           f1 = gData%rAttr(kui,nn1(n))  ! 4 input velocities
           f2 = gData%rAttr(kui,nn2(n))
           f3 =  cos(theta1)*npu + sin(theta1)*npv
           f4 =  cos(theta2)*npu + sin(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

           f1 = gData%rAttr(kvi,nn1(n))
           f2 = gData%rAttr(kvi,nn2(n))
           f3 = -sin(theta1)*npu + cos(theta1)*npv
           f4 = -sin(theta2)*npu + cos(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
           found = .true.
         endif
         n = n + 1     ! normal increment
       enddo
       if ( .not.found ) then
         write(6,*) subName,' ERROR: found = false ',found,m,olon,olat
       endif
     endif
   end do

   call shr_timer_stop(tnpf6)
   call shr_timer_start(tnpf7)

   call cpl_mct_aVect_clean(gData)
   deallocate(nn1)
   deallocate(nn2)
   deallocate(ilon1)
   deallocate(ilon2)
   deallocate(ilat1)
   deallocate(ilat2)

   call shr_timer_stop(tnpf7)

   call shr_sys_flush(6)
!  call shr_sys_abort()

end subroutine cpl_map_npFixNew2

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixNew3 - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered and equally spaced.
!
!    This version (New3) is the same as New except it saves the gGrid data
!    type from the first call.  This assumes the gGrid used in all calls
!    to npfix is the same.  It is different from New2 in that it doesn't
!    use gather to compute npu and npv.  This is bfb with New2 and New on 
!    9/1/2003.
!
! !REVISION HISTORY:
!    29Aug03 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixNew3(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod
   use shr_timer_mod       ! share timer routines

#if (! defined HIDE_MPI)
#include <mpif.h>         ! mpi library include file
#endif

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(out)  :: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP

   !--- local ---
   integer(IN)  :: n,m                   ! generic indices
   integer(IN)  :: n1,n2,n3              ! generic indices
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kin                   ! index index
   integer(IN)  :: nmin,nmax             ! indices of highest latitude in input
   integer(IN)  :: npts                  ! local number of points in an aV
   integer(IN)  :: num                   ! number of points at highest latitude
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   integer(IN)  :: index                 ! index value
   integer(IN)  :: nn1,nn2               ! local n to global n values
   real(R8)     :: rindex                ! index value
   real(R8)     :: latmax                ! value of highest latitude
   real(R8)     :: olon,olat             ! output bundle lon/lat
   real(R8)     :: ilon,ilat             ! input bundle lon/lat
   real(R8)     :: npu,npv               ! np velocity fields relative to lon
   real(R8)     :: theta1,theta2         ! angles for trig functions
   real(R8),allocatable :: ilon1(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat1(:)      ! lat of input grid at highest latitude
   real(R8),allocatable :: ilon2(:)      ! lon of input grid at highest latitude
   real(R8),allocatable :: ilat2(:)      ! lat of input grid at highest latitude
   real(R8)     :: w1,w2,w3,w4           ! weights
   real(R8)     :: f1,f2,f3,f4           ! function values
   real(R8)     :: alpha,beta            ! used to generate weights
   real(R8)     :: rtmp                  ! real temporary
   real(R8),allocatable :: lData(:,:)    ! last lat local input bundle data
                                         ! also compressed global data
   real(R8),allocatable :: gData(:,:,:)  ! last lat gathered input bundle data
   type(cpl_mct_aVect),save :: gGrid     ! global/gathered input bundle data
   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   logical      :: found                 ! search for new interpolation
   integer(IN)  :: rcode                 ! error code
   integer(IN)  :: np1                   ! n+1 or tmp
   real(R8)     :: ilon1x                ! tmp
   logical,save :: first_call = .true.   ! flags 1st invocation of routine
   integer(IN),save :: tnpf1,tnpf2,tnpf3,tnpf4,tnpf5,tnpf6,tnpf7,tnpf8,tnpf9

   !--- formats ---
   character(*),parameter :: subName = '(cpl_map_npFixNew3) '
   character(*),parameter :: F00 = "('(cpl_map_npFixNew3) ',8a)"

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

   if (first_call) then
     write(6,F00) " compute bilinear weights & indicies for NP region."
     call shr_timer_get(tnpf1,'tnpf1')
     call shr_timer_get(tnpf2,'tnpf2')
     call shr_timer_get(tnpf3,'tnpf3')
     call shr_timer_get(tnpf4,'tnpf4')
     call shr_timer_get(tnpf5,'tnpf5')
     call shr_timer_get(tnpf6,'tnpf6')
     call shr_timer_get(tnpf7,'tnpf7')
     call shr_timer_get(tnpf8,'tnpf8')
     call shr_timer_get(tnpf9,'tnpf9')
     call cpl_mct_aVect_gather(buni%dom%lGrid,gGrid,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
     if (cpl_comm_comp_pid /= pid0) call cpl_mct_aVect_clean(gGrid)
     call cpl_mct_aVect_bcast(gGrid,pid0,cpl_comm_comp,rcode)
     first_call = .false.
   endif

   call shr_timer_start(tnpf1)

   kui = cpl_mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = cpl_mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = cpl_mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = cpl_mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   kin   = cpl_mct_aVect_indexRA(buni%dom%lGrid,"index",perrWith=subName)
   klati = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   nmin = (buni%dom%ni)*(buni%dom%nj-1) + 1
   nmax = buni%dom%n
   num  = buni%dom%ni

   call shr_timer_stop(tnpf1)
   call shr_timer_start(tnpf2)

!  barrier not required but interesting for timing.
!  call shr_mpi_barrier(cpl_comm_comp,subName//" barrier")

   call shr_timer_stop(tnpf2)

   allocate(lData(3,num))
   allocate(gData(3,num,cpl_comm_comp_npe))
   lData = 0.
   gData = 0.
   npts = cpl_mct_aVect_lSize(buni%data)
   m = 0   
   do n=1,npts
     rindex = buni%dom%lGrid%rAttr(kin,n)
     if (rindex.ge.nmin) then
       m=m+1
       lData(1,m) = rindex
       lData(2,m) = buni%data%rAttr(kui,n)
       lData(3,m) = buni%data%rAttr(kvi,n)
     endif
   enddo

   call MPI_ALLGATHER(lData,3*num,MPI_REAL8, &
     gData,3*num,MPI_REAL8,cpl_comm_comp,rcode)

   if (rcode.ne.0) then
     write(6,*) trim(subName),' rcode error ',rcode
     call shr_sys_abort()
   endif

   call shr_timer_start(tnpf3)

   m = 0
   lData = 0.
   do n2=1,num
   do n3=1,cpl_comm_comp_npe
     if (gData(1,n2,n3).gt.0.1) then
       m = m+1
       index = nint(gData(1,n2,n3)) - nmin + 1
       lData(1:3,index) = gData(1:3,n2,n3)
     endif
   enddo
   enddo
   if (m.ne.num) write(6,*) trim(subName),' error allgather ',m,num
   do n2=1,num
     if (lData(1,n2).lt.0.1) then
       write(6,*) trim(subName),' error allgather2 ',n2
     endif
   enddo

   call shr_timer_stop(tnpf3)
   call shr_timer_start(tnpf4)

   allocate(ilon1(num))
   allocate(ilon2(num))
   allocate(ilat1(num))
   allocate(ilat2(num))

   call shr_timer_stop(tnpf4)
   call shr_timer_start(tnpf5)

   latmax = gGrid%rAttr(klati,nmin)
   npu = 0.
   npv = 0.
   do n = 1,num
     np1 = mod(n,num)+1
     nn1 = nmin + n - 1
     nn2 = nmin + np1 - 1 
     rtmp = gGrid%rAttr(kloni,nn1)
     ilon1(n) = mod(rtmp+360.,360.)
     rtmp = gGrid%rAttr(kloni,nn2)
     ilon2(n) = mod(rtmp+360.,360.)
     ilat1(n) =     gGrid%rAttr(klati,nn1)
     ilat2(n) =     gGrid%rAttr(klati,nn2)
     if (ilon2(n) < ilon1(n)) ilon2(n) = ilon2(n) + 360.

     latmax = max(latmax,ilat1(n))

     theta1 = ilon1(n)*cpl_const_deg2rad
     npu = npu + cos(theta1)*lData(2,n) &
               - sin(theta1)*lData(3,n)
     npv = npv + sin(theta1)*lData(2,n) &
               + cos(theta1)*lData(3,n)
   enddo
   npu = npu / float(num)
   npv = npv / float(num)

   call shr_timer_stop(tnpf5)
   call shr_timer_start(tnpf6)

   npts = cpl_mct_aVect_lSize(buno%data)
   do m = 1,npts
     olat = buno%dom%lGrid%rAttr(klato,m)
     if (olat >= latmax) then
       rtmp = buno%dom%lGrid%rAttr(klono,m)
       olon = mod(rtmp,360.)
       n = 1
       found = .false.
       do while (n <= num .and. .not.found )
         if (    olon >= ilon1(n) .and. olon < ilon2(n) .or.   &
            olon+360. >= ilon1(n) .and. olon < ilon2(n)) then
           np1 = mod(n,num)+1
           ilat = (ilat1(n) + ilat2(n)) * 0.5
           if (ilon2(n) == ilon1(n)) then
             alpha = 0.5
           else if (    olon >= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon - ilon1(n)) / (ilon2(n) - ilon1(n))
           else if (olon+360.>= ilon1(n) .and. olon < ilon2(n)) then
             alpha = (olon+360. - ilon1(n)) / (ilon2(n) - ilon1(n))
           else
             write(6,*) subName,' ERROR: olon ',olon,ilon1(n),ilon2(n)
           endif
           if (ilat >= 90.) then
             beta  = 1.0
           else
             beta  = (olat - ilat) / (90. - ilat)
           endif
           w1 = (1.0-alpha)*(1.0-beta)
           w2 = (    alpha)*(1.0-beta)
           w3 = (    alpha)*(    beta)
           w4 = (1.0-alpha)*(    beta)

           theta1 = ilon1(n)*cpl_const_deg2rad
           theta2 = ilon2(n)*cpl_const_deg2rad

           f1 = lData(2,n)
           f2 = lData(2,np1)
           f3 =  cos(theta1)*npu + sin(theta1)*npv
           f4 =  cos(theta2)*npu + sin(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kuo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4

           f1 = lData(3,n)
           f2 = lData(3,np1)
           f3 = -sin(theta1)*npu + cos(theta1)*npv
           f4 = -sin(theta2)*npu + cos(theta2)*npv
           rtmp = w1*f1 + w2*f2 + w3*f3 + w4*f4
           buno%data%rAttr(kvo,m) = w1*f1 + w2*f2 + w3*f3 + w4*f4
           found = .true.
         endif
         n = n + 1     ! normal increment
       enddo
       if ( .not.found ) then
         write(6,*) subName,' ERROR: found = false ',found,m,olon,olat
       endif
     endif
   end do

   call shr_timer_stop(tnpf6)
   call shr_timer_start(tnpf7)

   deallocate(gData)
   deallocate(lData)
   deallocate(ilon1)
   deallocate(ilon2)
   deallocate(ilat1)
   deallocate(ilat2)

   call shr_timer_stop(tnpf7)

end subroutine cpl_map_npFixNew3

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_map_npFixOld - correct the north pole mapping of velocity fields
!
! !DESCRIPTION:
!    Correct the north pole mapping of velocity fields from the atm to ocn
!    grids.  This assumes the input grid is a regular lat/lon with the north
!    pole surrounded by the last latitude line in the input array.  The
!    longitudes in the last latitude must be ordered.
!
! !REVISION HISTORY:
!    20Sep02 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_map_npFixOld(buni,buno,fld1,fld2)

! !USES:

   use cpl_const_mod

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout):: buni    ! input bundle
   type(cpl_bundle),intent(inout):: buno    ! output bundle
   character(*)    ,intent(in)   :: fld1    ! name of first input field
   character(*)    ,intent(in)   :: fld2    ! name of second input field

!EOP

   !----- local -----
   integer :: n           ! index relative to output array (subroutine arg)
   integer :: npts        ! number of output array points that need fixing
   integer :: nn          ! index that goes from 1 to npts
   integer :: m           ! index relative to input array (subroutine arg)
   integer :: mpts        ! number of grid points in top row of input grid
   integer :: mm          ! index that goes from 1 to mpts
   integer :: m1st        ! 1st m st y(m)=max input latitude

   real    :: npu         ! north-pole u-velocity
   real    :: npv         ! north-pole v-velocity

   real   ,allocatable ::  wght(:,:) ! 4 bilin map weights
   integer,allocatable ::    mw(:,:) ! index on input grid
   integer,allocatable ::    nw  (:) ! index on output grid
   real   ,allocatable :: Sa_in_u(:) ! local NP u-velocity "on" input grid
   real   ,allocatable :: Sa_in_v(:) ! local NP v-velocity "on" input grid
   real   ,allocatable ::  theta (:) ! rotation angle wrt 0 deg east (radians)

   real    :: alpha,beta    ! used to compute bilin weights
   real    :: dx,dy         ! used to compute biline weights
   real    :: xm,ym         ! x/y-coords of output data for bilin weights
   real    :: xl,xh,yl,yh   ! x/y-coords of  input data for bilin weights
   real    :: w1,w2,w3,w4   ! 4 bilin weights
   integer :: m1,m2,m3,m4   ! 4 bilin indicies wrt input grid
   integer :: mm3,mm4       ! bilin indicies wrt fabricated NP input grid row
   real    :: f1,f2,f3,f4   ! 4 bilin function values on input grid

   logical :: first_call = .true. ! flags 1st invocation of routine
   integer :: timer_n             ! system-clock-timer number

   type(cpl_mct_aVect) :: gDatai      ! global/gathered input bundle data
   type(cpl_mct_aVect) :: gGridi      ! global/gathered input bundle data
   type(cpl_mct_aVect) :: gDatao      ! global/gathered output bundle data
   type(cpl_mct_aVect) :: gGrido      ! global/gathered output bundle data

   integer(IN),parameter :: pid0 = 0     ! root process pid = zero
   integer(IN)  :: n_a,n_o               ! a/o input/output sizes
   integer(IN)  :: kui,kvi               ! field indices
   integer(IN)  :: kuo,kvo               ! field indices
   integer(IN)  :: kloni                 ! longitude index on input domain
   integer(IN)  :: klati                 ! latitude index on input domain
   integer(IN)  :: klono                 ! longitude index on output domain
   integer(IN)  :: klato                 ! latitude index on output domain
   integer(IN)  :: rcode                 ! error code

   save  

   !----- formats -----
   character(*),parameter :: subName = '(cpl_map_npFixOld) '
   character(*),parameter :: F00 = "('(cpl_map_npFixOld) ',8a)"

!-------------------------------------------------------------------------------
! PURPOSE:
!    Compute/correct atm wind velocties on ocn (or ice) grid.
!
! ASSUMPTIONS:
!  o the input (atm) grid has a "top row" with a constant latitude value
!  o the ocn & ice have identical grids, this routine works for both ocn & ice
!
! NOTATION:
!
!           ^  f4      f3
!       b   |
!       e  dy      f                interpolate to f from f1,f2,f3,f4
!       t   |                     
!       a   v  f1      f2
!              <--dx--->
!               alpha
!
!   f  is located at (xm,ym)
!   f1 is located at (xl,yl)
!   f2 is located at (xh,yl)
!   f3 is located at (xh,yh)
!   f4 is located at (xl,yh)
!
!-------------------------------------------------------------------------------

   call cpl_mct_aVect_gather(buni%data,gDatai,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call cpl_mct_aVect_gather(buni%dom%lGrid,gGridi,buni%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call cpl_mct_aVect_gather(buno%data,gDatao,buno%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call cpl_mct_aVect_gather(buno%dom%lGrid,gGrido,buno%dom%gsMap,pid0,cpl_comm_comp,rcode)


if (cpl_comm_comp_pid == pid0 ) then

   n_a = buni%dom%n
   n_o = buno%dom%n

   kui = cpl_mct_aVect_indexRA(buni%data,fld1,perrWith=subName)
   kvi = cpl_mct_aVect_indexRA(buni%data,fld2,perrWith=subName)
   kuo = cpl_mct_aVect_indexRA(buno%data,fld1,perrWith=subName)
   kvo = cpl_mct_aVect_indexRA(buno%data,fld2,perrWith=subName)

   klati = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lat",perrWith=subName)
   kloni = cpl_mct_aVect_indexRA(buni%dom%lGrid,"lon",perrWith=subName)
   klato = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lat",perrWith=subName)
   klono = cpl_mct_aVect_indexRA(buno%dom%lGrid,"lon",perrWith=subName)

   !----------------------------------------------------------------------------
   ! one-time setup computations
   !----------------------------------------------------------------------------
   if (first_call) then
      write(6,F00) " compute bilinear weights & indicies for NP region."

      !--- determine number of points in top row of atm (input) data ---
      mpts = 0
      m1st = 0
      do m=1,n_a
        if (gGridi%rAttr(klati,m) >= gGridi%rAttr(klati,n_a)) then
           m1st = m               ! 1st point in top row
           mpts = n_a - m1st + 1  ! number of points in top row
           exit
        end if 
      end do

      !--- allocate memory for a new, north-pole row ---
      allocate(Sa_in_u(mpts))  ! u-velocity
      allocate(Sa_in_v(mpts))  ! v-velocity
      allocate( theta (mpts))  ! rotation angle wrt prime-meridian

      !--- compute (u,v) rotation angles relative to prime-meridian ---
      do m=m1st,n_a
        mm = m - m1st + 1                   ! index into local NP row
        theta(mm) = gGridi%rAttr(kloni,m)*cpl_const_deg2rad ! units are radians
      end do

      !--- determine number of points that need correcting (output) data ---
      npts = 0
      do n=1,n_o
        if (gGrido%rAttr(klato,n) > gGridi%rAttr(klati,n_a)) npts = npts + 1
      end do

      !--- allocate bilin-map memory given number of output points involved ---
      allocate(   nw(  npts) ) ! indices to output points
      allocate(   mw(4,npts) ) ! for each output point, 4 input  indicies
      allocate( wght(4,npts) ) ! for each output point, 4 bilinear weights

      !--- compute bilin-map data that utilizes new/local NP row ---
      nn = 0
      do n=1,n_o
        if (gGrido%rAttr(klato,n) > gGridi%rAttr(klati,n_a)) then
           nn = nn + 1
           nw(nn) = n              ! output index associated with weights
           xm = gGrido%rAttr(klono,n)     ! x-coord of output
           ym = gGrido%rAttr(klato,n)     ! y-coord of output
           yl = gGridi%rAttr(klati,n_a)   ! constant y-coord of top atm row
           yh = 90.0               ! constant y-coord of north pole
           if ( xm < gGridi%rAttr(kloni,m1st)) then
             m1 = n_a              ! input indicies associated with weights
             m2 = m1st
             m3 = m2   + mpts      ! note: this index goes beyond actual data
             m4 = m1   + mpts      ! note: this index goes beyond actual data
             xl = gGridi%rAttr(kloni,m1) - 360.0
             xh = gGridi%rAttr(kloni,m2)
           else if ( xm < gGridi%rAttr(kloni,n_a) ) then
             do m=m1st,n_a - 1
               if ( xm < gGridi%rAttr(kloni,m+1) ) then
                 m1 = m            ! input indicies associated with weights
                 m2 = m+1
                 m3 = m2  + mpts   ! note: this index goes beyond actual data
                 m4 = m1  + mpts   ! note: this index goes beyond actual data
                 xl = gGridi%rAttr(kloni,m1)
                 xh = gGridi%rAttr(kloni,m2)
                 exit
               end if
             end do
           else
             m1 = n_a              ! input indicies associated with weights
             m2 = m1st
             m3 = m2   + mpts      ! note: this index goes beyond actual data
             m4 = m1   + mpts      ! note: this index goes beyond actual data
             xl = gGridi%rAttr(kloni,m1)
             xh = gGridi%rAttr(kloni,m2) + 360.0
           end if
           dx = xh - xl 
           dy = yh - yl 
           alpha = (xm - xl)/dx
           beta  = (ym - yl)/dy
           wght(1,nn) = (1.0-alpha)*(1.0-beta) ! 4 input bilinear weights
           wght(2,nn) = (    alpha)*(1.0-beta)
           wght(3,nn) = (    alpha)*(    beta)
           wght(4,nn) = (1.0-alpha)*(    beta) 
           mw  (1,nn) = m1                     ! input indicies wrt with weights
           mw  (2,nn) = m2   
           mw  (3,nn) = m3  
           mw  (4,nn) = m4 
        end if
      end do

      first_call = .false.
   end if

   !----------------------------------------------------------------------------
   ! fabricate multiple instances of atm north-pole data w/ various rotations
   !----------------------------------------------------------------------------

   npu = 0  ! avg NP value with y-basis vector pointing north along
   npv = 0  ! prime-meridian and x-basis vector pointing east at prime-meridian
   do mm=1,mpts
      m = n_a - mpts + mm
      npu = npu + cos(theta(mm))*gDatai%rAttr(kui,m) - sin(theta(mm))*gDatai%rAttr(kvi,m)
      npv = npv + sin(theta(mm))*gDatai%rAttr(kui,m) + cos(theta(mm))*gDatai%rAttr(kvi,m)
   end do
   npu = npu/mpts
   npv = npv/mpts

   !--- copies of NP value with basis-vectors aligned with cell's longitude
   do mm=1,mpts
      Sa_in_u(mm) =  cos(theta(mm))*npu + sin(theta(mm))*npv
      Sa_in_v(mm) = -sin(theta(mm))*npu + cos(theta(mm))*npv
   end do

   !----------------------------------------------------------------------------
   ! bilinear interpolation from atm to select set of ocn/ice grid points
   !----------------------------------------------------------------------------
   do nn=1,npts
      n   = nw(nn)            ! output grid index

      m1  = mw(1,nn)          ! input grid indicies on atm grid
      m2  = mw(2,nn)          ! input grid indicies on atm grid
      mm3 = mw(3,nn) - n_a    ! input grid indicies on local "NP grid" array
      mm4 = mw(4,nn) - n_a    ! input grid indicies on local "NP grid" array

      w1 = wght(1,nn)         ! 4 mapping weights
      w2 = wght(2,nn)
      w3 = wght(3,nn)
      w4 = wght(4,nn)

      f1 = gDatai%rAttr(kui,m1) ! 4 input velocities
      f2 = gDatai%rAttr(kui,m2)
      f3 = Sa_in_u(mm3)        ! local NP pole value
      f4 = Sa_in_u(mm4)        ! local NP pole value
      gDatao%rAttr(kuo,n) = w1*f1 + w2*f2 + w3*f3 + w4*f4

      f1 = gDatai%rAttr(kvi,m1)
      f2 = gDatai%rAttr(kvi,m2)
      f3 = Sa_in_v(mm3)        ! local NP pole value
      f4 = Sa_in_v(mm4)        ! local NP pole value
      gDatao%rAttr(kvo,n) = w1*f1 + w2*f2 + w3*f3 + w4*f4

   end do

endif

   call cpl_mct_aVect_clean(buno%data) 
   call cpl_mct_aVect_scatter(gDatao,buno%data,buno%dom%gsMap,pid0,cpl_comm_comp,rcode)

   if (cpl_comm_comp_pid == pid0 ) then
      call cpl_mct_aVect_clean(gDatai)
      call cpl_mct_aVect_clean(gGridi)
      call cpl_mct_aVect_clean(gDatao)
      call cpl_mct_aVect_clean(gGrido)
   endif

end subroutine cpl_map_npFixOld

!===============================================================================

end module cpl_map_mod

