!===============================================================================
! CVS: $Id: main.F90,v 1.57 2004/05/24 01:34:50 kauff Exp $
! CVS: $Source: /fs/cgd/csm/models/CVS.REPOS/cpl/cpl6/main.F90,v $
! CVS: $Name: ccsm3_0_rel04 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: main -- the main program
!
! !DESCRIPTION:
!    This is the cpl6 main program containing the main integration loop.
!
! !REMARKS:
!    This file provides a high-level, pseudo-code-like view of the coupler.
!
! !REVISION HISTORY:
!     2002-Dec-19 - T. Craig, code migrated from yak3 dev repo to CCSM repo
!     2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

program cpl

! !USES:

   use shr_date_mod       ! shr code: date/time data type & methods
   use shr_msg_mod        ! shr code: i/o redirect, chdir
   use shr_sys_mod        ! shr code: wrappers to system calls
   use shr_timer_mod      ! shr code: timing utilities
   use shr_mpi_mod        ! shr code: mpi layer

   use cpl_kind_mod       ! kinds for strong typing
   use cpl_mct_mod        ! wrapper/access to mct lib
   use cpl_comm_mod       ! MPI process & communicator group info
   use cpl_fields_mod     ! fields lists common to cpl & components
   use cpl_contract_mod   ! contract (message-passing) data type & methods
   ! steder: Directly call cpl_comm_init
   use cpl_comm_mod
   use cpl_interface_mod  ! contract/message-passing wrapper routines
   use cpl_infobuf_mod    ! information buffer module
   use cpl_control_mod    ! control flags/logic subsystem module
   use cpl_map_mod        ! mapping subsystem module
!  use cpl_iocdf_mod      ! netCDF file i/o routines (for debug data dumps)

   use data_mod           ! general          data declarations & routines
   use frac_mod           ! surface fraction data declarations & routines
   use areafact_mod       ! area correction  data declarations & routines
   use diag_mod           ! runtime diagnostic subsystem module
   use flux_mod           ! flux calculation   subsystem module
   use history_mod        ! history file i/o   subsystem module
   use restart_mod        ! restart file i/o   subsystem module
   use merge_mod          ! merging subsystem module
   use tStamp_mod         ! model date & wall clock time stamping
   use bitCheck_mod       ! useful to verify exact restart

!EOP

   implicit none
 
   !--- local ---
   integer(IN)          :: n            ! inner loop index
   integer(IN)          :: i,j,k        ! generic integers
   integer(IN)          :: rcode        ! routine return code
   integer(IN)          :: local_comm   ! local communicator ID
   integer(IN)          :: ncpl_a       ! number of atm communications per day
   integer(IN)          :: ncpl_i       ! number of ice communications per day
   integer(IN)          :: ncpl_l       ! number of lnd communications per day
   integer(IN)          :: ncpl_r       ! number of lnd/runoff communi per day
   integer(IN)          :: ncpl_o       ! number of ocn communications per day
   integer(IN)          :: lsize        ! size of local attrvect data
   integer(IN)          :: nfld         ! location of fld in attrvect data
   integer(IN)          :: nfld1        ! location of fld in attrvect data
   type(shr_date)       :: date         ! *the* official coupler date/time
   integer(IN)          :: year         ! model date's year
   integer(IN)          :: month        ! model date's month
   integer(IN)          :: day          ! model date's day 
   integer(IN)          :: sec          ! model date's seconds
   integer(IN)          :: cDate        ! model date's coded date (yymmss)
   character(8)         :: dstr         ! F90 wall clock date string yyyymmdd
   character(10)        :: tstr         ! F90 wall clock time string hhmmss.sss
   type(cpl_mct_aVect)  :: gData0       ! gathered data for optional land init
   type(cpl_mct_aVect)  :: gData1       ! gathered data for optional land init
   type(cpl_mct_aVect)  :: gData2       ! gathered data for optional land init
   real(R8),allocatable :: ifrac_i(:)   ! ice fraction
   integer(IN)          :: fid          ! debug netCDF file ID

   !--- for controlling vector-friendly mapping code
   logical              :: a2ovector, oi2avector, r2ovector

   !--- for atm-recv adjust for rain/snow, temporary until clm2.2 ---
   integer(IN) :: n1           ! loop index
   integer(IN) :: npts         ! size of AV
   real(R8)    :: rc,rl,sc,sl  ! rain & snow, convective & large-scale

   !--- for kludgy fix to problem with unmerged So_t values ---
   type(cpl_bundle):: fix_So2c_a ! alternate version of bun_So2c_a
   type(cpl_bundle):: fix_frac_a ! alternate version of bun_frac_a

   !--- for timers ---
   integer(IN)    ::     ti1,ti2,ti3,ti4,ti5,ti6,ti7,ti8,ti9 
   integer(IN)    ::     tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9 
   integer(IN)    :: t00,t01,t02,t03,t04,t05,t06,t07,t08,t09 
   integer(IN)    :: t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
   integer(IN)    :: t20,t21,t22,t23,t24,t25,t26,t27,t28,t29

   !--- formats ---
   character(*),parameter :: F00 = "('(main) ',8a)"
   character(*),parameter :: F01 = "('(main) ',a,5i6)"
   character(*),parameter :: F02 = "('(main) ',a,i4.4,2('-',i2.2),i6,'s')"
   character(*),parameter :: F03 = "('(main) ',a,i8.8)"
   character(*),parameter :: F12 = "('(main) date & time:',1x,&
   &                                     a4,2('-',a2),2x,a2,2(':',a2))"
   character(*),parameter :: F90 = "('(main) ',73('='))"
   character(*),parameter :: F91 = "('(main) ',73('-'))"
   character(*),parameter :: F92 = "('(main) ',73('-')/,'(main) ',a/,'(main) ',73('-'))"

   character(*),parameter :: modName = "(main) "
 
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !--- initialize timers ---
   call shr_timer_init
   call shr_timer_get(ti1,'ti1 - startup initialization')
   call shr_timer_get(t00,'t00 - main integration')
   call shr_timer_get(t01,'t01')
   call shr_timer_get(t02,'t02')
   call shr_timer_get(t03,'t03')
   call shr_timer_get(t04,'t04')
   call shr_timer_get(t05,'t05')
   call shr_timer_get(t06,'t06')
   call shr_timer_get(t07,'t07')
   call shr_timer_get(t08,'t08')
   call shr_timer_get(t09,'t09')
   call shr_timer_get(t10,'t10')
   call shr_timer_get(t11,'t11')
   call shr_timer_get(t12,'t12')
   call shr_timer_get(t13,'t13')
   call shr_timer_get(t14,'t14')
   call shr_timer_get(t15,'t15')
   call shr_timer_get(t16,'t16')
   call shr_timer_get(t17,'t17')
   call shr_timer_get(t18,'t18')
   call shr_timer_get(t19,'t19')
   call shr_timer_get(t20,'t20')
   call shr_timer_get(t21,'t21')
   call shr_timer_get(t22,'t22')
   call shr_timer_get(t23,'t23')
   call shr_timer_get(t24,'t24')
   call shr_timer_get(t25,'t25')
   call shr_timer_get(t26,'t26')
   call shr_timer_get(t27,'t27')
   call shr_timer_get(t28,'t28')
   call shr_timer_get(tm1,"map Sa2o")
   call shr_timer_get(tm2,"map Fa2o")
   call shr_timer_get(tm3,"map So2a")
   call shr_timer_get(tm4,"map Fo2a")
   call shr_timer_get(tm5,"map Xr2o")
   call shr_timer_start(ti1)

   !--- echo code/CVS ID to stdout --
   call date_and_time(dstr,tstr)
   write(6,F90) 
   write(6,F00) 'CCSM Coupler, version 6 (cpl6) '
   write(6,F00) 'CVS tag $Name: ccsm3_0_rel04 $'
   write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)
   write(6,F90) 

   !----------------------------------------------------------------------------
   ! determine cpl6's MPI communicator group, redirect stdin/out as necessary
   !----------------------------------------------------------------------------
   ! steder: call cpl_comm_init instead of cpl_interface_init
   !call cpl_interface_init(cpl_fields_cplname,local_comm)
   write(6,F00) 'Call to cpl_comm_init_wo_mph'
   call cpl_comm_init_wo_mph( cpl_fields_cplname, local_comm )
   !write(6,F00) 'Call to cpl_comm_init with MPH'
   !call cpl_comm_init( cpl_fields_cplname, local_comm )
!!!! shr_msg_stdio('cpl') ! consider MPH's needs wrt location of input files
   call shr_msg_chdir('cpl')
   call shr_msg_chStdIn ('cpl') ! open units 5 & 6 to named files
   !write(6,F00) 'After call to cpl_comm_init'
   write(6,F00) 'After call to cpl_comm_init_wo_mph'
   !if (cpl_comm_comp_pid == 0) call shr_msg_dirio('cpl')

   !----------------------------------------------------------------------------
   ! set flags for vector-friendly mapping
   !----------------------------------------------------------------------------
#ifdef CPP_VECTOR
   a2ovector  = .true.
   oi2avector = .true.
   r2ovector  = .true.
#else
   a2ovector  = .false.
   oi2avector = .false.
   r2ovector  = .false.
#endif

   !----------------------------------------------------------------------------
   write(6,F92) 'read input namelist file'
   !----------------------------------------------------------------------------
   call cpl_control_readNList()
   write(6,*)"debug value",cpl_control_infoDBug

!  !----------------------------------------------------------------------------
!  write(6,F92) 'create debug netCDF file'  ???
!  !----------------------------------------------------------------------------
!  call cpl_iocdf_create("DBUG.nc",desc="Add-hoc debug data file")
!  call cpl_iocdf_open  ("DBUG.nc",fid)
!  call cpl_iocdf_append(fid,date,bundle,"DBUG.nc")
!  call cpl_iocdf_close ("DBUG.nc",fid)

   !----------------------------------------------------------------------------
   write(6,F92) 'get simulation start date'
   !----------------------------------------------------------------------------
   call restart_readDate(cDate)
   write(6,F03) ' simulation start date is ',cDate

   !----------------------------------------------------------------------------
   write(6,F92) ' contract init: establishes domains & routers (excluding lnd)'
   !----------------------------------------------------------------------------

   !--- initialize info-buffers ---
   call cpl_infobuf_init(con_Xa2c%infobuf)
   call cpl_infobuf_init(con_Xi2c%infobuf)
   call cpl_infobuf_init(con_Xl2c%infobuf)
   call cpl_infobuf_init(con_Xo2c%infobuf)
   call cpl_infobuf_init(con_Xc2a%infobuf)
   call cpl_infobuf_init(con_Xc2i%infobuf)
   call cpl_infobuf_init(con_Xc2l%infobuf)
   call cpl_infobuf_init(con_Xc2o%infobuf)
   call cpl_infobuf_init(con_Dc2l%infobuf)

   !--- atm ---
   call cpl_interface_contractInit(con_Xa2c,cpl_fields_cplname, &
     cpl_fields_atmname,cpl_fields_a2c_fields,bunname="Xa2c_a", &
     decomp=cpl_control_decomp_a)
   call cpl_interface_contractInit(con_Xc2a,cpl_fields_cplname, &
     cpl_fields_atmname,cpl_fields_c2a_fields,bunname="Xc2a_a", &
     decomp=cpl_control_decomp_a)
   dom_a = con_Xa2c%domain
   con_Xc2a%domain = dom_a

   !--- ice ---
   call cpl_interface_contractInit(con_Xi2c,cpl_fields_cplname, &
     cpl_fields_icename,cpl_fields_i2c_fields,bunname="Xi2c_i", &
     decomp=cpl_control_decomp_i)
   call cpl_interface_contractInit(con_Xc2i,cpl_fields_cplname, &
     cpl_fields_icename,cpl_fields_c2i_fields,bunname="Xc2i_i", &
     decomp=cpl_control_decomp_i)
   dom_i = con_Xi2c%domain
   con_Xc2i%domain = dom_i

   !--- ocn ---
   call cpl_interface_contractInit(con_Xo2c,cpl_fields_cplname, &
     cpl_fields_ocnname,cpl_fields_o2c_fields,bunname="Xo2c_o", &
     decomp=cpl_control_decomp_o)
   call cpl_interface_contractInit(con_Xc2o,cpl_fields_cplname, &
     cpl_fields_ocnname,cpl_fields_c2o_fields,bunname="Xc2o_o", &
     decomp=cpl_control_decomp_o)
   dom_o = con_Xo2c%domain
   con_Xc2o%domain = dom_o

   !--- special domain data contract -- unique to lnd ---
   call cpl_interface_contractInit(con_Dc2l,cpl_fields_cplname, &
     cpl_fields_lndname,cpl_fields_c2lg_fields,bunname="Dc2l_l",&
     decomp=cpl_control_decomp_l)
   dom_l = con_Dc2l%domain
   if (con_Dc2l%infobuf%ibuf(cpl_fields_ibuf_inimask) == 0) then
      cpl_control_sendLndDom = .false.
   else
      cpl_control_sendLndDom = .true.
   endif

   !----------------------------------------------------------------------------
   write(6,F92) ' partial map data init, frac init '
   !----------------------------------------------------------------------------
   call cpl_map_init(map_Fo2a,dom_o,dom_a,"map_Fo2a",cpl_control_mapFn_o2aF,"dst")
   call frac_init(map_Fo2a,dom_a,dom_i,dom_l,dom_o)

   !----------------------------------------------------------------------------
   write(6,F92) ' send domain info to land model? (optional)'
   !----------------------------------------------------------------------------
   if (.not. cpl_control_sendLndDom ) then
      write(6,F00) ' lnd DOES NOT request optional domain data exchange'
   else
      write(6,F00) ' lnd requests optional domain data exchange'
      call shr_sys_flush(6)

      !--- initialize gData0,gData1,gData2 ---
      call cpl_bundle_zero(con_Dc2l%bundle)

      call cpl_mct_aVect_gather (con_Dc2l%bundle%data,gData0, &
           con_Dc2l%bundle%dom%gsMap,0,cpl_comm_comp,rcode)

      call cpl_mct_aVect_gather (dom_a%lGrid,gData1, &
                          dom_a%gsMap,0,cpl_comm_comp,rcode)

      call cpl_mct_aVect_gather (bun_frac_a%data,gData2, &
                          bun_frac_a%dom%gsMap,0,cpl_comm_comp,rcode)

      !--- set gData0 to be scattered as new data in con_Dc2l ---
      if (cpl_comm_comp_pid == 0) then
         nfld = cpl_mct_aVect_indexRA(gData2,"lfrac",perrWith='gData2 lfrac')
         lsize = cpl_mct_aVect_lsize(gData0)
         do n=1,lsize
            gData0%rAttr(cpl_fields_c2lg_alon,n)  = &
                   gData1%rAttr(cpl_fields_grid_lon,n)
            gData0%rAttr(cpl_fields_c2lg_alat,n)  = &
                   gData1%rAttr(cpl_fields_grid_lat,n)
            gData0%rAttr(cpl_fields_c2lg_aarea,n) = &
                   gData1%rAttr(cpl_fields_grid_area,n)
            gData0%rAttr(cpl_fields_c2lg_amask,n) = &
                   gData1%rAttr(cpl_fields_grid_mask,n)
            if (abs(gData2%rAttr(nfld,n)) < 1.0e-06) then
               gData0%rAttr(cpl_fields_c2lg_lmask,n) = 0.0
               gData0%rAttr(cpl_fields_c2lg_lfrac,n) = 0.0
               gData1%rAttr(cpl_fields_grid_mask,n)  = 0.0
            else
               gData0%rAttr(cpl_fields_c2lg_lmask,n) = 1.0
               gData0%rAttr(cpl_fields_c2lg_lfrac,n) = gData2%rAttr(nfld,n)
               gData1%rAttr(cpl_fields_grid_mask,n)  = 1.0
            endif
         enddo
      endif

      !--- reset dom_l based on dom_a ---
      !--- cpl_fields_grid_mask is not from dom_a (see above loop) ---
      call cpl_mct_aVect_scatter(gData1,dom_l%lGrid, &
        dom_l%gsMap,0,cpl_comm_comp,rcode)

      !--- scatter gData0 to con_Dc2l bundle ---
      call cpl_mct_aVect_scatter(gData0,con_Dc2l%bundle%data, &
        con_Dc2l%bundle%dom%gsMap,0,cpl_comm_comp,rcode)

      !--- clean up ---
      if (cpl_comm_comp_pid == 0) then
         call cpl_mct_aVect_clean(gData0)
         call cpl_mct_aVect_clean(gData1)
         call cpl_mct_aVect_clean(gData2)
      endif

      !--- send con_Dc2l ---
      call cpl_interface_contractSend(cpl_fields_lndname,con_Dc2l)
   endif

   !----------------------------------------------------------------------------
   write(6,F92) ' contract init: establish domain & router for lnd'
   !----------------------------------------------------------------------------
   call cpl_interface_contractInit(con_Xl2c,cpl_fields_cplname, &
     cpl_fields_lndname,cpl_fields_l2c_fields,bunname="Xl2c_l", &
     decomp=cpl_control_decomp_l)
   call cpl_interface_contractInit(con_Xc2l,cpl_fields_cplname, &
     cpl_fields_lndname,cpl_fields_c2l_fields,bunname="Xc2l_l", &
     decomp=cpl_control_decomp_l)
   dom_l%gsMap = con_Xl2c%domain%gsMap
   con_Xc2l%domain = dom_l
   con_Xl2c%domain = dom_l

   call cpl_interface_contractInit(con_Xr2c,cpl_fields_cplname, &
     cpl_fields_rtmname,cpl_fields_r2c_fields,bunname="Xr2c_r", &
     decomp=cpl_control_decomp_r)
   dom_r = con_Xr2c%bundle%dom

   !----------------------------------------------------------------------------
   write(6,F92) ' check atm/lnd and ocn/ice model domains for consistency'
   !----------------------------------------------------------------------------
   call cpl_domain_compare(con_Xa2c%domain,con_Xl2c%domain, &
   &                       enforce_grid=.true.,enforce_area=.true.)
   call cpl_domain_compare(con_Xo2c%domain,con_Xi2c%domain,enforce_all =.true.)

   !----------------------------------------------------------------------------
   write(6,F92) ' init data_mod bundles, tavg history bundles, maps'
   !----------------------------------------------------------------------------
   call data_bundleInit()
   call data_mapInit()
   call history_avbundleInit()

   !*** KLUDGE ******************************************************************
   ! fix_So2c_a is similar to bun_So2c_a but with a slightly different mapping
   ! fix_frac_a is used to re-normalize map_Fo2a mapping for fix_So2c_a only
   !*****************************************************************************
   call cpl_bundle_initv(fix_So2c_a,"fix_So2c_a",bun_So2c_a,bun_So2c_a%dom)
   call cpl_bundle_initv(fix_frac_a,"fix_frac_a",bun_frac_a,bun_frac_a%dom)
   !*****************************************************************************

   !----------------------------------------------------------------------------
   write(6,F92) 'msg send-init (send start date)'
   !----------------------------------------------------------------------------

   !--- set integer values ---
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_rcode)   = 0
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_stopeod) = 0
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_stopnow) = 0
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_cdate)   = cDate
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_sec)     = 0
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_infotim) = 0
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_infobug) = cpl_control_infodbug

   !--- set real values ---
   con_Xc2a%infobuf%rbuf(cpl_fields_rbuf_spval)   = cpl_const_spval
   con_Xc2a%infobuf%rbuf(cpl_fields_rbuf_eccen)   = cpl_control_orbEccen 
   con_Xc2a%infobuf%rbuf(cpl_fields_rbuf_obliqr)  = cpl_control_orbObliqr
   con_Xc2a%infobuf%rbuf(cpl_fields_rbuf_lambm0)  = cpl_control_orbLambm0
   con_Xc2a%infobuf%rbuf(cpl_fields_rbuf_mvelpp)  = cpl_control_orbMvelpp

   !--- all components get the same info-buffer ---
   con_Xc2i%infobuf = con_Xc2a%infobuf
   con_Xc2o%infobuf = con_Xc2a%infobuf
   con_Xc2l%infobuf = con_Xc2a%infobuf

   call cpl_interface_infobufSend(cpl_fields_atmname,con_Xc2a%infobuf%ibuf,con_Xc2a%infobuf%rbuf)
   call cpl_interface_infobufSend(cpl_fields_icename,con_Xc2i%infobuf%ibuf,con_Xc2i%infobuf%rbuf)
   call cpl_interface_infobufSend(cpl_fields_ocnname,con_Xc2o%infobuf%ibuf,con_Xc2o%infobuf%rbuf)
   call cpl_interface_infobufSend(cpl_fields_lndname,con_Xc2l%infobuf%ibuf,con_Xc2l%infobuf%rbuf)

   !----------------------------------------------------------------------------
   write(6,F92) 'verify acceptable coupling intervals, check for dead components'
   !----------------------------------------------------------------------------
   cpl_control_nCpl_a = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ncpl)
   cpl_control_nCpl_i = con_Xi2c%infobuf%ibuf(cpl_fields_ibuf_ncpl)
   cpl_control_nCpl_l = con_Xl2c%infobuf%ibuf(cpl_fields_ibuf_ncpl)
   cpl_control_nCpl_o = con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_ncpl)
   cpl_control_nCpl_r = con_Xr2c%infobuf%ibuf(cpl_fields_ibuf_ncpl)
   ncpl_a = cpl_control_nCpl_a
   ncpl_i = cpl_control_nCpl_i
   ncpl_l = cpl_control_nCpl_l
   ncpl_o = cpl_control_nCpl_o
   ncpl_r = cpl_control_nCpl_r
   if ( (ncpl_a < ncpl_o) .or. (mod(ncpl_a,ncpl_o)/=0)) then
      write(6,*) 'ERROR: unacceptable ncpl_a,ncpl_o=',ncpl_a,ncpl_o
      call shr_sys_flush(6)
      call shr_sys_abort()
   else if ( (ncpl_a /= ncpl_i) .or. ncpl_a /= ncpl_l ) then
      write(6,*) 'ERROR: unacceptable ncpl_[ail] = ',ncpl_a,ncpl_i,ncpl_l
      call shr_sys_flush(6)
      call shr_sys_abort()
   else
      write(6,F01) 'ncpl_[ailro] = ',ncpl_a,ncpl_i,ncpl_l,ncpl_r,ncpl_o
      call shr_sys_flush(6)
   end if

   !--- check for dead models ---
   cpl_control_dead_a  = .false.
   cpl_control_dead_i  = .false.
   cpl_control_dead_l  = .false.
   cpl_control_dead_o  = .false.
   cpl_control_dead_ao = .false.
   if (con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_dead) /= 0) then
      write(6,F00) 'dead atm component'
      cpl_control_dead_a  = .true.
      cpl_control_dead_ao = .true.
   end if
   if (con_Xi2c%infobuf%ibuf(cpl_fields_ibuf_dead) /= 0) then
      write(6,F00) 'dead ice component'
      cpl_control_dead_i  = .true.
   end if
   if (con_Xl2c%infobuf%ibuf(cpl_fields_ibuf_dead) /= 0) then
      write(6,F00) 'dead lnd component'
      cpl_control_dead_l  = .true.
   end if
   if (con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_dead) /= 0) then
      write(6,F00) 'dead ocn component'
      cpl_control_dead_o  = .true.
      cpl_control_dead_ao = .true.
   end if
   if ( .not. (cpl_control_dead_a .or. cpl_control_dead_i .or.  &
   &           cpl_control_dead_l .or. cpl_control_dead_o )     ) then
      write(6,F00) 'no dead components'
   end if

   !----------------------------------------------------------------------------
   write(6,F92) 'initialize model start date and control flags'
   ! now that we know what the start date is and what the time step is
   !----------------------------------------------------------------------------
   date = shr_date_initCDate(cDate,ncpl_a)
   call cpl_control_init(date)
   call cpl_control_update(date)

   !--- reset bundle domain areas equal to mapping domain areas ---

   nfld  = cpl_mct_aVect_indexRA(dom_a%lGrid     ,"aream"             , &
                                 perrWith='cpl main.F90')
   nfld1 = cpl_mct_aVect_indexRA(map_Fo2a%areasrc,cpl_map_areaAV_field, &
                                 perrWith='cpl main.F90')
   dom_a%lGrid%rAttr(nfld,:) = map_Fo2a%areadst%rAttr(nfld1,:)
   dom_l%lGrid%rAttr(nfld,:) = map_Fo2a%areadst%rAttr(nfld1,:)
   dom_o%lGrid%rAttr(nfld,:) = map_Fo2a%areasrc%rAttr(nfld1,:)
   dom_i%lGrid%rAttr(nfld,:) = map_Fo2a%areasrc%rAttr(nfld1,:)
   dom_r%lGrid%rAttr(nfld,:) = map_Xr2o%areasrc%rAttr(nfld1,:)
   con_Xc2a%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areadst%rAttr(nfld1,:)
   con_Xc2l%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areadst%rAttr(nfld1,:)
   con_Dc2l%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areadst%rAttr(nfld1,:)
   con_Xc2o%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areasrc%rAttr(nfld1,:)
   con_Xc2i%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areasrc%rAttr(nfld1,:)
   con_Xa2c%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areadst%rAttr(nfld1,:)
   con_Xl2c%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areadst%rAttr(nfld1,:)
   con_Xo2c%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areasrc%rAttr(nfld1,:)
   con_Xi2c%domain%lGrid%rAttr(nfld,:) = map_Fo2a%areasrc%rAttr(nfld1,:)
   con_Xr2c%domain%lGrid%rAttr(nfld,:) = map_Xr2o%areasrc%rAttr(nfld1,:)

   !--- initialize interface area correction factors ---
   call areafact_init(dom_a,dom_i,dom_l,dom_o,dom_r)

   !----------------------------------------------------------------------------
   write(6,F92) 'recv IC data from models'
   !----------------------------------------------------------------------------
   call cpl_interface_contractRecv(cpl_fields_atmname,con_Xa2c)
   call cpl_bundle_mult(con_Xa2c%bundle,bun_areafact_a,'comp2cpl',  &
                        bunlist=cpl_fields_a2c_fluxes)
   cpl_control_icData_a = .false.
   if (con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_userest) /= 0) cpl_control_icData_a = .true.

   call cpl_interface_contractRecv(cpl_fields_icename,con_Xi2c)
   call cpl_bundle_mult(con_Xi2c%bundle,bun_areafact_i,'comp2cpl',  &
                        bunlist=cpl_fields_i2c_fluxes)
   call flux_albi(date,con_Xi2c%bundle) ! cpl modifies this ice output state
   cpl_control_icData_i = .false.
   if (con_Xi2c%infobuf%ibuf(cpl_fields_ibuf_userest) /= 0) cpl_control_icData_i = .true.

   call cpl_interface_contractRecv(cpl_fields_lndname,con_Xl2c)
   call cpl_bundle_mult(con_Xl2c%bundle,bun_areafact_l,'comp2cpl',  &
                        bunlist=cpl_fields_l2c_fluxes)
   cpl_control_icData_l = .false.
   if (con_Xl2c%infobuf%ibuf(cpl_fields_ibuf_userest) /= 0) cpl_control_icData_l = .true.
   call cpl_interface_contractRecv(cpl_fields_lndname,con_Xr2c)
   call cpl_bundle_mult(con_Xr2c%bundle,bun_areafact_r,'comp2cpl',  &
                        bunlist=cpl_fields_r2c_fluxes)
!  --- use the land model infobuf flag, could use runoff infobuf flag ? ---
   cpl_control_icData_r =  cpl_control_icData_l 

   call cpl_interface_contractRecv(cpl_fields_ocnname,con_Xo2c)
   call cpl_bundle_mult(con_Xo2c%bundle,bun_areafact_o,'comp2cpl',  &
                        bunlist=cpl_fields_o2c_fluxes)
   call flux_albo(date,bun_oalbedo_o) ! cpl computes this ocn "output" state
   cpl_control_icData_o = .false.
   if (con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_userest) /= 0) cpl_control_icData_o = .true.

   !----------------------------------------------------------------------------
   write(6,F92) 'read IC data from restart file?'
   !----------------------------------------------------------------------------
   call restart_read(date)

   !----------------------------------------------------------------------------
   write(6,F92) 'process IC data'
   !----------------------------------------------------------------------------

   !--- process atm IC data ---
   call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,bun_Fa2c_a)
   call cpl_bundle_zero (bun_precip_a)
   call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainc')
   call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainl')
   call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowc')
   call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowl')
   cpl_control_fluxAShift = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ashift) 

   !--- process ice IC data ---
   lsize = cpl_mct_aVect_lsize(con_Xi2c%bundle%data)
   allocate(ifrac_i(lsize))
   call cpl_mct_aVect_getRAttr(con_Xi2c%bundle%data,"Si_ifrac",ifrac_i,rcode)
   call frac_set(ifrac_i,map_Fo2a,dom_a,dom_i,dom_l,dom_o) ! init all fracs
   call cpl_bundle_split(con_Xi2c%bundle,bun_Si2c_i,bun_Fi2c_i)
   call cpl_map_bun(bun_Si2c_i,bun_Si2c_a,map_So2a, &
                    bun_frac_i,'ifrac',bun_frac_a,'ifrac')
   call cpl_map_bun(bun_Fi2c_i,bun_Fi2c_a,map_Fo2a, &
                    bun_frac_i,'ifrac',bun_frac_a,'ifrac')

   !--- process lnd IC data ---
   call cpl_bundle_split(con_Xl2c%bundle,bun_Sl2c_l,bun_Fl2c_l )
   !--- map land fields to ocean, not allowed now.
!  call cpl_map_bun(bun_Sl2c_l,bun_Sl2c_o,map_Sa2o)
!  call cpl_map_bun(bun_Fl2c_l,bun_Fl2c_o,map_Fa2o)
   !--- map land fields to ocean, not allowed now.

   !--- process ocn IC data ---
   call cpl_bundle_split(con_Xo2c%bundle,bun_So2c_o,bun_Fo2c_o)
   call cpl_map_bun(bun_So2c_o,bun_So2c_a,map_So2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac')
   !*** KLUDGE - start *********************************************************
   call cpl_map_bun(bun_frac_o,fix_frac_a,map_So2a)               ! wrt KLUDGE
   call cpl_map_bun(bun_So2c_o,fix_So2c_a,map_So2a, &             ! wrt KLUDGE
                    bun_frac_o,'afrac',fix_frac_a,'afrac')        ! wrt KLUDGE
   !*** KLUDGE - end ***********************************************************
   call cpl_map_bun(bun_Fo2c_o,bun_Fo2c_a,map_Fo2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac')
   call cpl_map_bun(bun_aoflux_o,bun_aoflux_a,map_Fo2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac')
   call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac')

   !----------------------------------------------------------------------------
   write(6,F92) 'optional atm initialization: send albedos & recv new solar? '
   !----------------------------------------------------------------------------
   if ( con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_xalbic) == 0) then
      cpl_control_sendAtmAlb = .false.
   else
      cpl_control_sendAtmAlb = .true.
   end if
   if (.not. cpl_control_sendAtmAlb ) then
      write(6,F00) '* atm component requests NO recalculation of initial solar'
   else
      write(6,F00) '* atm component requests recalculation of initial solar'

      !--- map ocn & ice albedos ondo atm domain ---
      ! already done above

      !--- merge atm inputs ---
      !--- only albedo, surface temp, snow, ifrac & ofrac need to be valid ---
!     call merge_atm()
      call merge_atm(fix_So2c_a)  ! KLUDGE

      !--- send message ---
      write(6,F00) 'send albedos to atm, recv new atm IC''s'
      call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'cpl2comp',  &
                           bunlist=cpl_fields_c2a_fluxes)
      call cpl_interface_contractSend(cpl_fields_atmname,con_Xc2a)
      call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'comp2cpl',  &
                           bunlist=cpl_fields_c2a_fluxes)

      !--- wait for new atm IC data ---
      call cpl_interface_contractRecv(cpl_fields_atmname,con_Xa2c)
      call cpl_bundle_mult(con_Xa2c%bundle,bun_areafact_a,'comp2cpl',  &
                           bunlist=cpl_fields_a2c_fluxes)

      !--- process atm IC data ---
      call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,bun_Fa2c_a)
      call cpl_bundle_zero (bun_precip_a)
      call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainc')
      call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainl')
      call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowc')
      call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowl')
      cpl_control_fluxAShift = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ashift)
   end if

   !----------------------------------------------------------------------------
   write(6,F92) 'create data as necessary for 1st iteration of main event loop'
   !----------------------------------------------------------------------------

   !--- map ocn & ice albedos onto atm domain ---
   call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
                    bun_frac_o,'afrac',bun_frac_a,'ofrac')
   call cpl_map_bun(con_Xr2c%bundle,bun_Xr2c_o,map_Xr2o)

   !--- prepare cpl->atm bundle (as necessary for 1st diag calc) ---
   !--- merge atm inputs (necessary for 1st diag calc) ---
!  call merge_atm()
   call merge_atm(fix_So2c_a)  ! KLUDGE

   !----------------------------------------------------------------------------
   write(6,F92) 'start of main integration loop' 
   !----------------------------------------------------------------------------

   call shr_date_getYMD(date,year,month,day,sec)
   call tStamp_write("cpl",year,month,day,sec)

   call shr_timer_stop(ti1)
   call shr_timer_zero_all()
   call shr_timer_start(t00)

   do while ( .not. cpl_control_stopNow )
      do n=1,ncpl_a

         call shr_timer_start(t01)
 
         !--- send msg to ocn ---
         if(mod(n-1,ncpl_a/ncpl_o) == 0 ) then
            if (.not. cpl_control_lagOcn) then
              call cpl_bundle_mult(con_Xc2o%bundle,bun_areafact_o,'cpl2comp',  &
                                   bunlist=cpl_fields_c2o_fluxes)
              call cpl_interface_contractSend(cpl_fields_ocnname,con_Xc2o)
              call cpl_bundle_mult(con_Xc2o%bundle,bun_areafact_o,'comp2cpl',  &
                                   bunlist=cpl_fields_c2o_fluxes)
            endif
            call cpl_bundle_zero(bun_Xc2oPSUM_o) ! zero-out partial sum
         endif

         call shr_timer_stop (t01) ; call shr_timer_start(t02)

         !--- send msg to lnd ---
         if (mod(n-1,ncpl_a/ncpl_l) == 0 ) then
            call cpl_bundle_gather(con_Xc2l%bundle, bun_Sa2c_a, bun_Fa2c_a, &
            &                                       bun_Sl2c_l, bun_Fl2c_l, &
            &                                       bun_So2c_a, bun_Fo2c_a, &
            &                                       bun_Si2c_a, bun_Fi2c_a  ) 
            call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,'cpl2comp', &
                                 bunlist=cpl_fields_c2l_fluxes)
            call cpl_interface_contractSend(cpl_fields_lndname,con_Xc2l)
            call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,'comp2cpl', &
                                 bunlist=cpl_fields_c2l_fluxes)
         endif

         call shr_timer_stop (t02) ; call shr_timer_start(t03)

         !--- compute, map, merge ice inputs ----
         call shr_timer_start(tm1)
         call cpl_map_bun(bun_Sa2c_a,bun_Sa2c_o,map_Sa2o,mvector=a2ovector)
         call shr_timer_stop(tm1) ; call shr_timer_start(tm2)
         call cpl_map_bun(bun_Fa2c_a,bun_Fa2c_o,map_Fa2o,mvector=a2ovector)
         call cpl_map_bun(bun_precip_a,bun_precip_o,map_Fa2o,mvector=a2ovector)
         call shr_timer_stop(tm2) ; 

         call shr_timer_stop (t03) ; call shr_timer_start(t04)

         !--- force zero net water flux into ocn+ice ? ----
         if (cpl_control_fluxEpbal(1:3) /= 'off') then
             cpl_control_fluxEPfac = con_Xo2c%infobuf%ibuf(cpl_fields_ibuf_precAdj)
             cpl_control_fluxEPfac = cpl_control_fluxEPfac * 1.0e-6_r8
             
             call flux_epbal(date,bun_aoflux_o,con_Xi2c%bundle, &
                  &                    bun_precip_o,bun_Xr2c_o,bun_frac_o)
         endif

         call shr_timer_stop (t04) ; call shr_timer_start(t05)

         !--- merge total snow and precip for ice input ---
         call cpl_bundle_zero (con_Xc2i%bundle,'Faxc_rain')
         call cpl_bundle_zero (con_Xc2i%bundle,'Faxc_snow')
         call cpl_bundle_copy(bun_precip_o,bunrList='Faxc_rain',&
            bunTrList='Faxc_rain',outbun=con_Xc2i%bundle)
         call cpl_bundle_copy(bun_precip_o,bunrList='Faxc_snow',&
            bunTrList='Faxc_snow',outbun=con_Xc2i%bundle)

         call shr_timer_stop (t05) ; call shr_timer_start(t06)

         !--- correct a->o vector mapping near NP ---
         call cpl_map_npfix(bun_Sa2c_a,bun_Sa2c_o,'Sa_u','Sa_v')

         call shr_timer_stop (t06) ; call shr_timer_start(t07)

         !--- send msg to ice ---
         if (mod(n-1,ncpl_a/ncpl_i) == 0 ) then
            call cpl_bundle_gather(con_Xc2i%bundle, bun_Sa2c_o, bun_Fa2c_o, &
!--- not allowed now ------------------------------ bun_Sl2c_o, bun_Fl2c_o, &
            &                                       bun_So2c_o, bun_Fo2c_o, &
            &                                       bun_Si2c_i, bun_Fi2c_i, &
            &                                       fcopy = .true.  )
            call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,'cpl2comp', &
                                 bunlist=cpl_fields_c2i_fluxes)
            call cpl_interface_contractSend(cpl_fields_icename,con_Xc2i)
            call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,'comp2cpl', &
                                 bunlist=cpl_fields_c2i_fluxes)
         endif

         call shr_timer_stop (t07) ; call shr_timer_start(t08)

         !--- compute net solar flux into ocn ---
         call flux_solar(bun_Fa2c_o, bun_oalbedo_o, bun_aoflux_o )

         call shr_timer_stop (t08) ; call shr_timer_start(t09)

         !--- merge ocn inputs ---
         call merge_ocn()

         call shr_timer_stop (t09) ; call shr_timer_start(t10)

         !--- form partial sum of tavg ocn inputs (virtual "send" to ocn) ---
         call cpl_bundle_accum(bun_Xc2oSNAP_o,outbun=bun_Xc2oPSUM_o)

         call shr_timer_stop (t10) ; call shr_timer_start(t11)

         !--- write bit-for-bit check info ---
         if ( cpl_control_infoBCheck ) then
            call bitCheck_write(date,cpl_control_caseName,con_Xa2c%bundle,"Sa_tbot")
            call bitCheck_write(date,cpl_control_caseName,con_Xi2c%bundle,"Si_ifrac")
            call bitCheck_write(date,cpl_control_caseName,con_Xl2c%bundle,"Sl_t")
            call bitCheck_write(date,cpl_control_caseName,con_Xo2c%bundle,"So_t")
         end if

         !--- do diagnostics ---
         if ((.not.cpl_control_lagOcn) .or. n /= 1)  then
            call diag_dodiag(date,                                             &
            & con_Xa2c%bundle,con_Xc2a%bundle,con_Xl2c%bundle,con_Xc2l%bundle, &
            & con_Xr2c%bundle,con_Xi2c%bundle,con_Xc2i%bundle,con_Xo2c%bundle, &
            & bun_Xc2oSNAP_o ,bun_Fa2c_o     ,bun_oalbedo_o,                   &
            & bun_frac_l     ,bun_frac_i     ,bun_frac_o   )
         end if

         call shr_timer_stop (t11) ; call shr_timer_start(t12)

         !--- compute ocn albedos (virtual "recv" from ocn) ---
          call flux_albo(date,bun_oalbedo_o)

         call shr_timer_stop (t12) ; call shr_timer_start(t13)

         !--- compute atm/ocn fluxes ---
          call flux_atmOcn(con_Xo2c%bundle,bun_Sa2c_o,cpl_control_dead_ao,bun_aoflux_o )

         call shr_timer_stop (t13) ; call shr_timer_start(t14)

         !--- recv msg from ice ---
         if (mod(n-1,ncpl_a/ncpl_i) == 0 ) then
            call cpl_interface_contractRecv(cpl_fields_icename,con_Xi2c)
            call cpl_bundle_mult(con_Xi2c%bundle,bun_areafact_i,'comp2cpl',  &
              bunlist=cpl_fields_i2c_fluxes)

            !--- update surface fracs wrt new ice frac ---
            call cpl_mct_aVect_getRAttr(con_Xi2c%bundle%data,"Si_ifrac",ifrac_i,rcode)
             call frac_set(ifrac_i,map_Fo2a,dom_a,dom_i,dom_l,dom_o)
         endif

         call shr_timer_stop (t14) ; call shr_timer_start(t15)

         !--- add diurnal cycle to ice albedos (?) ---
          call flux_albi(date,con_Xi2c%bundle)
          call cpl_bundle_split(con_Xi2c%bundle,bun_Si2c_i,bun_Fi2c_i)

         call shr_timer_stop (t15) ; call shr_timer_start(t16)

         !--- map ocn states/fluxes to atm ---
         call shr_timer_start(tm3)
          call cpl_map_bun(bun_Si2c_i,bun_Si2c_a,map_So2a, &
                           bun_frac_i,'ifrac',bun_frac_a,'ifrac',oi2avector)
          call cpl_map_bun(bun_So2c_o,bun_So2c_a,map_So2a,  &
                           bun_frac_o,'afrac',bun_frac_a,'ofrac',oi2avector)
   !*** KLUDGE - start *********************************************************
    call cpl_map_bun(bun_frac_o,fix_frac_a,map_So2a,mvector=oi2avector)  ! wrt KLUDGE
    call cpl_map_bun(bun_So2c_o,fix_So2c_a,map_So2a, &                  ! wrt KLUDGE
                     bun_frac_o,'afrac',fix_frac_a,'afrac',oi2avector)  ! wrt KLUDGE
   !*** KLUDGE - end ***********************************************************
          call shr_timer_stop(tm3) ; call shr_timer_start(tm4)
          call cpl_map_bun(bun_Fi2c_i,bun_Fi2c_a,map_Fo2a, &
                           bun_frac_i,'ifrac',bun_frac_a,'ifrac',oi2avector)
          call cpl_map_bun(bun_Fo2c_o,bun_Fo2c_a,map_Fo2a, &
                           bun_frac_o,'afrac',bun_frac_a,'ofrac',oi2avector)
          call cpl_map_bun(bun_aoflux_o,bun_aoflux_a,map_Fo2a, &
                           bun_frac_o,'afrac',bun_frac_a,'ofrac',oi2avector)
          call shr_timer_stop(tm4) ; call shr_timer_start(tm3)
          call cpl_map_bun(bun_oalbedo_o,bun_oalbedo_a,map_So2a, &
                           bun_frac_o,'afrac',bun_frac_a,'ofrac',oi2avector)
         call shr_timer_stop(tm3)

         call shr_timer_stop (t16) ; call shr_timer_start(t17)

         !--- recv msg from lnd ---
         if (mod(n-1,ncpl_a/ncpl_l) == 0 ) then
            call cpl_interface_contractRecv(cpl_fields_lndname,con_Xl2c)
            call cpl_bundle_mult(con_Xl2c%bundle,bun_areafact_l,'comp2cpl',  &
                                 bunlist=cpl_fields_l2c_fluxes)
            call cpl_bundle_split(con_Xl2c%bundle,bun_Sl2c_l,bun_Fl2c_l )
         endif

         call shr_timer_stop (t17) ; call shr_timer_start(t18)

         if (mod(n,ncpl_a/ncpl_r) == 0 ) then
            call cpl_interface_contractRecv(cpl_fields_lndname,con_Xr2c)
            call cpl_bundle_mult(con_Xr2c%bundle,bun_areafact_r,'comp2cpl',  &
                                 bunlist=cpl_fields_r2c_fluxes)
         endif

         call shr_timer_stop (t18) ; call shr_timer_start(t19)

         !--- diagnostics: verify net solar calcs are coordinated ---
         if(cpl_control_diagNow .and. .not.cpl_control_lagOcn)  &
          call diag_solar (con_Xa2c%bundle,con_Xl2c%bundle,con_Xi2c%bundle,bun_frac_l,bun_frac_i)

         call shr_timer_stop (t19) ; call shr_timer_start(t20)

         !--- merge atm states & fluxes ---
!        call merge_atm()
          call merge_atm(fix_So2c_a)  ! KLUDGE
         
         call shr_timer_stop (t20) ; call shr_timer_start(t21)

         !--- send msg to atm ---
         if (mod(n-1,ncpl_a/ncpl_a) == 0 ) then
            call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'cpl2comp',  &
                                 bunlist=cpl_fields_c2a_fluxes)
            call cpl_interface_contractSend(cpl_fields_atmname,con_Xc2a)
            call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'comp2cpl',  &
                                 bunlist=cpl_fields_c2a_fluxes)
         endif

         call shr_timer_stop (t21) ; call shr_timer_start(t22)

         !--- map land to ocn, not allowed now ---
!        call shr_timer_start(tm1)
!        call cpl_map_bun(bun_Sl2c_l,bun_Sl2c_o,map_Sa2o)
!        call shr_timer_stop(tm1) ; call shr_timer_start(tm2)
!        call cpl_map_bun(bun_Fl2c_l,bun_Fl2c_o,map_Fa2o)
!        call shr_timer_stop(tm2)
         !--- map land to ocn, not allowed now ---

         call shr_timer_start(tm5)
          call cpl_map_bun(con_Xr2c%bundle,bun_Xr2c_o,map_Xr2o,mvector=r2ovector)
         call shr_timer_stop(tm5)

         call shr_timer_stop (t22) ; call shr_timer_start(t23)

         !--- create history files ---
         call history_write(date)
         call history_avwrite(date)

         call shr_timer_stop (t23) ; call shr_timer_start(t24)

         !--- recv msg from ocn ---
         if (mod(n,ncpl_a/ncpl_o) == 0 ) then
            !--- form tavg of ocn inputs ---
            call cpl_bundle_avg (bun_Xc2oPSUM_o)
            call cpl_bundle_copy(bun_Xc2oPSUM_o,outbun=con_Xc2o%bundle)

            
            if (.not. cpl_control_lagOcn) then
               call cpl_interface_contractRecv(cpl_fields_ocnname,con_Xo2c)
               call cpl_bundle_mult(con_Xo2c%bundle,bun_areafact_o,'comp2cpl',  &
                                    bunlist=cpl_fields_o2c_fluxes)
            endif
            call cpl_bundle_split(con_Xo2c%bundle,bun_So2c_o,bun_Fo2c_o)

            !--- start normal interaction with ocn ---
            if (cpl_control_lagOcn) then
               write(6,*)
               write(6,*) ' Start of time coordinated integration'
               write(6,*)
               cpl_control_lagOcn = .false.
            end if
         endif

         call shr_timer_stop (t24) ; call shr_timer_start(t25)

         !--- recv msg from atm ---
         if (mod(n-1,ncpl_a/ncpl_a) == 0 ) then
            call cpl_interface_contractRecv(cpl_fields_atmname,con_Xa2c)
            call cpl_bundle_mult(con_Xa2c%bundle,bun_areafact_a,'comp2cpl',  &
                                 bunlist=cpl_fields_a2c_fluxes)
            call cpl_bundle_split(con_Xa2c%bundle,bun_Sa2c_a,bun_Fa2c_a)
            cpl_control_fluxAShift = con_Xa2c%infobuf%ibuf(cpl_fields_ibuf_ashift)
            !---form total rain and snow
            call cpl_bundle_zero (bun_precip_a)
            call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainc')
            call cpl_bundle_add(bun_precip_a,'Faxc_rain',bun_Fa2c_a,'Faxa_rainl')
            call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowc')
            call cpl_bundle_add(bun_precip_a,'Faxc_snow',bun_Fa2c_a,'Faxa_snowl')
         endif

         call shr_timer_stop (t25) ; call shr_timer_start(t26)

         !--- advance date & update control flags ---
         call shr_date_adv1step(date)
         call shr_date_getCDate(date,cDate,sec)
         call shr_date_getYMD(date,year,month,day,sec)
         if (cpl_control_infodbug >= 2 .or. sec == 0) then
            call tStamp_write("cpl",year,month,day,sec)
         endif
         call cpl_control_update(date)

         !--- set infobuf flags --- KLUDGE: need to move this elsewhere
         con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restEOD) = 0
         con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restNow) = 0
         if (cpl_control_restEOD) con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restEOD) = 1
         if (cpl_control_stopEOD) con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_stopEOD) = 1
         if (cpl_control_restNow) con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_restNow) = 1
         con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_cdate)   = cDate
         con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_sec)     = sec
         con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_infobug) = cpl_control_infodbug
         con_Xc2l%infobuf = con_Xc2a%infobuf
         con_Xc2i%infobuf = con_Xc2a%infobuf
         con_Xc2o%infobuf = con_Xc2a%infobuf
         call shr_sys_flush(6)

         call shr_timer_stop (t26) ; call shr_timer_start(t27)

         !--- create a restart file ---
         call restart_write(date)

         call shr_timer_stop (t27)

      enddo

      !--- verify all models are coordinated in time ---
      if (cpl_control_dead_a .or. cpl_control_dead_i .or. &
      &   cpl_control_dead_l .or. cpl_control_dead_o     ) then
         !--- dead models have erroneous date? - don't enforce coordination ---
         call timeCheck(date,.false.,.false.) ! date,print,enforce
      else
         call timeCheck(date,.false.,.true. ) ! date,print,enforce
      end if
      call shr_sys_flush(6)
   enddo
   call shr_timer_stop(t00)

   !----------------------------------------------------------------------------
   write(6,F92) 'end of main integration loop' 
   !----------------------------------------------------------------------------

   !--- last chance to create t-avg history file ---
   call history_avwrite(date)

   !--- send last message with a stopnow signal ---
   con_Xc2a%infobuf%ibuf(cpl_fields_ibuf_stopnow) = 1
   con_Xc2i%infobuf%ibuf(cpl_fields_ibuf_stopnow) = 1
   con_Xc2l%infobuf%ibuf(cpl_fields_ibuf_stopnow) = 1
   con_Xc2o%infobuf%ibuf(cpl_fields_ibuf_stopnow) = 1

   call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'cpl2comp',  &
                        bunlist=cpl_fields_c2a_fluxes)
   call cpl_interface_contractSend(cpl_fields_atmname,con_Xc2a)
   call cpl_bundle_mult(con_Xc2a%bundle,bun_areafact_a,'comp2cpl',  &
                        bunlist=cpl_fields_c2a_fluxes)

   call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,'cpl2comp',  &
                        bunlist=cpl_fields_c2i_fluxes)
   call cpl_interface_contractSend(cpl_fields_icename,con_Xc2i)
   call cpl_bundle_mult(con_Xc2i%bundle,bun_areafact_i,'comp2cpl',  &
                        bunlist=cpl_fields_c2i_fluxes)

   call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,'cpl2comp',  &
                        bunlist=cpl_fields_c2l_fluxes)
   call cpl_interface_contractSend(cpl_fields_lndname,con_Xc2l)
   call cpl_bundle_mult(con_Xc2l%bundle,bun_areafact_l,'comp2cpl',  &
                        bunlist=cpl_fields_c2l_fluxes)

   call cpl_bundle_mult(con_Xc2o%bundle,bun_areafact_o,'cpl2comp',  &
                        bunlist=cpl_fields_c2o_fluxes)
   call cpl_interface_contractSend(cpl_fields_ocnname,con_Xc2o)
   call cpl_bundle_mult(con_Xc2o%bundle,bun_areafact_o,'comp2cpl',  &
                        bunlist=cpl_fields_c2o_fluxes)

   call shr_timer_print_all
   call shr_timer_free_all

   !----------------------------------------------------------------------------
   write(6,F92) 'end of main program' 
   !----------------------------------------------------------------------------

!  call cpl_iocdf_close (fid)  ! DEBUG.nc file

   call cpl_interface_finalize(cpl_fields_cplname)

stop
end program cpl
!===============================================================================
