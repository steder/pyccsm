!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      module forcing_coupled

!***********************************************************************
!
!  CVS:$Id: forcing_coupled.F,v 1.2.2.27 2004/04/26 15:28:18 njn01 Exp $
!  CVS:$Name: ccsm3_0_rel04 $
!
!-----------------------------------------------------------------------


      use kinds_mod
      use model_size
      use domain
      use communicate
      use constants
      use io
      use stencils
      use time_management
      use grid
      use prognostic
      use exit_mod
      use ice
      use forcing_shf
      use forcing_sfwf
      use qflux_mod
      use ms_balance
      use timers
      use cpl_contract_mod
      use cpl_interface_mod
      use cpl_comm_mod
      use cpl_fields_mod
      use shr_sys_mod
      use registry
      
      implicit none
      save
      
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

      logical (kind=log_kind) ::
     &  lcoupled              ! flag for coupled forcing
     &, ldiag_cpl = .false.


      integer (kind=int_kind) ::
     &  coupled_freq_iopt,    ! coupler frequency option
     &  coupled_freq          ! frequency of coupling

!-----------------------------------------------------------------------
!
!  diurnal cycle switch for the net shortwave heat flux
!
!     qsw_diurnal_cycle = .T.  diurnal cycle is ON
!                       = .F.  diurnal cycle is OFF
!
!-----------------------------------------------------------------------

      logical (kind=log_kind) :: qsw_diurnal_cycle

      real (kind=dbl_kind), dimension(:), allocatable ::
     &  diurnal_cycle_factor



#if coupled
      integer (kind=int_kind) :: 
     &  cpl_task              ! task id for coupler

      integer (kind=int_kind) ::
     &  timer_send_to_cpl  
     &, timer_recv_from_cpl
     &, timer_recv_to_send 
     &, timer_send_to_recv 
 
      integer (kind=int_kind), private :: 
     &  cpl_stop_now          ! flag id for stop_now flag
     &, cpl_ts                ! flag id for coupled_ts flag
     &, cpl_write_restart     ! flag id for write restart
     &, cpl_write_history     ! flag id for write history
     &, cpl_write_tavg        ! flag id for write tavg      
     &, cpl_diag_global       ! flag id for computing diagnostics
     &, cpl_diag_transp       ! flag id for computing diagnostics

      integer (kind=int_kind), dimension(cpl_fields_ibuf_total) ::
     &   isbuf                ! integer control buffer for sends
     &  ,irbuf                ! integer control buffer for receives
      type(cpl_contract) ::
     &   contractS            ! contract for sends to coupler
     &  ,contractR            ! contract for receives from coupler
      real (kind=dbl_kind), dimension(:,:), allocatable ::
     &   sbuf                 ! temporary send/recv buffer

!-----------------------------------------------------------------------
!  The following variables are used in the exchange of 
!  information between cpl6 and the ocean code.
!
!  ocn --> cpl6
!  ============
!    cpl_fields_ibuf_total -- length of integer ocean "send buffer" vector (isbuf)
!    cpl_fields_o2c_total --  total number of 2D fields sent to cpl6 from ocn
!
!    integer send buffer indices (isbuf in subroutine init_coupled):  
!
!     o  cpl_fields_ibuf_cdate   -- ocean's character date string (yyyymmdd)
!     o  cpl_fields_ibuf_sec     -- ocean's character time string (seconds)
!     o  cpl_fields_ibuf_precadj -- precipitation adjustment factor*1e6
!     o  cpl_fields_ibuf_lsize   -- (iphys_e-iphys_b+1)*(jphys_e-jphys_b+1)
!     o  cpl_fields_ibuf_lisize  -- (iphys_e-iphys_b+1)
!     o  cpl_fields_ibuf_ljsize  -- (jphys_e-jphys_b+1)
!     o  cpl_fields_ibuf_gsize   -- imt_global*jmt_global
!     o  cpl_fields_ibuf_gisize  -- imt_global
!     o  cpl_fields_ibuf_gjsize  -- jmt_global
!     o  cpl_fields_ibuf_ncpl    -- ncouple_per_day
!     o  cpl_fields_ibuf_nfields -- cpl_fields_grid_total
!     o  cpl_fields_ibuf_dead    --  0 ==>  not a "dead" model
!
!    real send buffer indices (sbuf in subroutine init_coupled):
!
!     o  cpl_fields_grid_lon   -- radian*TLONG(i,j)
!     o  cpl_fields_grid_lat   -- radian*TLAT(i,j)
!     o  cpl_fields_grid_area  -- TAREA(i,j)/(radius*radius)
!     o  cpl_fields_grid_mask  -- float(REGION_MASK(i,j))
!     o  cpl_fields_grid_index -- (j_global(j)-1)*(imt_global)+i_global(i)
!
!    real send buffer indices (sbuf in subroutine send_to_coupler):
!
!      (1)  cpl_fields_o2c_u     -- surface u velocity
!      (2)  cpl_fields_o2c_v     -- surface v velocity
!      (3)  cpl_fields_o2c_t     -- surface temperature
!      (4)  cpl_fields_o2c_s     -- surface salinity
!      (5)  cpl_fields_o2c_dhdx  -- e,w surface slope
!      (6)  cpl_fields_o2c_dhdy  -- n,s surface slope
!      (7)  cpl_fields_o2c_q     -- qflux
!
!
!    cpl6 --> ocn  
!    ============
!
!      cpl_fields_ibuf_total -- length of integer ocean "receive buffer" vector (irbuf)
!      cpl_fields_c2o_total --  total number of 2D fields sent to ocn from cpl6
!
!      integer receive buffer indices (irbuf in subroutine recv_from_coupler):
!
!       o  cpl_fields_ibuf_stopnow  -- stop ocean integration now
!       o  cpl_fields_ibuf_infobug  -- write ocean/coupler diagnostics now  
!       o  cpl_fields_ibuf_resteod  -- write ocean restart files at end of day
!       o  cpl_fields_ibuf_histeod  -- write ocean history files at end of day
!       o  cpl_fields_ibuf_histtavg -- write ocean "tavg"  files at end of day
!       o  cpl_fields_ibuf_diageod  -- write ocean diagnostics   at end of day
!
!      real receive buffer indices (sbuf in subroutine recv_from_coupler):
!
!       o   cpl_fields_c2o_taux  -- zonal wind stress (taux)
!       o   cpl_fields_c2o_tauy  -- meridonal wind stress (tauy)
!       o   cpl_fields_c2o_snow  -- water flux due to snow
!       o   cpl_fields_c2o_rain  -- water flux due to rain
!       o   cpl_fields_c2o_evap  -- evaporation flux
!       o   cpl_fields_c2o_meltw -- snow melt flux
!       o   cpl_fields_c2o_roff  -- river runoff flux
!       o   cpl_fields_c2o_salt  -- salt
!       o   cpl_fields_c2o_swnet -- net short-wave heat flux
!       o   cpl_fields_c2o_sen   -- sensible heat flux
!       o   cpl_fields_c2o_lwup  -- longwave radiation (up)
!       o   cpl_fields_c2o_lwdn  -- longwave radiation (down)
!       o   cpl_fields_c2o_melth -- heat flux from snow&ice melt
!       o   cpl_fields_c2o_ifrac -- ice fraction
!       o   cpl_fields_c2o_press -- sea-level pressure
!       o   cpl_fields_c2o_duu10 -- 10m wind speed squared
!
!
!-----------------------------------------------------------------------


      real (kind=dbl_kind) ::
     &  tlast_coupled

      real (kind=dbl_kind), dimension(imt,jmt,cpl_fields_o2c_total) ::
     &  SBUFF_SUM           ! accumulated sum of send buffer quantities
                            ! for averaging before being sent

      real (kind=dbl_kind), dimension(imt,jmt) ::
     &  EVAP_F = c0         ! evaporation   flux    from cpl (kg/m2/s)
     &, PREC_F = c0         ! precipitation flux    from cpl (kg/m2/s)
                            ! (rain + snow)
     &, SNOW_F = c0         ! snow          flux    from cpl (kg/m2/s)
     &, MELT_F = c0         ! melt          flux    from cpl (kg/m2/s)
     &, ROFF_F = c0         ! river runoff  flux    from cpl (kg/m2/s)
     &, SALT_F = c0         ! salt          flux    from cpl (kg(salt)/m2/s)
     &, SENH_F = c0         ! sensible heat flux    from cpl (W/m2   )
     &, LWUP_F = c0         ! longwave heat flux up from cpl (W/m2   )
     &, LWDN_F = c0         ! longwave heat flux dn from cpl (W/m2   )
     &, MELTH_F= c0         ! melt     heat flux    from cpl (W/m2   )
 
#endif
!***********************************************************************

      contains

!***********************************************************************

      subroutine init_coupled(SMF, SMFT, STF, SHF_QSW, lsmft_avail)

!-----------------------------------------------------------------------
!
!     this routine sets up everything necessary for coupling with
!     the CCSM2 flux coupler, version 6 (cpl6)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(imt,jmt,2), intent(out) ::
     &   SMF                !  surface momentum fluxes (wind stress)
     &,  SMFT               !  surface momentum fluxes at T points

      real (kind=dbl_kind), dimension(imt,jmt,nt), intent(out) ::
     &   STF                !  surface tracer fluxes

      real (kind=dbl_kind), dimension(imt,jmt), intent(out) ::
     &   SHF_QSW            !  penetrative solar heat flux

      logical (kind=log_kind), intent(out) ::
     &   lsmft_avail        ! true if SMFT is an available field

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len) ::
     &  coupled_freq_opt

      namelist /coupled_nml/ coupled_freq_opt, coupled_freq,
     &                       qsw_diurnal_cycle

      integer (kind=int_kind) :: k,
     &  ncouple_per_day,      ! num of coupler comms per day
     &  nml_error             ! namelist i/o error flag
     &, allocate_error_status ! allocate error flag

!-----------------------------------------------------------------------
!
!     variables associated with the solar diurnal cycle
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind) ::
     &   time_for_forcing,    ! time of day for surface forcing
     &   frac_day_forcing,    ! fraction of day based on time_for_forcing
     &   cycle_function,      ! intermediate result of the diurnal cycle function
     &   weight_forcing,      ! forcing weights
     &   sum_forcing          ! sum of forcing weights

      integer (kind=int_kind) ::
     &   count_forcing        ! time step counter (== nsteps_this_interval+1)

      integer (kind=int_kind) ::
     &   i,j,n

!-----------------------------------------------------------------------
!
!     read coupled_nml namelist to start coupling and determine
!     coupling frequency
!
!-----------------------------------------------------------------------
      
      lcoupled          = .false.
      coupled_freq_opt  = 'never'
      coupled_freq_iopt = freq_opt_never
      coupled_freq      = 100000
      qsw_diurnal_cycle = .false.
      
      if (my_task == master_task) then
        nml_error = -1
        open (nml_in, file=nml_filename, status='old')
   10   continue  !*** keep reading until find right namelist
        read(nml_in, nml=coupled_nml, err=10, end=20)
        close(nml_in)
        nml_error = 0
   20   continue
      endif

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
        call exit_POP('ERROR: reading coupled_nml')
      endif


      if (my_task == master_task) then
          write(stdout,*) ' '
          write(stdout,*) ' Document Namelist Parameters:'
          write(stdout,*) ' ============================ '
          write(stdout,*) ' '
          write(stdout, coupled_nml)
          write(stdout,*) ' '
      endif

      if (my_task == master_task) then
        select case (coupled_freq_opt)

        case ('nyear')
          coupled_freq_iopt = -1000

        case ('nmonth')
          coupled_freq_iopt = -1000

        case ('nday')
          if (coupled_freq == 1) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nday
            ncouple_per_day = 1
          else
            coupled_freq_iopt = -1000
          endif

        case ('nhour')
          if (coupled_freq <= 24) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nhour
            ncouple_per_day = 24/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('nsecond')
          if (coupled_freq <= seconds_in_day) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nsecond
            ncouple_per_day = seconds_in_day/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('nstep')
          if (coupled_freq <= nsteps_per_day) then
            lcoupled = .true.
            coupled_freq_iopt = freq_opt_nstep
            ncouple_per_day = nsteps_per_day/coupled_freq
          else
            coupled_freq_iopt = -1000
          endif

        case ('never')
          lcoupled = .false.

        case default
          coupled_freq_iopt = -2000
        end select

      endif
            
      call broadcast_scalar(lcoupled,          master_task)
      call broadcast_scalar(coupled_freq_iopt, master_task)
      call broadcast_scalar(coupled_freq     , master_task)
      call broadcast_scalar(qsw_diurnal_cycle, master_task)

      if (coupled_freq_iopt == -1000) then
        call exit_POP(
     &       'ERROR: Coupling frequency must be at least once per day')
      else if (coupled_freq_iopt == -2000) then
        call exit_POP('ERROR: Unknown option for coupling frequency')
      endif

!-----------------------------------------------------------------------
!
!     register lcoupled if running with the flux coupler
!
!-----------------------------------------------------------------------

      if (lcoupled) call register_string('lcoupled')

!-----------------------------------------------------------------------
!
!     check consistency of the qsw_diurnal_cycle option with various
!     time manager options
!
!-----------------------------------------------------------------------

      if ( qsw_diurnal_cycle ) then
        if ( tmix_iopt /= tmix_avgfit )
     &    call exit_POP(
     &         'ERROR: time_mix_opt must be set to avgfit '//
     &         '       for solar diurnal cycle' ) 

        if ( dttxcel(1) /= c1  .or.  dtuxcel /= c1 ) 
     &    call exit_POP(
     &         'ERROR: using the specified accelerated integration '//
     &         '       technique may not be appropriate for solar '//
     &         '       diurnal cycle' )
      endif

!-----------------------------------------------------------------------
!
!     allocate and compute the short wave heat flux multiplier for the 
!     diurnal cycle
!
!-----------------------------------------------------------------------

      allocate ( diurnal_cycle_factor(nsteps_per_interval), 
     &           stat=allocate_error_status )
      call check_allocate ( allocate_error_status,
     &       exit_message = 'diurnal_cycle in forcing_coupled.F' )

      diurnal_cycle_factor = c1
      if ( qsw_diurnal_cycle ) then

!       mimic a day

        time_for_forcing = c0 
        count_forcing    =  1
        sum_forcing      = c0

        do n=1,nsteps_per_interval
          frac_day_forcing = time_for_forcing / seconds_in_day 
          cycle_function = cos( pi * ( c2 * frac_day_forcing - c1 ) )
          diurnal_cycle_factor(n) = c2 * ( cycle_function 
     &                               + abs(cycle_function) ) 
     &                               * cycle_function
          weight_forcing = c1
          if (  count_forcing == 2  .or.
     &      mod(count_forcing,time_mix_freq) == 0 )
     &      weight_forcing = p5
          time_for_forcing = time_for_forcing + weight_forcing * dt(1)
          sum_forcing = sum_forcing 
     &              + weight_forcing * dt(1) * diurnal_cycle_factor(n)
          count_forcing = count_forcing + 1
        enddo

        diurnal_cycle_factor = diurnal_cycle_factor * seconds_in_day 
     &                        / sum_forcing

!       check the final integral

        count_forcing =  1
        sum_forcing   = c0

        do n=1,nsteps_per_interval
          weight_forcing = c1
          if (  count_forcing == 2  .or.
     &      mod(count_forcing,time_mix_freq) == 0 )
     &      weight_forcing = p5
          sum_forcing = sum_forcing
     &              + weight_forcing * dt(1) * diurnal_cycle_factor(n)
          count_forcing = count_forcing + 1
        enddo

        if ( sum_forcing < (seconds_in_day - 1.0e-5_dbl_kind)  .or.
     &       sum_forcing > (seconds_in_day + 1.0e-5_dbl_kind) )
     &    call exit_POP ('ERROR: qsw diurnal cycle temporal '//
     &                   'integral is incorrect')

      endif
      
#if coupled
      if (.not. lcoupled) then
        call exit_POP(
     &       'ERROR: Coupled ifdef option enabled but lcoupled=false')
      endif

!-----------------------------------------------------------------------
!
!     Initialize flags and shortwave absorption profile
!     Note that the cpl_write_xxx flags have _no_ default value;
!     therefore, they must be explicitly set .true. and .false.
!     at the appropriate times
!
!-----------------------------------------------------------------------

      cpl_stop_now      = init_time_flag('stop_now',default=.false.)
      cpl_ts            = init_time_flag('coupled_ts',
     &                                   freq_opt = coupled_freq_iopt,
     &                                   freq     = coupled_freq)
      cpl_write_restart = init_time_flag('cpl_write_restart')
      cpl_write_history = init_time_flag('cpl_write_history')
      cpl_write_tavg    = init_time_flag('cpl_write_tavg'   )
      cpl_diag_global   = init_time_flag('cpl_diag_global')
      cpl_diag_transp   = init_time_flag('cpl_diag_transp')

      lsmft_avail = .true.
      tlast_coupled = c0


!-----------------------------------------------------------------------
!
!       initialize and send buffer
!
!-----------------------------------------------------------------------

      isbuf = 0

      isbuf(cpl_fields_ibuf_cdate  ) = iyear*10000 + imonth*100 + iday
      isbuf(cpl_fields_ibuf_sec    ) = 
     &   ihour*seconds_in_hour + iminute*seconds_in_minute + isecond
      isbuf(cpl_fields_ibuf_lsize  ) = (iphys_e-iphys_b+1)*
     &                                 (jphys_e-jphys_b+1)
      isbuf(cpl_fields_ibuf_lisize ) = (iphys_e-iphys_b+1)
      isbuf(cpl_fields_ibuf_ljsize ) = (jphys_e-jphys_b+1)
      isbuf(cpl_fields_ibuf_gsize  ) = imt_global*jmt_global
      isbuf(cpl_fields_ibuf_gisize ) = imt_global
      isbuf(cpl_fields_ibuf_gjsize ) = jmt_global
      isbuf(cpl_fields_ibuf_ncpl   ) = ncouple_per_day
      isbuf(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
      isbuf(cpl_fields_ibuf_dead   ) = 0           ! not a dead model

      allocate(sbuf((iphys_e-iphys_b+1)*(jphys_e-jphys_b+1)
     &,        cpl_fields_grid_total))
      sbuf = -888.0
      n=0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n=n+1
         sbuf(n,cpl_fields_grid_lon  ) = radian*TLONG(i,j)
         sbuf(n,cpl_fields_grid_lat  ) = radian*TLAT(i,j)
         sbuf(n,cpl_fields_grid_area ) = TAREA(i,j)/(radius*radius)
         sbuf(n,cpl_fields_grid_mask ) = float(REGION_MASK(i,j))
         sbuf(n,cpl_fields_grid_index) =   
     &      (j_global(j)-1)*(imt_global) + i_global(i)
      enddo
      enddo


!-----------------------------------------------------------------------
!     initialize the contracts
!-----------------------------------------------------------------------

      call cpl_interface_contractInit(contractS,cpl_fields_ocnname,
     &   cpl_fields_cplname,cpl_fields_o2c_fields,isbuf,sbuf)
      call cpl_interface_contractInit(contractR,cpl_fields_ocnname,
     &   cpl_fields_cplname,cpl_fields_c2o_fields,isbuf,sbuf)

      write(stdout,*) '(ocn) Initialized contracts with coupler'
      call shr_sys_flush(stdout)

      deallocate(sbuf)

      !--- receive initial message from coupler
      call cpl_interface_ibufRecv(cpl_fields_cplname,irbuf)

!-----------------------------------------------------------------------
!
!     send initial state info to coupler
!
!-----------------------------------------------------------------------

      call sum_buffer

      call send_to_coupler

!-----------------------------------------------------------------------
!
!     initialize timers for coupled model
!
!-----------------------------------------------------------------------
      
      call get_timer (timer_send_to_cpl  , 'SEND'        )
      call get_timer (timer_recv_from_cpl, 'RECV'        )
      call get_timer (timer_recv_to_send , 'RECV to SEND')
      call get_timer (timer_send_to_recv , 'SEND to RECV')
#endif
!-----------------------------------------------------------------------

      end subroutine init_coupled

!***********************************************************************

      subroutine set_coupled_forcing(SMF,SMFT,STF,SHF_QSW,IFRAC,
     &  ATM_PRESS, U10_SQR)


!-----------------------------------------------------------------------
!
!     this routine call coupler communication routines to set
!     surface forcing data
!
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     Note: We are using intent "inout" for SMF,SMFT, STF, SHF_QSW
!           and IFRAC in order to preserve their values inbetween
!           coupling timesteps.
!-----------------------------------------------------------------------
 
 
      real (kind=dbl_kind), dimension(imt,jmt,2), intent(inout) ::
     &   SMF                !  surface momentum fluxes (wind stress)
     &,  SMFT               !  surface momentum fluxes at T points

      real (kind=dbl_kind), dimension(imt,jmt,nt), intent(inout) ::
     &   STF                !  surface tracer fluxes

      real (kind=dbl_kind), dimension(imt,jmt), intent(inout) ::
     &   SHF_QSW            !  penetrative solar heat flux
     &,  IFRAC              !  fractional ice coverage
     &,  ATM_PRESS          !  atmospheric pressure forcing
     &,  U10_SQR            !  10m wind speed squared

     
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n

           
 
#if coupled
!-----------------------------------------------------------------------
!
!     if it is time to couple, exchange data with flux coupler
!     be sure to trigger communication on very first time step
!
!-----------------------------------------------------------------------

      if (nsteps_run /= 0) call sum_buffer

      if (check_time_flag(cpl_ts) .or. nsteps_run == 0) then

        !*** send state variables at end of coupling interval 
 

        if (nsteps_run /= 0) then
           call timer_stop  (timer_recv_to_send)
 
           call timer_start (timer_send_to_cpl)
           call send_to_coupler
           call timer_stop  (timer_send_to_cpl)
        endif

        call timer_start (timer_send_to_recv)
        call timer_stop  (timer_send_to_recv) 
 
        !*** recv data to advance next time step

        call timer_start (timer_recv_from_cpl)
        call recv_from_coupler(SMF,SMFT,STF,SHF_QSW,IFRAC,ATM_PRESS,
     &                         U10_SQR)
        call timer_stop  (timer_recv_from_cpl)

        call timer_start (timer_recv_to_send)
 

        if ( shf_formulation == 'partially-coupled' ) then
          SHF_COMP(:,:,shf_comp_cpl) = STF(:,:,1) 
          if ( .not. lms_balance ) then
            SHF_COMP(:,:,shf_comp_cpl) = 
     &              SHF_COMP(:,:,shf_comp_cpl) * MASK_SR
            SHF_QSW = SHF_QSW * MASK_SR
          endif
        endif
 
        SHF_COMP(:,:,shf_comp_qsw) = SHF_QSW

        if ( sfwf_formulation == 'partially-coupled' ) then

          if (sfc_layer_type == sfc_layer_varthick .and.
     &        .not. lfw_as_salt_flx) then
            SFWF_COMP(:,:,sfwf_comp_cpl) =  FW * MASK_SR
            do n=1,nt
              TFW_COMP(:,:,n,tfw_comp_cpl) = TFW(:,:,n) * MASK_SR
            enddo
          else
            SFWF_COMP(:,:,sfwf_comp_cpl) = STF(:,:,2) * MASK_SR
          endif

        else

          if ( sfc_layer_type == sfc_layer_varthick .and.
     &         .not. lfw_as_salt_flx .and. liceform ) then
            SFWF_COMP(:,:,sfwf_comp_cpl) = FW
            TFW_COMP(:,:,:,tfw_comp_cpl) = TFW
          endif

        endif
 
        if ( luse_cpl_ifrac ) then
          OCN_WGT = (c1-IFRAC) * RCALCT
        endif
 
      endif

#endif
!-----------------------------------------------------------------------

      end subroutine set_coupled_forcing

!***********************************************************************

      subroutine set_combined_forcing (STF)

!-----------------------------------------------------------------------
!
!     this routine combines terms when the "partially-coupled"
!     has been selected
!
!-----------------------------------------------------------------------
 
      real (kind=dbl_kind), dimension(imt,jmt,nt) ::
     &   STF    !  surface tracer fluxes
 
      real (kind=dbl_kind), dimension(imt,jmt) :: WORK1, WORK2
 
      
#if coupled

      if ( shf_formulation == 'partially-coupled' ) then
        STF(:,:,1) =  SHF_COMP(:,:,shf_comp_wrest)
     &              + SHF_COMP(:,:,shf_comp_srest)
     &              + SHF_COMP(:,:,shf_comp_cpl)
      endif

      if ( sfwf_formulation == 'partially-coupled' ) then
        if (sfc_layer_type == sfc_layer_varthick .and.
     &      .not. lfw_as_salt_flx) then
           STF(:,:,2) =  SFWF_COMP(:,:,sfwf_comp_wrest)
     &                 + SFWF_COMP(:,:,sfwf_comp_srest)
           FW         =  SFWF_COMP(:,:,sfwf_comp_cpl)
     &                 + SFWF_COMP(:,:,sfwf_comp_flxio)
           TFW        =  TFW_COMP(:,:,:,tfw_comp_cpl)
     &                 + TFW_COMP(:,:,:,tfw_comp_flxio)
        else
          if ( lms_balance ) then

            WORK1 = SFWF_COMP(:,:,sfwf_comp_flxio) / salinity_factor
            WORK2 = SFWF_COMP(:,:,sfwf_comp_cpl)

            call ms_balancing ( WORK2, EVAP_F, PREC_F, MELT_F,
     &                          ROFF_F, SALT_F, 'salt',
     &                          ICEOCN_F=WORK1 )

            STF(:,:,2) =  SFWF_COMP(:,:,sfwf_comp_wrest)
     &                  + SFWF_COMP(:,:,sfwf_comp_srest)
     &                  + WORK2 
     &                  + SFWF_COMP(:,:,sfwf_comp_flxio) * MASK_SR
          else     

            STF(:,:,2) =  SFWF_COMP(:,:,sfwf_comp_wrest)
     &                  + SFWF_COMP(:,:,sfwf_comp_srest)
     &                  + SFWF_COMP(:,:,sfwf_comp_cpl)
     &                  + SFWF_COMP(:,:,sfwf_comp_flxio) 
 
          endif
        endif
      endif

#endif
!-----------------------------------------------------------------------

      end subroutine set_combined_forcing

#if coupled
!***********************************************************************

      subroutine recv_from_coupler(SMF,SMFT,STF,SHF_QSW,IFRAC,ATM_PRESS,
     &  U10_SQR)

!-----------------------------------------------------------------------
!
!     this routine receives message from coupler with surface flux
!     data
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(imt,jmt,2), intent(out) ::
     &   SMF                !  surface momentum fluxes (wind stress)
     &,  SMFT               !  surface momentum fluxes at T points

      real (kind=dbl_kind), dimension(imt,jmt,nt), intent(out) ::
     &   STF                !  surface tracer fluxes

      real (kind=dbl_kind), dimension(imt,jmt), intent(out) ::
     &   SHF_QSW            !  penetrative solar heat flux
     &,  IFRAC              !  fractional ice coverage
     &,  ATM_PRESS          !  atmospheric pressure forcing
     &,  U10_SQR            !  10m wind speed squared

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len)    :: label
 
      integer (kind=int_kind) ::
     &  i,j,k,n

      real (kind=dbl_kind), dimension(imt,jmt) :: 
     &  WORK1, WORK2        ! local work space
     &, WORK3, WORK4        ! local work space
     &, WORK5               ! local work space

      real (kind=dbl_kind) :: 
     &   m2percm2
     &,  gsum

!-----------------------------------------------------------------------
!
!     receive message from coupler and check for terminate signal
!
!-----------------------------------------------------------------------

      allocate(sbuf((iphys_e-iphys_b+1)*(jphys_e-jphys_b+1)
     &,       cpl_fields_c2o_total))
      call cpl_interface_contractRecv(cpl_fields_cplname,contractR
     &,       irbuf,sbuf)

!-----------------------------------------------------------------------
!
!     check all coupler flags and respond appropriately
!
!-----------------------------------------------------------------------

      if (irbuf(cpl_fields_ibuf_stopnow) == 1) then
        call set_time_flag(cpl_stop_now,.true.)
        if (my_task == master_task) then
          call int_to_char (4,iyear   , cyear  )
          call int_to_char (2,imonth  , cmonth )
          call int_to_char (2,iday    , cday   )
          call int_to_char (2,ihour   , chour  )
          call int_to_char (2,iminute , cminute)
          call int_to_char (2,isecond , csecond)
          write(stdout,*) ' cpl requests termination now:'             
     &,                     ' ', cyear,'/',cmonth,'/',cday
     &,                     ' ', chour,':',cminute,':',csecond
        endif
       !call exit_POP('Received terminate signal from coupler')
        RETURN
      endif


      if (irbuf(cpl_fields_ibuf_infobug) >= 2) then
         ldiag_cpl = .true. 
      else
         ldiag_cpl = .false.
      endif

      if (irbuf(cpl_fields_ibuf_resteod) == 1) then
         call set_time_flag(cpl_write_restart,.true.)
         if (my_task == master_task) then
           write(stdout,*) ' cpl requests restart file at eod '
     &,                     ' ', cyear,'/',cmonth,'/',cday, '  '
         endif
      endif
  
!      if (irbuf(cpl_fields_ibuf_histeod) == 1) then
!       ignore for now
!        call set_time_flag(cpl_write_history,.true.)
!        call int_to_char (4,iyear   , cyear )
!        call int_to_char (2,imonth  ,cmonth )
!        call int_to_char (2,iday    ,cday   )
!        call int_to_char (2,ihour   ,chour  )
!        call int_to_char (2,iminute ,cminute)
!        call int_to_char (2,isecond ,csecond)
!        if (my_task == master_task) then
!        write(stdout,*) ' cpl requests history file at eod '
!    &,                     ' ', cyear,'/',cmonth,'/',cday, '  '
!        endif
!      endif
  
  
!      if (irbuf(cpl_fields_ibuf_histtavg) == 1) then
!       ignore for now
!        call set_time_flag(cpl_write_tavg, .true.)
!        call int_to_char (4,iyear   , cyear )
!        call int_to_char (2,imonth  ,cmonth )
!        call int_to_char (2,iday    ,cday   )
!        call int_to_char (2,ihour   ,chour  )
!        call int_to_char (2,iminute ,cminute)
!        call int_to_char (2,isecond ,csecond)
!        if (my_task == master_task) then
!        write(stdout,*) ' cpl requests tavg file at eod '
!    &,                     ' ', cyear,'/',cmonth,'/',cday, '  '
!        endif
!      endif
  
       if (irbuf(cpl_fields_ibuf_diageod) == 1) then
         call set_time_flag(cpl_diag_global,.true.)
         call set_time_flag(cpl_diag_transp,.true.)
 
         call int_to_char (4,iyear   ,cyear  )
         call int_to_char (2,imonth  ,cmonth )
         call int_to_char (2,iday    ,cday   )
         call int_to_char (2,ihour   ,chour  )
         call int_to_char (2,iminute ,cminute)
         call int_to_char (2,isecond ,csecond)
 
         if (my_task == master_task) then
         write(stdout,*) ' cpl requests diagnostics at eod '
     &,                     ' ', cyear,'/',cmonth,'/',cday, '  '
         endif
       endif

!-----------------------------------------------------------------------
!
!     unpack and distribute wind stress, then convert to correct units
!     if necessary, rotate components to logical coordinates
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         WORK1(i,j) = sbuf(n,cpl_fields_c2o_taux)
         WORK2(i,j) = sbuf(n,cpl_fields_c2o_tauy)
      enddo
      enddo
      call boundary_2d(WORK1)
      call boundary_2d(WORK2)

      !***
      !*** Rotate true zonal/meridional wind stress into local
      !*** coordinates and convert to dyne/cm**2
      !***

      SMFT(:,:,1) = (WORK1*cos(ANGLET) + WORK2*sin(ANGLET))*
     &              RCALCT*momentum_factor
      SMFT(:,:,2) = (WORK2*cos(ANGLET) - WORK1*sin(ANGLET))*
     &              RCALCT*momentum_factor
 
      !***
      !*** Shift SMFT to U grid
      !***

      call ne_4pt_nobndy(SMF (:,:,1),AU0,AUN,AUE,AUNE,SMFT(:,:,1))
      call ne_4pt_nobndy(SMF (:,:,2),AU0,AUN,AUE,AUNE,SMFT(:,:,2))


!-----------------------------------------------------------------------
!
!     unpack and distribute fresh water flux and salt flux
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         SNOW_F(i,j) = sbuf(n,cpl_fields_c2o_snow)
         WORK1(i,j)  = sbuf(n,cpl_fields_c2o_rain)
         EVAP_F(i,j) = sbuf(n,cpl_fields_c2o_evap)
         MELT_F(i,j) = sbuf(n,cpl_fields_c2o_meltw)
         ROFF_F(i,j) = sbuf(n,cpl_fields_c2o_roff)
         SALT_F(i,j) = sbuf(n,cpl_fields_c2o_salt)
      enddo
      enddo

      PREC_F = WORK1 + SNOW_F         ! rain + snow

      call boundary_2d(SNOW_F)
      call boundary_2d(PREC_F)
      call boundary_2d(EVAP_F)
      call boundary_2d(MELT_F)
      call boundary_2d(ROFF_F)
      call boundary_2d(SALT_F)

!-----------------------------------------------------------------------
!
!     for variable thickness surface layer, compute fresh water and
!     salt fluxes
!
!-----------------------------------------------------------------------

      if (sfc_layer_type == sfc_layer_varthick .and. 
     &    .not. lfw_as_salt_flx) then

        !*** compute fresh water flux (cm/s)

        FW = RCALCT*(PREC_F+EVAP_F+ROFF_F)*fwmass_to_fwflux

        WORK1 = RCALCT * MELT_F * fwmass_to_fwflux

        !*** compute tracer concentration in fresh water
        !*** in principle, temperature of each water flux
        !*** could be different. e.g.
        !TFW(:,:,1) = RCALCT*(PREC_F*TEMP_PREC +
        !                     EVAP_F*TEMP_EVAP +
        !                     MELT_F*TEMP_MELT +
        !                     ROFF_F*TEMP_ROFF)*fwmass_to_fwflux
        !*** currently assume water comes in at sea surface temp

        TFW(:,:,1) = FW*TRACER(:,:,1,1,curtime)
     &              + WORK1 * TMLT(TRACER(:,:,1,2,curtime))

        FW = FW + WORK1

        !*** compute salt flux
        !*** again, salinity could be different for each
        !***   component of water flux
        !TFW(:,:,2) = RCALCT*(PREC_F*SALT_PREC +
        !                     EVAP_F*SALT_EVAP +
        !                     MELT_F*SALT_MELT +
        !                     ROFF_F*SALT_ROFF)*fwmass_to_fwflux
        !*** currently assume prec, evap and roff are fresh
        !*** and all salt come from ice melt

        where (MELT_F /= c0)
          WORK1 = SALT_F/MELT_F ! salinity (msu) of melt water
        elsewhere
          WORK1 = c0
        end where

        TFW(:,:,2) =  RCALCT*MELT_F*fwmass_to_fwflux*WORK1
                      ! + PREC_F*c0 + EVAP_F*c0 + ROFF_F*c0

        do n=3,nt
          TFW(:,:,n) = c0  ! no additional tracers in fresh water
        end do

!-----------------------------------------------------------------------
!
!     if not a variable thickness surface layer or if fw_as_salt_flx
!     flag is on, convert fresh and salt inputs to a virtual salinity
!     flux
!
!-----------------------------------------------------------------------

      else  ! convert fresh water to virtual salinity flux

        STF(:,:,2) = RCALCT*(
     &               (PREC_F+EVAP_F+MELT_F+ROFF_F)*salinity_factor 
     &              + SALT_F*sflux_factor)  
 
!-----------------------------------------------------------------------
!
!     balance salt/freshwater in marginal seas
!
!-----------------------------------------------------------------------
 
      if  (lms_balance .and. sfwf_formulation /= 'partially-coupled' )
     & call ms_balancing (STF(:,:,2),EVAP_F, PREC_F, MELT_F,ROFF_F
     &,                   SALT_F, 'salt')
 
      endif
 
!-----------------------------------------------------------------------
!
!     unpack and distribute net shortwave and convert units
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         WORK1(i,j) = sbuf(n,cpl_fields_c2o_swnet)
      enddo
      enddo
      call boundary_2d(WORK1)
      SHF_QSW = WORK1*RCALCT*hflux_factor  !  convert from W/m**2

!-----------------------------------------------------------------------
!
!     unpack and distribute non-sw heat flux (lat+sen+lw+melt+snow)
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         WORK1(i,j) = sbuf(n,cpl_fields_c2o_sen)
         WORK2(i,j) = sbuf(n,cpl_fields_c2o_lwup)
         WORK3(i,j) = sbuf(n,cpl_fields_c2o_lwdn)
         WORK4(i,j) = sbuf(n,cpl_fields_c2o_melth)
      enddo
      enddo
      call boundary_2d(WORK1)
      call boundary_2d(WORK2)
      call boundary_2d(WORK3)
      call boundary_2d(WORK4)

      SENH_F  = WORK1
      LWUP_F  = WORK2
      LWDN_F  = WORK3
      MELTH_F = WORK4

      WORK5 = - SNOW_F * latent_heat_fusion_mks

!-----------------------------------------------------------------------
!
!     combine and convert from W/m**2
!           (note: latent heat flux = evaporation*latent_heat_vapor)
!
!-----------------------------------------------------------------------

      STF(:,:,1) = (EVAP_F*latent_heat_vapor + WORK1 + WORK2 
     &             + WORK3 + WORK4 + WORK5)
     &             *RCALCT*hflux_factor 

!-----------------------------------------------------------------------
!
!     unpack and distribute fractional ice coverage
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         WORK1(i,j) = sbuf(n,cpl_fields_c2o_ifrac)
      enddo
      enddo
      call boundary_2d(WORK1)

      IFRAC = WORK1 * RCALCT

!-----------------------------------------------------------------------
!
!     unpack and distribute atmospheric pressure forcing
!     converting from Pa to dynes/cm**2
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         WORK1(i,j) = sbuf(n,cpl_fields_c2o_press)
      enddo
      enddo
      call boundary_2d(WORK1)

      ATM_PRESS = c10 * WORK1 * RCALCT

!-----------------------------------------------------------------------
!
!     unpack and distribute 10m wind speed squared
!     converting from m**2/s**2 to cm**2/s**2
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         WORK1(i,j) = sbuf(n,cpl_fields_c2o_duu10)
      enddo
      enddo
      call boundary_2d(WORK1)

      U10_SQR = cmperm * cmperm * WORK1 * RCALCT

!-----------------------------------------------------------------------
!
!     diagnostics
!
!-----------------------------------------------------------------------

      if (ldiag_cpl) then
        if  (my_task == master_task) then
          call int_to_char (4,iyear   ,cyear  )
          call int_to_char (2,imonth  ,cmonth )
          call int_to_char (2,iday    ,cday   )
          call int_to_char (2,ihour   ,chour  )
          call int_to_char (2,iminute ,cminute)
          call int_to_char (2,isecond ,csecond)
          write(stdout,*)' Global averages of fluxes received from cpl'
     &,                  ' at ', cyear,'/',cmonth ,'/',cday  
     &,                     ' ', chour,':',cminute,':',csecond
          call shr_sys_flush(stdout)
       endif
 
         m2percm2  = mpercm*mpercm
         do k = 1,cpl_fields_c2o_total
            n = 0
            do j=jphys_b,jphys_e
            do i=iphys_b,iphys_e
               n = n + 1
               WORK1(i,j) = sbuf(n,k)*TAREA(i,j)
            enddo
            enddo
            call boundary_2d(WORK1)
            gsum = global_sum(WORK1,RCALCT)*m2percm2
            if (my_task == master_task) then
               call cpl_fields_getField(label,k,cpl_fields_c2o_fields)
               write(stdout,1100)'ocn','recv', label ,gsum
            endif
         enddo
         if (my_task == master_task) call shr_sys_flush(stdout)
      endif

1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

      deallocate(sbuf)

!-----------------------------------------------------------------------

      end subroutine recv_from_coupler

!***********************************************************************

      subroutine send_to_coupler

!-----------------------------------------------------------------------
!
!     this routine packs fields into a message buffer and sends the
!     message to the flux coupler
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      character (char_len)    :: label
 
      integer (kind=int_kind) ::
     &  i,j,k,n

      real (kind=dbl_kind), dimension(imt,jmt) :: 
     &  WORK1, WORK2        ! local work space
     &, WORK3, WORK4

      real (kind=dbl_kind) :: 
     &   m2percm2
     &,  gsum

!-----------------------------------------------------------------------
!
!     initialize control buffer
!
!-----------------------------------------------------------------------

      allocate(sbuf((iphys_e-iphys_b+1)*(jphys_e-jphys_b+1)
     &,       cpl_fields_o2c_total))

      isbuf = 0

      if (check_time_flag(cpl_stop_now)) then
        isbuf(cpl_fields_ibuf_stopnow) = 1
      endif


      isbuf(cpl_fields_ibuf_cdate) = iyear*10000 + imonth*100 + iday
      isbuf(cpl_fields_ibuf_sec) = 
     &   ihour*seconds_in_hour + iminute*seconds_in_minute + isecond

      if ( lsend_precip_fact )
     &  isbuf(cpl_fields_ibuf_precadj) = precip_fact * 1.0e6_dbl_kind   ! send real as integer
 
!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points and rotate on T grid
!
!-----------------------------------------------------------------------

      call sw_4pt_nobndy(WORK3,AT0,ATS,ATW,ATSW,
     &                   SBUFF_SUM(:,:,cpl_fields_o2c_u))
      call sw_4pt_nobndy(WORK4,AT0,ATS,ATW,ATSW,
     &                   SBUFF_SUM(:,:,cpl_fields_o2c_v))

      WORK1 = (WORK3*cos(ANGLET)+WORK4*sin(-ANGLET))
     &       * mpercm/tlast_coupled
      WORK2 = (WORK4*cos(ANGLET)-WORK3*sin(-ANGLET))
     &       * mpercm/tlast_coupled

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         sbuf(n,cpl_fields_o2c_u) = WORK1(i,j)
         sbuf(n,cpl_fields_o2c_v) = WORK2(i,j)
      enddo
      enddo

!-----------------------------------------------------------------------
!
!     convert and pack surface temperature
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         sbuf(n,cpl_fields_o2c_t) = 
     &       SBUFF_SUM(i,j,cpl_fields_o2c_t)/tlast_coupled + T0_Kelvin
      enddo
      enddo

!-----------------------------------------------------------------------
!
!     convert and pack salinity
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         sbuf(n,cpl_fields_o2c_s) = 
     &       SBUFF_SUM(i,j,cpl_fields_o2c_s)*salt_to_ppt/tlast_coupled
      enddo
      enddo

!-----------------------------------------------------------------------
!
!     interpolate onto T-grid points, then rotate on T grid
!
!-----------------------------------------------------------------------

      call sw_4pt_nobndy(WORK3,AT0,ATS,ATW,ATSW,
     &                   SBUFF_SUM(:,:,cpl_fields_o2c_dhdx))
 
      call sw_4pt_nobndy(WORK4,AT0,ATS,ATW,ATSW,
     &                   SBUFF_SUM(:,:,cpl_fields_o2c_dhdy))
 
      WORK1 = (WORK3*cos(ANGLET) + WORK4*sin(-ANGLET))
     &        /grav/tlast_coupled
      WORK2 = (WORK4*cos(ANGLET) - WORK3*sin(-ANGLET))
     &        /grav/tlast_coupled


      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         sbuf(n,cpl_fields_o2c_dhdx) = WORK1(i,j)
         sbuf(n,cpl_fields_o2c_dhdy) = WORK2(i,j)
      enddo
      enddo

!-----------------------------------------------------------------------
!
!     pack heat flux due to freezing/melting (W/m^2)
!     QFLUX computation and units conversion occurs in ice.F
!
!-----------------------------------------------------------------------

      n = 0
      do j=jphys_b,jphys_e
      do i=iphys_b,iphys_e
         n = n + 1
         sbuf(n,cpl_fields_o2c_q) = QFLUX(i,j)
      enddo
      enddo

      tlast_ice = c0
      AQICE     = c0
      QICE      = c0

!-----------------------------------------------------------------------
!
!     send fields to coupler
!
!-----------------------------------------------------------------------

      call cpl_interface_contractSend(cpl_fields_cplname,contractS,isbuf,sbuf)
 
!-----------------------------------------------------------------------
!
!     diagnostics
!
!-----------------------------------------------------------------------

      if (ldiag_cpl) then
        if (my_task == master_task) then
          call int_to_char (4,iyear   ,cyear  )
          call int_to_char (2,imonth  ,cmonth )
          call int_to_char (2,iday    ,cday   )
          call int_to_char (2,ihour   ,chour  )
          call int_to_char (2,iminute ,cminute)
          call int_to_char (2,isecond ,csecond)
          write(stdout,*) ' Global averages of fluxes sent to cpl at '
     &,                   ' ', cyear,'/',cmonth, '/',cday  
     &,                   ' ', chour,':',cminute,':',csecond
          call shr_sys_flush(stdout)
        endif
 
         m2percm2  = mpercm*mpercm
         do k = 1,cpl_fields_o2c_total
            n = 0
            do j=jphys_b,jphys_e
            do i=iphys_b,iphys_e
               n = n + 1
               WORK1(i,j) = sbuf(n,k)*TAREA(i,j)
            enddo
            enddo
            call boundary_2d(WORK1)
            gsum = global_sum(WORK1,RCALCT)*m2percm2
            if (my_task == master_task) then
               call cpl_fields_getField(label,k,cpl_fields_o2c_fields)
               write(stdout,1100)'ocn','send', label ,gsum
            endif
         enddo
         if (my_task == master_task) call shr_sys_flush(stdout)
      endif

1100  format ('comm_diag ', a3, 1x, a4, 1x, a8, 1x, es26.19:, 1x, a6)

      tlast_coupled = c0

      deallocate(sbuf)

!-----------------------------------------------------------------------

      end subroutine send_to_coupler

!***********************************************************************

      subroutine sum_buffer

!-----------------------------------------------------------------------
!
!     this routine accumulates sums for averaging fields to
!     be sent to the coupler
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind) :: 
     &  delt                ! time interval since last step

!-----------------------------------------------------------------------
!
!     zero buffer if this is the first time after a coupling interval
!
!-----------------------------------------------------------------------

      if (tlast_coupled == c0) SBUFF_SUM = c0

!-----------------------------------------------------------------------
!
!     update time since last coupling
!
!-----------------------------------------------------------------------

      if (avg_ts .or. back_to_back) then
        delt = p5*dtt
      else
        delt =    dtt
      endif
      tlast_coupled = tlast_coupled + delt

!-----------------------------------------------------------------------
!
!     accumulate sums of U,V,T,S and GRADP
!     ice formation flux is handled separately in ice routine
!
!-----------------------------------------------------------------------

      SBUFF_SUM(:,:,cpl_fields_o2c_u) = SBUFF_SUM(:,:,cpl_fields_o2c_u) + delt*
     &                             UVEL(:,:,1,curtime)

      SBUFF_SUM(:,:,cpl_fields_o2c_v) = SBUFF_SUM(:,:,cpl_fields_o2c_v) + delt*
     &                             VVEL(:,:,1,curtime)

      SBUFF_SUM(:,:,cpl_fields_o2c_t ) = SBUFF_SUM(:,:,cpl_fields_o2c_t ) + delt*
     &                             TRACER(:,:,1,1,curtime)

      SBUFF_SUM(:,:,cpl_fields_o2c_s ) = SBUFF_SUM(:,:,cpl_fields_o2c_s ) + delt*
     &                             TRACER(:,:,1,2,curtime)

      SBUFF_SUM(:,:,cpl_fields_o2c_dhdx) = SBUFF_SUM(:,:,cpl_fields_o2c_dhdx) + delt*
     &                             GRADPX(:,:,curtime)

      SBUFF_SUM(:,:,cpl_fields_o2c_dhdy) = SBUFF_SUM(:,:,cpl_fields_o2c_dhdy) + delt*
     &                             GRADPY(:,:,curtime)

      end subroutine sum_buffer

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ocn_coupling_setup - sets mpi communicators and task ids
!
! !INTERFACE:
!
      subroutine ocn_coupling_setup
!
! !DESCRIPTION:
!
! This routine sets the model communicator from ccsm share code
!
! !REVISION HISTORY:
!
! author: T. Craig, NCAR, Jan 07, 2002: for cpl6
!
!
!EOP

      !call cpl_interface_init(cpl_fields_ocnname,MPI_COMM_OCN)
      call cpl_comm_init_wo_mph(cpl_fields_ocnname,MPI_COMM_OCN)

      end subroutine ocn_coupling_setup

!-----------------------------------------------------------------------
#endif
 
 
!***********************************************************************

      end module forcing_coupled

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

