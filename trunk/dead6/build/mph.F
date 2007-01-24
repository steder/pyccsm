! Copyright ©2002 The Regents of the University of California.
! All Rights Reserved. 

! Authors: 
! Yun (Helen) He, Lawrence Berkeley National Laboratory, yhe@lbl.gov
! Chris H.Q. Ding,  Lawrence Berkeley National Laboratory, chqding@lbl.gov

! THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE COPYRIGHTED BY 
! THE REGENTS OF THE UNIVERSITY OF CALIFORNIA. THE SOFTWARE PROGRAM 
! AND DOCUMENTATION ARE SUPPLIED "AS IS", WITHOUT ANY ACCOMPANYING 
! SERVICES FROM THE REGENTS. THE REGENTS DOES NOT WARRANT THAT THE 
! OPERATION OF THE PROGRAM WILL BE UNINTERRUPTED OR ERROR-FREE. 
! THE END-USER UNDERSTANDS THAT THE PROGRAM WAS DEVELOPED FOR 
! RESEARCH PURPOSES AND IS ADVISED NOT TO RELY EXCLUSIVELY ON THE 
! PROGRAM FOR ANY REASON. 

!---------------------------------------------------------------------
!BOI

! !TITLE: Multi Program-Components Handshaking (MPH) Utility Version 3 User's Manual
      
! !AUTHORS: Yun He and Chris Ding 

! !AFFILIATION: Lawrence Berkeley National Laboratory (yhe@lbl.gov, chqding@lbl.gov)                 

! !DATE: June 20, 2002


! !COPYRIGHT:

! Copyright ©2002 The Regents of the University of California.
! All Rights Reserved. 

! Permission to use, copy, modify, and distribute this software and
! its documentation for educational, research and non-profit purposes,
! without fee, and without a written agreement is hereby granted,
! provided that the above copyright notice, this paragraph and the
! following three paragraphs appear in all copies. 

! Permission to incorporate this software into commercial products
! may be obtained by contacting the University of California, Charles
! Rzeszutko, Campus Liaison Officer Office of Technology Transfer,
! 1111 Franklin Street, 5th Floor, Oakland, CA 94607, 
! (510) 587-6063, charles.rzeszutko@ucop.edu. 

! THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE COPYRIGHTED BY 
! THE REGENTS OF THE UNIVERSITY OF CALIFORNIA. THE SOFTWARE PROGRAM 
! AND DOCUMENTATION ARE SUPPLIED "AS IS", WITHOUT ANY ACCOMPANYING 
! SERVICES FROM THE REGENTS. THE REGENTS DOES NOT WARRANT THAT THE 
! OPERATION OF THE PROGRAM WILL BE UNINTERRUPTED OR ERROR-FREE. 
! THE END-USER UNDERSTANDS THAT THE PROGRAM WAS DEVELOPED FOR 
! RESEARCH PURPOSES AND IS ADVISED NOT TO RELY EXCLUSIVELY ON THE 
! PROGRAM FOR ANY REASON. 

! IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY 
! PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL 
! DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS 
! SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA 
! HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY 
! OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT 
! NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
! FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON 
! AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS
! TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
! MODIFICATIONS. 


! !INTRODUCTION:

! \begin{verbatim}
! MPH Version 3 combines all features of MPH version 1, unifies the
! interfaces, and provides more flexible components integration/execution 
! modes.  
!
! In a distributed multi-component environment, each executable resides
! on a set of SMP nodes. Components within an executable may overlap on
! different nodes or processors. 
!
! MPH Version 3 contains the following functionality: 
! 
!   o component name registration 
!   o resource allocation 
!   o multi-component single executable, multi-component 
!     multi-executable, etc. 
!   o inter-component communication 
!   o inquiry on the multi-component environment 
!   o standard in/out redirect 
!
! Please see more information at 
!    http://www.nersc.gov/research/SCG/acpi/MPH
! and please list the following in your reference if useful:      
!      "MPH: a Library for Distributed Multi-Component Environment" 
!       Chris Ding and Yun He, Lawrence Berkeley Nat'l Lab Tech 
!       Report 47930, May 2001. 
!
! Consider the entire simulation system (CCSM) consists of many
! executables, each executable containing one or more components.  
! This architecture offers complete flexibility, and is consistent
! with CORBA, DCE, CCA et al.
!
! 1) Every executable starts with
!
!    mpi_exec_world =                       &
!          MPH_components(name1='ocean', name2='atmosphere',...)
!
!    You may have only one component in this executable, or up to
!    10 components in this executable. Component names are nametags
!    (place holder) and are completely arbitrary.  They must be 
!    self-consistently used, and match the "mph_processors_map.in" 
!    registration file.  This setup subroutine replaces MPH_setup()
!    in MPH version 1.  All other MPH functionality remains identical.
!
! 2) Some usages:
!
!    a) CCSM example.  Ice & Land share one executable.
! 
!           coupler - one executable 
!        atmosphere - one executable 
!             ocean - one executable 
!        ice & land - one executable with 2 components, 
!                     they may overlap on processors 
!    
!    b) CCSM example.  Multiple instances of atmosphere.
!
!           coupler - one executable 
!        atmosphere - one executable of 3 components 
!                     each is a CCM instance of a different Dycore. 
!             ocean - one executable 
!              land - one executable with 3 components for CCMs.
!                     each is a land model to match the CCM 
!               ice - one executable 
!  
!   c) PCM example. 
!
!                 couple - one executable 
!      atmosphere & land - one executable 
!            ocean & ice - one executable 
!
! 3) "mph_processors_map.in" registration file
!
!    The following example contains 3 executables:
!      1st executable has a single component: coupler 
!      2nd executable has 2 components: ocean, ice 
!      3rd executable has 3 components: atmosphere, land, chemistry 
!   
!        BEGIN
!        coupler
!        Multi_Comp_Start
!        2
!        ocean      0  3
!        ice        4 10
!        Multi_Comp_End
!        Multi_Comp_Start
!        3
!        atmosphere 0  10
!        land       11 13
!        chemistry  14 25
!        Multi_Comp_End
!        END
!   
!   a) Allocation of processors for each executable is controlled 
!      by job launching process (different on IBM, SGI, Compaq).
!
!   b) Processor ranges for each components are defined local 
!      to the exeutable.
! \end{verbatim}

! \section{How to Use}

! \begin{verbatim}
! Users need to "use MPH_module" in the application codes, and 
! invoke the appropriate MPH_components function for the multiple 
! components in each executable. For example, ICE_LAND_World = 
! MPH_components (name1="ice", name2="land"). You could use 
! MPH_debug call to determine the output message amount, the 
! default level is 0. "MPH_help" call provides you the available 
! inquiry functions for that mode. 
!
! An input file called "mph_processors_map.in" to give detailed 
! information of component nametags and processor ranges. See more 
! detail about how this file looks like in Section 1.

! Each component maintains its own output in a separate file (file 
! name defined by environment variable either in command line or 
! in batch run script), assuming the local processor 0 of each 
! component being responsible for most output, other occasional 
! writes from all the components are stored in one combined standand 
! output file. 
!
! This is accomplished by processor rank 0 of each component call 
! subroutine "MPH_redirect_output" with the model name as argument. 
! IBM and SGI could do the output redirect with the help of system 
! function "getenv" or "pxfgetenv". Compaq cannot do this.  And 
! T3E is able to get the correct output files created using 
! "pxfgetenv", but only output with those "write(6,*)" could be 
! redirected, but not those with "write(*,*)", since * is equal to
! unit 101, and permanently related to the non-redirectable stdout. 
! \end{verbatim}

! \section{How to Compile and Run}

! \begin{verbatim}
! 
! The shared "Makefile" detects the machine architecture and 
! compiles appropriately for IBM, SGI and Compaq. For test case 1,
! type "make test1", and for test case 2, type "make test2". 
! or "gmake ..." depends on your machine). 
!
! After compile, you will have executables generated ("ice_land",
! "cpl", "pop_atm" for test1, and "ice_land", "cpl" for test2) in
! the corresponding subdirectory. Each sample subdirectoy also 
! includes batch scripts and sample output. 
!
! Go to that directory first (here we use test2 as an example), 
! and then: 
!
! 1) To run on NERSC and NCAR IBM SP interactively: 
!    a) % unsetenv MP_TASKS_PER_NODE 
!    b) % setenv ice_out_env ice.log 
!       % setenv land_out_env land.log 
!       % setenv cpl_out_env cpl.log 
!    c) Make sure the following command in ONE LINE: 
!       % poe -pgmmodel mpmd -cmdfimle tasklist -nodes 3 -procs 6 
!             -stdoutmode ordered -infolevel 2 > & output & 
!      
!    This is to run the executables listed in user supplied "tasklist"
!    in the mpmd mode on total of 3 nodes and 6 procs.   
!    
!    And "tasklist" looks like this:
!          ice_land
!          ice_land
!          cpl
!          cpl
!          ice_land
!          ice_land
!            
!    To run on IBM SP with batch script: 
!    % llsubmit runscript.ibm
 
!    And "runscript.ibmc" looks like this: 
!        #! /usr/bin/csh -f 
!        # @ output  = poe.stdout.$(jobid).$(stepid)
!        # @ error  = poe.stderr.$(jobid).$(stepid)
!        # @ class = debug
!        # @ job_type = parallel
!        # @ task_geometry = {(0,2)(1,3)(4,5)}
!        # @ total_tasks=6
!        # @ network.MPI = css0, not_shared, us 
!        # @ queue
!        setenv MP_PGMMODEL       mpmd  
!        setenv MP_CMDFILE        tasklist
!        setenv ice_out_env ice.log
!        setenv land_out_env land.log
!        setenv cpl_out_env cpl.log
!        poe  

!    Again, it needs a user supplied "tasklist", and it runs in mpmd mode.
!    The task_geometry keyword specifies which tasks run in the same node.

! 2) We could not run it on NERSC CRAY T3E since there is no mpmd
!    mechnism. 
!
! 3) To run on NCAR SGI interactively: 
!    a) % setenv ice_out_env ice.log 
!       % setenv land_out_env land.log 
!       % setenv cpl_out_env cpl.log 
!    b) % mpirun -p "[%g]" -np 4 ice_land : -np 2 cpl > output 

!    This is to run ice_land on 4 procs and cpl on 2 procs.   
!    [%g] is to print the global id as a prefix for each output line.  
!
! 4) To run on NCAR Compaq with batch script:    
!    % prun -n6 -t runscript.dec
!
!    And "runscript.dec" looks like this: 
!          #! /bin/csh 
!          if ($RMS_RANK >= 0 && $RMS_RANK <= 3) ice-land &                        
!          if ($RMS_RANK >= 4 && $RMS_RANK <= 5) cpl &                      
!          exit         
!\end{verbatim}

!\section{Acknowledgement}
! MPH is developed in collaboration with Tony Craig, Brian Kauffman,
! Vince Wayland and Tom Bettge of National Center of Atmospheric Research, 
! and Rob Jacobs and Jay Larson of Argonne National Laborotory. 
! Vince Wayland of NCAR contributes for Makefiles on SGI and Compaq.
! Dan Anderson and Bill Celmaster of NCAR also help in the batch run
! script of Compaq. 
! This work is supported by the Office of Biological and
! Environmental Research, Climate Change Prediction Program, under 
! ACPI Avant Garde project, and by the Office of Computational and 
! Technology Research, Division of Mathematical, Information, and 
! Computational Sciences, both of the U.S. Department of
! Energy under contract number DE-AC03-76SF00098. 
!EOI
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!BOP
!
! !MODULE: MPH_module -- Multi Program-Components Handshaking 
!
! !DESCRIPTION:
!     This module multiple executables with multiple components in
!     each executable.  This module multiple executables with multiple
!     components in.

! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      module MPH_module

! !USES:
      implicit none
      include 'mpif.h'
      private    ! except
  
! !PUBLIC MEMBER FUNCTIONS:
      public :: MPH_components
      public :: PE_in_component
      public :: PE_in_num_comps
      public :: MPH_global_id
      public :: MPH_comm_join
!     public :: MPH_redirect_output
      public :: MPH_help
      public :: MPH_debug
      public :: MPH_timer
      public :: MPH_total_components
      public :: MPH_comp_name
      public :: MPH_comp_id
      public :: MPH_local_world
      public :: MPH_exe_id
      public :: MPH_total_num_exe
      public :: MPH_num_comps
      public :: MPH_local_proc_id
      public :: MPH_local_totProcs
      public :: MPH_global_proc_id
      public :: MPH_global_totProcs
      public :: MPH_exe_world
      public :: MPH_exe_low_proc_limit
      public :: MPH_exe_up_proc_limit

! !PUBLIC DATA MEMBERS:
      integer, public :: istatus(MPI_STATUS_SIZE), ierr  
      integer, public :: MPH_Global_World  ! total processor for the whole world

! !DEFINED PARAMETERS:
      integer, parameter :: max_num_comps=10     ! maximum number of components
      integer, parameter :: maxProcs_comp=128    ! maximum number of procs per comp
      integer, parameter :: max_num_exes=10      ! maximum number of executables    
      integer, parameter :: N_CHANNELS=10        ! number of channels for timing

! !LOCAL VARIABLES:
      type Acomponent
         character (len=80) :: name         ! component name
         integer :: num_process             ! number of processors 
         integer :: process_list (maxProcs_comp) 
                                    ! global processor_id, increasing order
      end type Acomponent    
  
      type (Acomponent) :: components (max_num_comps) ! allocate components
      integer :: MPI_Acomponent
  
      integer :: local_world (max_num_comps)    ! communicator for each component
      integer :: local_proc_id (max_num_comps)  ! proc id in each component
      integer :: local_totProcs (max_num_comps) 
                              ! total number of processors in each component
      integer :: global_proc_id     ! proc id in the whole world
      integer :: global_totProcs    ! total number of processors

      integer :: COMM_master    ! communicator for submaster of each component
      integer :: total_components  ! total number of components
      character (len=80) :: component_names (max_num_comps)  ! component names
      character (len=80) :: name (max_num_comps)  ! name array used in setup
      integer comp_id (max_num_comps)        ! component id of each component

      integer :: num_comps (max_num_comps)  
                                  ! number of components in each executable
      integer :: exe_low_proc_limit (max_num_comps)
                                  ! lower processor limit of each component
                                  ! in each executable world
      integer :: exe_up_proc_limit (max_num_comps)
                                  ! upper processor limit of each component
                                  ! in each executable world
      integer :: exe_world_proc_id (max_num_exes)
                                  ! processor id in the executable world
      integer :: exe_world_totProcs (max_num_exes) 
                                  ! number of processors in each executable
      integer :: exe_world (max_num_exes)  
                                  ! communicator for each executable 
      integer :: exe_ids (max_num_comps)  ! executable ids
      integer :: total_num_exe    ! total number of executables
      integer :: exe_id           ! executable id

      integer :: debug_level = 0  ! level of debug
!     .. for timer ..
      real (kind=8) :: init_time = -1.0
      real (kind=8) :: last_time, tot_time (0:N_CHANNELS)

!EOP
!-------------------------------------------------------------------------

      contains 
!==========================================================================
!      integer function MPH_components(name1, name2, name3, name4, name5)
!==========================================================================
!BOP 
!
! !IROUTINE: MPH_components -- main MPH setup function
!
! !DESCRIPTION:
!    This is the main function for each of the executable to call to 
!    setup the distributed multi-component environment.  For example,
!    if ocean and atmosphere sits in one executable, the source code
!    will contain:
!    \begin{verbatim}
!    mpi_exec_world =                       &
!          MPH_components(name1='ocean', name2='atmosphere',...)
!    \end{verbatim}
!    This function returns the MPI communicator of local executable world. 

! !REVISION HISTORY:
!     2001-Dec-03 -- reduce from 10 arguments to 5
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_components(name1, name2, name3, name4, name5)

! !USES:
      implicit none

! !INPUT PARAMETERS:
!     These are component names
      character(len=*),intent(in) :: name1
      character(len=*),intent(in),optional :: name2, name3, name4, name5

! !OUTPUT PARAMETERS:
!     ! This function returns the MPI communicator of this executable 
!     ! world: exe_world (exe_id)

! !SEE ALSO: 
!     MPH_init, MPH_local, MPH_global, MPH_debug

! !LOCAL VARIABLES:
      integer :: k

!EOP 
!----------------------------------------------------------------------

      name(1) = trim(name1)
      do k = 2, max_num_comps
         name(k) = ''
      enddo
      if (present(name5)) name(5) = trim(name5)
      if (present(name4)) name(4) = trim(name4)
      if (present(name3)) name(3) = trim(name3)
      if (present(name2)) name(2) = trim(name2)

      write(*,*)'==================================='
      if (debug_level == 2) then
         do k = 1, max_num_comps
            write(*,*)'name(',k,')=', name(k)
         enddo
      endif

      call MPH_init ()
      call MPH_local ()
      call MPH_global ()

      MPH_components = exe_world (exe_id)
      end function MPH_components


!=========================================================================
!      subroutine MPH_init ()
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_init -- initialize MPI and read the processors map info
!
! !DESCRIPTION:
!    This routine calls mpi\_init, obtains global processor id.  It
!    reads and processes the ``processors\_map.in'' file.  It also defines 
!    an MPI\_Acomponent sturcture (includes component name, number of
!    processors and processor list for each component) for easy 
!    gather and scatter.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      subroutine MPH_init ()

! !USES:
      implicit none

! !SEE ALSO:
!     MPH_read_list, MPH_local, MPH_global

! !LOCAL VARIABLES:      
      integer :: iblock(3), idisp(3), itype(3)

!EOP
!---------------------------------------------------------------------

      call MPI_INIT (ierr)
      call MPI_COMM_DUP (MPI_COMM_WORLD, MPH_Global_World, ierr)
      call MPI_COMM_RANK (MPH_Global_World, global_proc_id, ierr)
      call MPI_COMM_SIZE (MPH_Global_World, global_totProcs, ierr)

      total_components = MPH_read_list("mph_processors_map.in",
     &        "PROCESSORS_MAP", component_names,   
     &        exe_low_proc_limit, exe_up_proc_limit,
     &        local_totProcs, exe_ids, num_comps, total_num_exe)       

! create a new MPI data type MPI_Acomponent

      iblock(1) = 80
      iblock(2) = 1
      iblock(3) = maxProcs_comp
      idisp(1) = 0
      idisp(2) = 80
      idisp(3) = 84
      itype(1) = MPI_CHARACTER
      itype(2) = MPI_INTEGER
      itype(3) = MPI_INTEGER
      call MPI_TYPE_STRUCT(3,iblock,idisp,itype,MPI_Acomponent,ierr)
      call MPI_TYPE_COMMIT (MPI_Acomponent, ierr)

      end subroutine MPH_init


!=========================================================================
!      subroutine MPH_local ()
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_local -- local handshaking
!
! !DESCRIPTION:
!    This routine first defines exe\_id , and creates local exe\_world 
! for each executable.  It then gathers gobal processor ids onto
! submaster (whose rank is 0 in exe\_world).  And then it creates 
! local\_world for each component within exe\_world based on its
! upper and lower processor limits.  Finally it collects name, number
! of processors, and processor list of each component onto submaster of 
! each executable world.

! !REVISION HISTORY:
!     2001-Dec-13 -- add warning for overlapping processors
!     2001-Nov-27 -- add local_totProcs for single component executables
!     2001-Nov-19 -- add PROTEX convention, use new MPH_read_list interface 
!     2001-May-20 -- first prototype

! !INTERFACE: 
      subroutine MPH_local ()
      
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_init, MPH_global, MPH_find_name
      
! !LOCAL VARIABLES:      
      integer :: color, key
      integer :: id, comp_id_end, i, k

!EOP
!---------------------------------------------------------------------

! define exe_id, create exe_world for each executable
      comp_id(1) = MPH_find_name (name(1), component_names,
     &                              total_components)
      key = 0
      exe_id = exe_ids (comp_id(1))
      call MPI_COMM_SPLIT (MPH_Global_World, exe_id, key,
     &                     exe_world(exe_id),ierr)

! setup exe_world_proc_id, exe_world_totProcs
      call MPI_COMM_RANK(exe_world(exe_id), exe_world_proc_id(exe_id),
     &                   ierr)
      call MPI_COMM_SIZE(exe_world(exe_id), exe_world_totProcs(exe_id),
     &                   ierr)

! initialize process_list to be -1 
      do id = 1, total_components
         do i = 1, exe_world_totProcs(exe_id)
           components(id)%process_list(i) = -1
         enddo
      enddo

      components(comp_id(1))%name = name(1)
      components(comp_id(1))%num_process=exe_world_totProcs(exe_id)
      
! gather processor ids to 0th proc in this component
      call MPI_GATHER (global_proc_id, 1, MPI_INTEGER,
     &                 components(comp_id(1))%process_list, 1, 
     &                 MPI_INTEGER, 0, exe_world(exe_id), ierr)

      if (name(2) .ne. '') then
         ! more than one component in this executable
         comp_id(2) = MPH_find_name (name(2), component_names,
     &                              total_components)
         comp_id_end = comp_id(2)
         do k = 3, max_num_comps
            if (name(k) .ne. '') then
               comp_id(k) = MPH_find_name (name(k), component_names,
     &                                 total_components)
               comp_id_end = comp_id(k)
            endif
         enddo

! create local world for each component within exe_world
         do id = comp_id(1), comp_id_end
            if (exe_low_proc_limit(id) .le. exe_world_proc_id(exe_id)
     &             .and. exe_world_proc_id(exe_id) .le.
     &             exe_up_proc_limit(id)) then
               color = 1
            else
               color = 2
            endif
            key = 0
            call MPI_COMM_SPLIT (exe_world(exe_id), color, key,
     &                           local_world(id), ierr)
            call MPI_COMM_RANK(local_world(id),local_proc_id(id),ierr)
         enddo
      else
         ! only one component in this executable
         call MPI_COMM_DUP(exe_world(exe_id), local_world(comp_id(1)),
     &                     ierr)
         call MPI_COMM_RANK(local_world(comp_id(1)),
     &                      local_proc_id(comp_id(1)),ierr)
         call MPI_COMM_SIZE(local_world(comp_id(1)),
     &                      local_totProcs(comp_id(1)),ierr)
         exe_low_proc_limit (exe_id) = 0
         exe_up_proc_limit (exe_id) = local_totProcs(comp_id(1))-1  
      endif

! collect name, number of processors, and processor list of each component
! onto submaster of each executable world
      if (exe_world_proc_id(exe_id)==0) then
         if (name(2) .ne. '') then
            components(comp_id(2))%name = name(2)
            components(comp_id(2))%num_process =
     &                             local_totProcs(comp_id(2))

            do i = exe_low_proc_limit(comp_id(2)),
     &             exe_up_proc_limit(comp_id(2))
               components(comp_id(2))%process_list(i-
     &                    exe_low_proc_limit(comp_id(2))+1)
     &       = components(comp_id(1))%process_list(i+1)
            enddo

            do k = 3, max_num_comps
               if (name(k) .ne. '') then
                  components(comp_id(k))%name = name(k)
                  components(comp_id(k))%num_process =
     &                                   local_totProcs(comp_id(k))

                  do i = exe_low_proc_limit(comp_id(k)),
     &                   exe_up_proc_limit(comp_id(k))
                     components(comp_id(k))%process_list(i-
     &                          exe_low_proc_limit(comp_id(k))+1)
     &               = components(comp_id(1))%process_list(i+1)
                  enddo
               endif
            enddo

            components(comp_id(1))%num_process =
     &                             local_totProcs(comp_id(1))

            do i = exe_low_proc_limit(comp_id(1)),
     &             exe_up_proc_limit(comp_id(1))
 
               components(comp_id(1))%process_list(i-
     &                    exe_low_proc_limit(comp_id(1))+1)
     &         = components(comp_id(1))%process_list(i+1)
            enddo
            do i = exe_up_proc_limit(comp_id(1))+1,
     &             exe_world_totProcs(exe_id)
               components(comp_id(1))%process_list(i+1)=-1
            enddo

! check overlap in this multi-component executable
           do i = 1, num_comps(exe_id)
           do k = i+1, num_comps(exe_id)
              if (exe_up_proc_limit(comp_id(i)) >=
     &            exe_low_proc_limit(comp_id(k))) then
                 write(*,*)'warning: local procs ',
     & exe_low_proc_limit(comp_id(k)), ' to ',
     & exe_up_proc_limit(comp_id(i)), ' is overlapped in ',
     & 'these two components of the executable: ',
     & components(comp_id(i))%name, components(comp_id(k))%name 
              endif
           enddo
           enddo
        endif

      endif

      end subroutine MPH_local


!=========================================================================
!      subroutine MPH_global ()
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_global -- global handshaking
!
! !DESCRIPTION:
! This routine first creates an MPI communicator COMM\_master for all
! submasters (whose rank is 0 in the executable world).  It then does
! an MPI\_allgatherv in COMM\_master to collect all the components 
! information from each submaster.  Then each submaster broadcasts 
! AComponents to all PEs in its local exe\_world.  Finally, every 
! processor lists the complete info of all the components.
!
! !REVISION HISTORY:
!     2002-Apr-15 -- correct a bug in declaring sendbuf
!     2001-Dec-13 -- add warning for overlapping processors
!     2001-Nov-19 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      subroutine MPH_global ()
      
! !USES:      
      implicit none
      
! !SEE ALSO:
!     MPH_init, MPH_local, MPH_debug

! !LOCAL VARIABLES:    
      integer :: id, i, color, key
      type (Acomponent) :: sendbuf(max_num_comps)
      integer :: sendcount
      integer :: recvcounts(0:total_num_exe-1)
      integer :: displs(0:total_num_exe-1)

!EOP
!-------------------------------------------------------------------

! create a MPI communicator COMM_master for all submasters
! arrange the rank of the submasters in COMM_master by their exe_id
! i.e., their rank of the executables in "mph_processors_map.in"
      if (exe_world_proc_id (exe_id) == 0) then
         color = 1
      else
         color = 2
      endif
      key = exe_id
      call MPI_COMM_SPLIT (MPH_Global_World,color,key,COMM_master,ierr)

! do a MPH_ALLGATHERV in COMM_master

      if (exe_world_proc_id (exe_id) == 0) then
         sendcount = num_comps(exe_id)
         displs(0) = 0
         recvcounts(0) = num_comps(1)
         do id = 2, total_num_exe
            recvcounts(id-1) = num_comps(id)
            displs(id-1) = displs(id-2) + num_comps(id-1)
         enddo
         do i = 1, num_comps(exe_id)
            sendbuf(i) = components(displs(exe_id-1)+i)
         enddo
      
         call MPI_ALLGATHERV(sendbuf, sendcount, MPI_AComponent,
     &          components, recvcounts, displs, MPI_AComponent,
     &          COMM_master, ierr)
      endif

! submaster broadcast AComponents to all PEs in local exe_world
      call MPI_BCAST (components, total_components,
     &                MPI_Acomponent, 0, exe_world(exe_id), ierr)
             
! everybody lists the complete info
      if (debug_level >= 1) then
         write(*,*)'global_proc_id=', global_proc_id
      endif

      if (debug_level == 2) then
         write(*,*)'infos I have for all proc of all components are:'
         do i = 1, total_components
            write(*,*)'   ', components(i)%name
            write(*,*)'   ', components(i)%num_process
            write(*,*)'   ', components(i)%process_list(1:
     & components(i)%num_process)
         enddo
      endif

! print a warning message for overlapping processors 
      if (PE_in_num_comps() > 1) then
         write(*,*)'warning: global processor ', global_proc_id, '
     & is overlapped'
      endif 

      end subroutine MPH_global


!=========================================================================
!      logical function PE_in_component (name, comm)
!=========================================================================

!BOP 
!
! !IROUTINE: PE_in_component -- check if a processor is in a component
!
! !DESCRIPTION:
!    This is a logical function to check if a processor is in a component.

! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      logical function PE_in_component (name, comm)
      
! !USES:      
      implicit none

! !SEE ALSO:      
!     MPH_find_name

! !INPUT PARAMETERS:
      character(len=*), intent(in) :: name   ! component name
     
! !OUTPUT PARAMETERS:
!     ! the local communicator of that component is written in comm.
      integer, intent(out) :: comm   ! communicator for the component
     
! !LOCAL VARIABLES:          
      integer :: id, i
 
!EOP
!-------------------------------------------------------------------

      id = MPH_find_name(name, component_names, total_components)
      do i = 1, components(id)%num_process
         if (global_proc_id == components(id)%process_list(i)) then
            comm = local_world (id)
            PE_in_component = .true.
            return
         else
            continue
         endif
      enddo
      
      PE_in_component = .false.

      end function PE_in_component 


!=========================================================================
!      integer function PE_in_num_comps ()
!=========================================================================

!BOP 
!
! !IROUTINE: PE_in_num_comps -- return the number of components a processor
! in
!
! !DESCRIPTION:
!    This function returns the number of components a processor is in.

! !REVISION HISTORY:
!     2001-Dec-13 -- first prototype

! !INTERFACE: 
      integer function PE_in_num_comps ()

! !USES:      
      implicit none

! !SEE ALSO:      
!     PE_in_component

! !LOCAL VARIABLES:          
      integer :: id, i
 
!EOP
!-------------------------------------------------------------------

      PE_in_num_comps = 0
      do id = 1, total_components
         do i = 1, components(id)%num_process
            if (global_proc_id == components(id)%process_list(i)) then
               PE_in_num_comps = PE_in_num_comps + 1
            endif
         enddo
      enddo
      
      end function PE_in_num_comps 


!=========================================================================
!      integer function MPH_global_id (cname, lid) 
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_global_id -- find global processor id
!
! !DESCRIPTION:
!    This function returns global processor id given the component 
!    name and local processor id in that component.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_global_id (cname, lid)
      
! !USES:      
      implicit none

! !SEE ALSO:      
!     MPH_find_name        
      
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: cname   ! component name
      integer, intent(in) :: lid  
                               ! local processor id in the component
      
! !OUTPUT PARAMETERS:
!     ! This function returns global_proc_id given the component 
!     ! name and local_proc_id in that component.
           
! !LOCAL VARIABLES:         
      integer :: temp        
      
!EOP
!-------------------------------------------------------------------

! then find out the component rank
      temp = MPH_find_name (cname,component_names,total_components)

! process_list starts from 1, while proc rank starts from 0
      MPH_global_id = components(temp) % process_list(lid+1)

      end function MPH_global_id


!=========================================================================
!      subroutine MPH_comm_join (name1, name2, comm_joined) 
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_comm_join -- join two components
!
! !DESCRIPTION:
!    This routine creates a joined MPI communicators for any two 
!    components.  The order of these two components appeared in 
!    the subroutine parameter argumment has an effect on the 
!    local process id within the joined communicator.

! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      subroutine MPH_comm_join (name1, name2, comm_joined)
      
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_find_name, PE_in_component

! !INPUT PARAMETERS:
      character(len=*), intent(in) :: name1, name2   ! two component names

! !OUTPUT PARAMETERS:
      integer, intent(out) :: comm_joined   
                                 ! joined communicator for two components
                   
! !LOCAL VARIABLES: 
      integer :: id1, id2      
      integer :: color, key
      logical :: con1, con2
      integer :: comm1, comm2  
 
      
!EOP
!-------------------------------------------------------------------

      id1 = MPH_find_name (name1,component_names,total_components)
      id2 = MPH_find_name (name2,component_names,total_components)

! the order of two components does matter: first one has lower ranks  
! in the new joined communicator, and second one has higher ranks.

      con1 = PE_in_component (name1, comm1)
      con2 = PE_in_component (name2, comm2)
 
      if (con1 .or. con2) then
         color = 1
         if (con1) then 
            key = local_proc_id (id1)
         else
            key = global_totProcs + local_proc_id (id2)
         endif
      else
         color = 2
         key = 0
      endif

      call MPI_COMM_SPLIT (MPH_Global_World,color,key,comm_joined,ierr)

      end subroutine MPH_comm_join


!=========================================================================
!      subroutine MPH_redirect_output (name) 
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_redirect_output -- redirect output from each component
!
! !DESCRIPTION:
!    This routine redirects output to a log file defined by an environment
!    variable.  System functions ("getenv" for IBM and "pxfgetenv" for
!    SGI and T3E) are used to retrieve the environmet variable.
!
! !REMARKS:
!    In order to redirect component output to a separate file, 
!    a user will setup something like the following in the run script:
!         setenv ice_out_env ice.log
!         setenv land_out_env land.log
!         setenv cpl_out_env cpl.log 

! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      subroutine MPH_redirect_output (name)
      
! !USES:      
      implicit none
      
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: name   ! part of the log file name

! !LOCAL VARIABLES: 
      integer :: lenname, lenval, rcode
      character(len=80) :: output_name_env
      character(len=80) :: output_name, temp_value
       
!EOP
!-------------------------------------------------------------------

      output_name = ' '
      output_name_env = trim (name) // "_out_env"

#if (defined AIX)
      call getenv (trim(output_name_env), temp_value)
      output_name = trim (temp_value)
      if (len_trim(output_name) == 0) then
         write(*,*)'output file names not preset by env varibales'
         write(*,*)'so output not redirected'
      else
         open (unit=6, file=output_name, position='append')
         call flush_(6)
      endif
#endif

#if (defined IRIX64 || defined CRAY || defined sn6711)
      lenname = len_trim (output_name_env)
      call pxfgetenv (output_name_env,lenname,output_name,lenval,rcode)
      if (len_trim(output_name) == 0) then
         write(*,*)'output file names not preset by env varibales'
         write(*,*)'so output not redirected'
      else
         open (unit=6, file=output_name, position='append')
         call flush(6)
      endif
#endif

#if (!defined AIX && !defined IRIX64 && !defined CRAY && !defined sn6711)
      write(*,*) 'No implementation for this architecture'
      write(*,*) 'output redirect is not performed by getenv'
#endif

      end subroutine MPH_redirect_output


!=========================================================================
!      integer function MPH_read_list (filename, filetag, namelist,  
!     &           low, up, local_num, id_exe, num_comp, total_exe)        
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_read_list -- read and process info from "mph_processors_map.in"
!
! !DESCRIPTION:
!    This routine reads and processes info from "processors\_map.in".   
!    Please see a sample input file in Introduction.  

! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!                    remove two arguments: max_num_comp, max_num_exe
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_read_list (filename, filetag, namelist,  
     &           low, up, local_num, id_exe, num_comp, total_exe)       
      
! !USES:      +
      implicit none

! !SEE ALSO:
!     MPH_find_name, MPH_init
       
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: filename    ! the input file name
      character(len=*), intent(in) :: filetag     ! "PROCESSORS_MAP"
     
! !OUTPUT PARAMETERS:
      character (len=80), intent(out) :: namelist(max_num_comps)
                                  ! component names
      integer, intent(out) :: low(max_num_comps)
                                  ! lower processor limit of each component
                                  ! in each executable world
      integer, intent(out) :: up(max_num_comps)
                                  ! upper processor limit of each component
                                  ! in each executable world
      integer, intent(out) :: local_num(max_num_comps)
                           ! total number of processors for each component 
      integer, intent(out) :: id_exe(max_num_comps)   ! executable ids
      integer, intent(out) :: num_comp(max_num_comps) 
                                 ! number of components in each executable
      integer, intent(out) :: total_exe   ! total number of executables
        
! !LOCAL VARIABLES:    
      integer :: i, k, multi_num, id    
      character (len=80) :: firstline, temp
      integer :: itemp1, itemp2
     
!EOP
!-------------------------------------------------------------------

      open(10, file=filename, status='unknown')
      read(10, *, end=100) firstline
      if (firstline .ne. filetag) then
         print *, 'ERROR: filetag inconsistent', filename
         print *, 'ERROR: ', filetag, '!=', firstline
         stop
      endif

      read(10, *, end=200) temp
      if (temp .ne. "BEGIN") then
         print *, 'ERROR: no BEGIN in ', filename
         stop
      endif

      i = 0
      id = 0
1000  read(10, *, err=300, end=400) temp
      if (temp .eq. "Multi_Comp_Start") goto 1100
      if (temp .eq. "END") goto 500
      i = i + 1
      id = id + 1
      if (i . gt. max_num_comps) goto 600
      namelist(i) = temp
      low(i) = -1
      up(i) = -1
      local_num(i) = -1
      id_exe(i) = id  
      num_comp(id) = 1  
      goto 1000

1100  read(10, *, err=310, end=400) multi_num
         id = id + 1
         if (id . gt. max_num_exes) goto 650
      do k = 1, multi_num
         read(10, *, err=320, end=400) temp, itemp1, itemp2
         i = i + 1
         if (i . gt. max_num_comps) goto 600
         namelist(i) = temp
         low(i) = itemp1
         up(i) = itemp2
         local_num(i) = itemp2 - itemp1 + 1
         id_exe(i) = id
         num_comp(id) = multi_num
      enddo
      read(10, *, err=300, end=400) temp
      if (temp .ne. "Multi_Comp_End") goto 700
      goto 1000
  
100   print *, 'ERROR: no filetag in ', filename
      stop

200   print *, 'ERROR: no BEGIN in ', filename
      stop

300   print *, 'ERROR: not a character'
      stop

310   print *, 'ERROR: not an integer: number of comps'
      stop

320   print *, 'ERROR: not in format: char, int, int'
      stop

400   print *, 'ERROR: no END in ', filename
      stop

600   print *, 'ERROR: exceed maximum number of components allowed'
      stop

650   print *, 'ERROR: exceed maximum number of executables allowed'
      stop

700   print *, 'ERROD: no Multicompend for corresponding Multicompstart'
      stop

500   total_exe = id
      MPH_read_list = i 

      return
      end function MPH_read_list

 
!=========================================================================
!      integer function MPH_find_name (name, namelist, num)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_find_name -- find name in a namelist
!
! !DESCRIPTION:
!    This routine finds if a certain name exists in an array of
!    namelist and returns the rank if it does or -1 if it does not.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_find_name (name, namelist, num)
      
! !USES:      
      implicit none

! !INPUT PARAMETERS:
      character(len=*), intent(in) :: name   ! name to be found
      integer :: num                      ! length of namelist array
      character (len=80) namelist(num)    ! name list array

! !OUTPUT PARAMETERS:
!     ! the rank of a name in an array of namelist or -1 if not exist
      
! !LOCAL VARIABLES:          
      integer :: i
    
!EOP
!-------------------------------------------------------------------

      do i = 1, num
         if (name == namelist(i)) then
!            print *, i, name, namelist(i)
            goto 100
         endif
      enddo

! name is not found
      MPH_find_name = -1
      print *, "ERROR: ", name, " not found in mph_processors_map.in"
      stop

100   MPH_find_name = i
      return
      end function MPH_find_name 


!=========================================================================
!      subroutine MPH_help (arg)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_help -- display help info
!
! !DESCRIPTION:
!    This routine displays some help info for the MPH setup interface
!    and some inquiry functions.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      subroutine MPH_help (arg)
     
! !USES:      
      implicit none

! !INPUT PARAMETERS:
      character(len=*), intent(in) :: arg   ! either 'on' or 'off'
     
!EOP
!-------------------------------------------------------------------

      write(*,*)'Message from MPH_help:'
      if (arg .eq. 'off') then
         write(*,*)'off'
         return
      else if (arg .eq. 'on') then
         write(*,*)'default debug level=0, choose 1 or 2 for more
     & debug messages' 
         write(*,*)'MPH: multiple executables with multiple components
     & in each executable'
         write(*,*)'Sample setup functions are: '
         write(*,*)'   ICE_LAND_World = MPH_components (name1="ice",
     & name2="land")'
         write(*,*)'   CPL_World = MPH_components ("name1=coupler")'
         write(*,*)'Required input file is "mph_processors_map.in"'
         write(*,*)'Subroutine call to join two communicators is:'
         write(*,*)'   MPH_comm_join (name1,name2,comm_joined)'

         write(*,*)'Available inquiry functions are:'
         write(*,*)'(arguments are optional except for MPH_global_id)'
         write(*,*)'   integer MPH_global_id(cname,lid)'
         write(*,*)'   integer MPH_total_components()'
         write(*,*)'   character (len=80) MPH_comp_name(cid)'
         write(*,*)'   integer MPH_comp_id(cname)'
         write(*,*)'   integer MPH_local_world(cname)'
         write(*,*)'   integer MPH_exe_id(cname)'
         write(*,*)'   integer MPH_total_num_exe()'
         write(*,*)'   integer MPH_num_comps(eid)'
         write(*,*)'   integer MPH_local_proc_id(cid)'
         write(*,*)'   integer MPH_local_totProcs(cid)'
         write(*,*)'   integer MPH_global_proc_id()'
         write(*,*)'   integer MPH_global_totProcs()'
         write(*,*)'   integer MPH_exe_world(eid)'
         write(*,*)'   integer MPH_exe_low_proc_limit(eid)'
         write(*,*)'   integer MPH_exe_up_proc_limit(eid)'
      else
         write(*,*)'wrong argument for MPH_help'
      endif
      write(*,*)

      end subroutine MPH_help


!=========================================================================
!      subroutine MPH_debug (level)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_debug -- define debug level
!
! !DESCRIPTION:
!    This routine defines the debug level. The higher the level is,
!    the more debug information the code will display.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      subroutine MPH_debug (level)

! !USES:      
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: level  !  0 (default), 1 or 2

!EOP
!-------------------------------------------------------------------

      debug_level = level
      end subroutine MPH_debug


!=========================================================================
!      integer function MPH_timer (flag, channel)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_timer -- collect timing info in different channels.
!
! !DESCRIPTION:  This function collects timing info in different
!  channels.
!  \begin{verbatim}
!  Usage:
!  channel 0 is the default channel, using init_time.
!  ---------------------------------------------------------
!  timer calls to walk-clock dclock(), and do the following:
!  ---------------------------------------------------------
!  flag=0  : Sets initial time; init all channels.
!  flag =1 : Calculates the most recent time interval; accrues it 
!            to the specified channel (default 0);  Returns it to 
!            calling process.
!  flag =2 : Calculates the most recent time interval; accrues it 
!            to the specified channel (default 0);  Returns the 
!            curent total time in the specified channel.
!  ---------------------------------------------------------
!  \end{verbatim}

! !REVISION HISTORY:
!     2001-Nov-19 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      real (kind=8) function MPH_timer (flag, channel)
     
! !USES:      
      implicit none

! !INPUT PARAMETERS:
      integer :: flag, channel
           
! !OUTPUT PARAMETERS:
!     ! real (kind=8) MPH_timer (flag, channel)
      
! !LOCAL VARIABLES:    
      real (kind=8) :: new_time, delta_time, MPI_Wtime

!EOP
!-------------------------------------------------------------------
    
      new_time = MPI_Wtime()

      if (flag == 0) then
         init_time = new_time
         last_time = new_time
         tot_time = 0.0
         MPH_timer = new_time - init_time
      else if (init_time == -1.0) then
!        Error Condition
         MPH_timer = init_time
      endif

! Timer is initialized and flag != 0

      delta_time = new_time - last_time
      last_time = new_time

! For channel=0 or other undefined channels which is treated as 0
      if ( channel < 0  .or. channel > N_CHANNELS) then
         write(*,*) 'Timer channel is not properly specified!'
      endif

! channel != 0

      if (flag == 1) then
         tot_time(channel) = tot_time(channel) + delta_time
         MPH_timer = delta_time
      else if (flag == 2) then
         tot_time(channel) = tot_time(channel) + delta_time
         MPH_timer = tot_time(channel)
      else
!        Error Condition
         MPH_timer = -1.0
      endif

      end function MPH_timer


! ---- some inquiry functions----------------------------------------------

!=========================================================================
!      integer function MPH_total_components ()
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_total_components -- find number of total components
!
! !DESCRIPTION:
!    This function returns the number of total components.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_total_components ()

! !USES:      
      implicit none
  
! !OUTPUT PARAMETERS:
!     ! total_components
 
!EOP
!-------------------------------------------------------------------

      MPH_total_components = total_components
      end function MPH_total_components


!=========================================================================
!      character (len=80) function MPH_comp_name (cid)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_comp_name -- find component name given component id
!
! !DESCRIPTION:
!    This function returns component name given component id.
!
! !REVISION HISTORY:
!     2001-Dec-13 -- use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      character (len=80) function MPH_comp_name (cid)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_find_name, MPH_comp_id, MPH_comm

! !INPUT PARAMETERS:
      integer, intent(in), optional :: cid   ! component id
      
! !OUTPUT PARAMETERS:
!     ! component_names (cid)

! !LOCAL VARIABLES:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cid)) then 
         MPH_comp_name = component_names (cid)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cid is required but missing for
     & function MPH_comp_name(cid)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_comp_name = component_names(id)
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_comp_name


!=========================================================================
!      integer function MPH_comp_id (cname)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_comp_id -- find component id given component name
! !DESCRIPTION:
!    This routine returns component id given component name.
!
! !REVISION HISTORY:
!     2001-Dec-13 -- use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_comp_id (cname)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_find_name, MPH_comp_name, MPH_comm

! !INPUT PARAMETERS:
      character(len=*), intent(in), optional :: cname   ! component name
     
! !OUTPUT PARAMETERS:
!     ! MPH_comp_id       

! !LOCAL VARIABLES:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cname)) then 
         MPH_comp_id = MPH_find_name (cname, component_names,
     &                                       total_components)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cname is required but missing for
     & function MPH_comp_id(cname)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_comp_id = id
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_comp_id


!=========================================================================
!      integer function MPH_local_world (cname)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_local_world -- find local communicator given component name
! !DESCRIPTION:
!    This routine returns local communicator given component name.
!
! !REVISION HISTORY:
!     2001-Dec-13 -- first prototype

! !INTERFACE: 
      integer function MPH_local_world (cname)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_find_name, MPH_comp_id

! !INPUT PARAMETERS:
      character(len=*), intent(in), optional :: cname   ! component name
     
! !OUTPUT PARAMETERS:
!     ! MPH_local_world       

! !LOCAL VARIABLES:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cname)) then 
         id = MPH_find_name (cname, component_names,
     &                                   total_components)
         MPH_local_world = local_world(id)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cname is required but missing for
     & function MPH_local_world(cname)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_local_world = local_world(id)
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_local_world


!=========================================================================
!      integer function MPH_exe_id (cname)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_exe_id -- find executable id given component name
!
! !DESCRIPTION:
!    This function returns the executable id given component name.
!
! !REVISION HISTORY:
!     2001-Dec-13 -- use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_exe_id (cname)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_find_name

! !INPUT PARAMETERS:
      character(len=*), intent(in), optional :: cname   ! component name
     
! !OUTPUT PARAMETERS:
!     ! exe_ids (id)

! !LOCAL PARAMETERS:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cname)) then 
         id = MPH_find_name (cname, component_names,
     &                                total_components)
         MPH_exe_id = exe_ids (id)
      else if (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cname is required but missing for
     & function MPH_exe_id(cname)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_exe_id = exe_ids (id)
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_exe_id


!=========================================================================
!      integer function MPH_total_num_exe ()
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_total_num_exe -- find total number of executables
!
! !DESCRIPTION:
!    This fuction returns the total number of executables.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_total_num_exe ()
     
! !USES:      
      implicit none
    
! !OUTPUT PARAMETERS:
!     ! total_num_exe

!EOP
!-------------------------------------------------------------------

      MPH_total_num_exe = total_num_exe
      end function MPH_total_num_exe


!=========================================================================
!      integer function MPH_num_comps (eid)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_num_comps -- find number of components in an executable
!
! !DESCRIPTION:
!    This fuction returns number of components in an executable
!    given the executable id.

! !REVISION HISTORY:
!     2001-Dec-13 -- use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_num_comps (eid)
     
! !USES:      
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in), optional :: eid   ! executable id
           
! !OUTPUT PARAMETERS:
!     ! num_comps (eid)
            
!EOP
!-------------------------------------------------------------------

      if (present(eid)) then 
         MPH_num_comps = num_comps (eid)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument eid is required but missing for
     & function MPH_num_comps(eid)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         MPH_num_comps = num_comps (eid)
      endif

100   return
      end function MPH_num_comps


!=========================================================================
!      integer function MPH_local_proc_id (cid)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_local_proc_id -- find local processor id in a component
!
! !DESCRIPTION:
!    This function returns the local processor id given the component id.
!
! !REVISION HISTORY:
!     2001-Dec-13 -- use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_local_proc_id (cid)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_global_proc_id

! !INPUT PARAMETERS:
      integer, intent(in), optional :: cid   ! component id
      
! !OUTPUT PARAMETERS:
!     ! local_proc_id (cid) 

! !LOCAL VARIABLES:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cid)) then 
         MPH_local_proc_id = local_proc_id (cid)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cid is required but missing for
     & function MPH_local_proc_id(cid)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_local_proc_id = local_proc_id (id)
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_local_proc_id


!=========================================================================
!      integer function MPH_local_totProcs (cid)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_local_totProcs -- find total number of processors 
! in a component.
!
! !DESCRIPTION:
!    This function returns the total number of processors in a 
!    component given the component id.
!
! !REVISION HISTORY:
!     2001-Dec-13 -- use optional argument
!     2001-Dec-13 -- use optional argument
!     2001-Nov-27 -- first prototype

! !INTERFACE: 
      integer function MPH_local_totProcs (cid)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_global_totProcs

! !INPUT PARAMETERS:
      integer, intent(in), optional :: cid   ! component id
      
! !OUTPUT PARAMETERS:
!     ! local_totProcs (cid) 

! !LOCAL VARIABLES:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cid)) then 
         MPH_local_totProcs = local_totProcs (cid)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cid is required but missing for
     & function MPH_local_totProcs(cid)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_local_totProcs = local_totProcs (id)
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_local_totProcs


!=========================================================================
!      integer function MPH_global_proc_id ()
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_global_proc_id -- find global processor id
!
! !DESCRIPTION:
!    This function returns the global processor id.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_global_proc_id ()
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_local_proc_id

! !OUTPUT PARAMETERS:
!     ! global_proc_id
      
!EOP
!-------------------------------------------------------------------

      MPH_global_proc_id = global_proc_id
      end function MPH_global_proc_id


!=========================================================================
!      integer function MPH_global_totProcs ()
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_global_proc_id -- find total number of processors

! !DESCRIPTION:
!    This function returns the total number of processors in MPH world.
!
! !REVISION HISTORY:
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_global_totProcs ()
     
! !USES:      
      implicit none

! !OUTPUT PARAMETERS:
!     ! global_totProcs
      
!EOP
!-------------------------------------------------------------------

      MPH_global_totProcs = global_totProcs
      end function MPH_global_totProcs


!=========================================================================
!      integer function MPH_exe_world (eid)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_local_world -- find local communicator of an executable
!
! !DESCRIPTION:
!    This function returns the local MPI communicator of an executable
!    given the executable id. 
!
! !REVISION HISTORY:
!     2001-Dec-13 -- change function name from MPH_local_world to MPH_exe_world,
!                    use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_exe_world (eid)
     
! !USES:      
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in), optional :: eid   ! executable id
     
! !OUTPUT PARAMETERS:
!     ! exe_world (eid) 

!EOP
!-------------------------------------------------------------------

      if (present(eid)) then 
         MPH_exe_world = exe_world (eid)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument eid is required but missing for
     & function MPH_exe_world(eid)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         MPH_exe_world = exe_world (exe_id)
      endif

100   return
      end function MPH_exe_world


!=========================================================================
!      integer function MPH_exe_low_proc_limit (cid)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_exe_low_proc_limit - find lower processor limit of a component
!
! !DESCRIPTION:
!    This function returns the relative lower processor limit of a component 
!    in the executable world.  
!
! !REVISION HISTORY:
!     2002-Jun-20 -- correct the argument from eid to cid
!     2001-Dec-13 -- use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_exe_low_proc_limit (cid)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_exe_up_proc_limit

! !INPUT PARAMETERS:
      integer, intent(in), optional :: cid   ! component id
      
! !OUTPUT PARAMETERS:
!     ! exe_low_proc_limit (cid) 

! !LOCAL VARIABLES:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cid)) then 
         MPH_exe_low_proc_limit = exe_low_proc_limit (cid)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cid is required but missing for
     & function MPH_exe_low_proc_limit (cid)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_exe_low_proc_limit = exe_low_proc_limit (id)
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_exe_low_proc_limit


!=========================================================================
!      integer function MPH_exe_up_proc_limit (cid)
!=========================================================================

!BOP 
!
! !IROUTINE: MPH_exe_up_proc_limit - find upper processor limit of a component
!
! !DESCRIPTION:
!    This function returns the relative upper processor limit of a component 
!    in the executable world.  
!
! !REVISION HISTORY:
!     2002-Jun-20 -- correct the argument from eid to cid
!     2001-Dec-13 -- use optional argument
!     2001-Nov-15 -- add PROTEX convention
!     2001-May-20 -- first prototype

! !INTERFACE: 
      integer function MPH_exe_up_proc_limit (cid)
     
! !USES:      
      implicit none

! !SEE ALSO:
!     MPH_exe_low_proc_limit

! !INPUT PARAMETERS:
      integer, intent(in), optional :: cid   ! component id
      
! !OUTPUT PARAMETERS:
!     ! exe_up_proc_limit (cid) 

! !LOCAL VARIABLES:
      integer :: id, comm

!EOP
!-------------------------------------------------------------------

      if (present(cid)) then 
         MPH_exe_up_proc_limit = exe_up_proc_limit (cid)
      elseif (PE_in_num_comps() > 1) then
         write(*,*)'ERROR: argument cid is required but missing for
     & function MPH_exe_up_proc_limit (cid)'
         write(*,*)'       due to ambiguity caused by processor overlap,
     & quitting...'
         stop
      else
         do id = 1, total_components
            if (PE_in_component(component_names(id), comm)) then
               MPH_exe_up_proc_limit = exe_up_proc_limit (id)
               goto 100
            endif
         enddo
      endif

100   return
      end function MPH_exe_up_proc_limit

      end module MPH_module


