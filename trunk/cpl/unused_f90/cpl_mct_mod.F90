!===============================================================================
! CVS: $Id: cpl_mct_mod.F90,v 1.1 2005/02/14 20:06:32 steder Exp $
! CVS: $Source: /home/cvsroot/steder/pyCPL/cpl/unused_f90/cpl_mct_mod.F90,v $
! CVS: $Name:  $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_mct_mod -- provides a standard API naming convention for MCT code
!
! !DESCRIPTION:
!    This module should be used instead of accessing mct modules directly.  
!    This module:
!    \begin{itemize}
!    \item Renames MCT routines and data types so that all have an mct\_ prefix
!          and so that related data types and routines have related names. 
!    \item Provides easy and uniform access to 
!          all MCT routines and data types that must be accessed.
!    \item Provides a convienient list of 
!          all MCT routines and data types that can be accessed.
!    \item Blocks access to MCT routines that are not needed.
!    \end{itemize}
!
! !REVISION HISTORY:
!     2001-Aug-14 - B. Kauffman - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_mct_mod

! !USES:

   use shr_sys_mod          ! share system routines
   use shr_mpi_mod          ! mpi layer
   use cpl_kind_mod         ! kinds
   use cpl_const_mod        ! constant module

   use m_MCTWorld           ,only: cpl_mct_world_init     => init

   use m_AttrVect           ,only: cpl_mct_aVect          => AttrVect
   use m_AttrVect           ,only: cpl_mct_aVect_init     => init
   use m_AttrVect           ,only: cpl_mct_aVect_clean    => clean
   use m_AttrVect           ,only: cpl_mct_aVect_zero     => zero
   use m_AttrVect           ,only: cpl_mct_aVect_lsize    => lsize
   use m_AttrVect           ,only: cpl_mct_aVect_indexIA  => indexIA
   use m_AttrVect           ,only: cpl_mct_aVect_indexRA  => indexRA
   use m_AttrVect           ,only: cpl_mct_aVect_getIList => getIList
   use m_AttrVect           ,only: cpl_mct_aVect_getRList => getRList
   use m_AttrVect           ,only: cpl_mct_aVect_nIAttr   => nIAttr
   use m_AttrVect           ,only: cpl_mct_aVect_nRAttr   => nRAttr
   use m_AttrVect           ,only: cpl_mct_aVect_copy     => Copy
   use m_AttrVectComms      ,only: cpl_mct_aVect_scatter  => scatter
   use m_AttrVectComms      ,only: cpl_mct_aVect_gather   => gather 
   use m_AttrVectComms      ,only: cpl_mct_aVect_bcast    => bcast  

   use m_Transfer           ,only: cpl_mct_send          => Send
   use m_Transfer           ,only: cpl_mct_recv          => Recv

   use m_GlobalSegMap       ,only: cpl_mct_gsMap         => GlobalSegMap
   use m_GlobalSegMap       ,only: cpl_mct_gsMap_init    => init
   use m_GlobalSegMap       ,only: cpl_mct_gsMap_clean   => clean
   use m_GlobalSegMap       ,only: cpl_mct_gsMap_lsize   => lsize
   use m_GlobalSegMap       ,only: cpl_mct_gsMap_gsize   => gsize

   use m_Rearranger         ,only: cpl_mct_rearr           => Rearranger
   use m_Rearranger         ,only: cpl_mct_rearr_init      => init
   use m_Rearranger         ,only: cpl_mct_rearr_clean     => clean
   use m_Rearranger         ,only: cpl_mct_rearr_rearrange => rearrange

   use m_Router             ,only: cpl_mct_router        => Router
   use m_Router             ,only: cpl_mct_router_init   => init
   use m_Router             ,only: cpl_mct_router_clean  => clean
   use m_SparseMatrixToMaps ,only: cpl_mct_sMat_2XgsMap  => SparseMatrixToXGlobalSegMap
   use m_SparseMatrixToMaps ,only: cpl_mct_sMat_2YgsMap  => SparseMatrixToYGlobalSegMap

   use m_SparseMatrix       ,only: cpl_mct_sMat              => SparseMatrix
   use m_SparseMatrix       ,only: cpl_mct_sMat_Init         => init
   use m_SparseMatrix       ,only: cpl_mct_sMat_Vecinit      => vecinit
   use m_SparseMatrix       ,only: cpl_mct_sMat_Clean        => clean
   use m_SparseMatrix       ,only: cpl_mct_sMat_indexIA      => indexIA
   use m_SparseMatrix       ,only: cpl_mct_sMat_indexRA      => indexRA
   use m_SparseMatrix       ,only: cpl_mct_sMat_lsize        => lsize
   use m_SparseMatrix       ,only: cpl_mct_sMat_nrows        => nRows
   use m_SparseMatrix       ,only: cpl_mct_sMat_ncols        => nCols
   use m_SparseMatrix       ,only: cpl_mct_sMat_SortPermute  => SortPermute
   use m_SparseMatrix       ,only: cpl_mct_sMat_GNumEl       => GlobalNumElements
   use m_SparseMatrixComms  ,only: cpl_mct_sMat_ScatterByRow => ScatterByRow
   use m_SparseMatrixComms  ,only: cpl_mct_sMat_ScatterByCol => ScatterByColumn
   use m_MatAttrVectMul     ,only: cpl_mct_sMat_avMult       => sMatAvMult
   use m_GlobalToLocal      ,only: cpl_mct_sMat_g2lMat       => GlobalToLocalMatrix

   use m_List               ,only: cpl_mct_list          => list     
   use m_List               ,only: cpl_mct_list_init     => init
   use m_List               ,only: cpl_mct_list_get      => get 
   use m_List               ,only: cpl_mct_list_nitem    => nitem 
   use m_List               ,only: cpl_mct_list_clean    => clean
   use m_string             ,only: cpl_mct_string        => string 
   use m_string             ,only: cpl_mct_string_clean  => clean
   use m_string             ,only: cpl_mct_string_toChar => toChar 
   use m_die                ,only: cpl_mct_perr_die      => mp_perr_die

   implicit none

!  This was added specifically for ES compiler bug:
   public :: cpl_mct_list, cpl_mct_gsmap, cpl_mct_router

!EOP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_mct_aVect_info - print out aVect info for debugging
!
! !DESCRIPTION:
!     print out attribute vector information
!
! !REVISION HISTORY:
!     2003 Jul 01 - B. Kauffman, T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_mct_aVect_info(flag,aVect,comm,pe,fld,istr)

! !USES:
  
! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in)           :: flag  ! info level flag
   type(cpl_mct_aVect),intent(in)   :: aVect ! Attribute vector
   integer(IN),intent(in),optional  :: comm  ! MPI communicator
   integer(IN),intent(in),optional  :: pe    ! processor number
   character(*),intent(in),optional :: fld   ! fld
   character(*),intent(in),optional :: istr  ! string for print

!EOP

   !--- local ---
   integer(IN)          :: i,j,k,n      ! generic indicies
   integer(IN)          :: ks,ke        ! start and stop k indices
   integer(IN)          :: nflds        ! number of flds in AV to diagnose
   integer(IN)          :: nsize        ! grid point size of AV
   type(cpl_mct_string) :: item         ! mct string
   character(CL)        :: itemc        ! item converted to char
   integer(IN)          :: comm_loc     ! local variable for comm
   integer(IN)          :: pe_loc       ! local variable for pe
   logical              :: commOK       ! is comm available
   logical              :: peOK         ! is pe available
   real(R8),allocatable :: minl(:)      ! local  min
   real(R8),allocatable :: ming(:)      ! global min
   real(R8),allocatable :: maxl(:)      ! local  max
   real(R8),allocatable :: maxg(:)      ! global max

   !--- formats ---
   character(*),parameter :: subName = '(cpl_mct_aVect_info) '
   character(*),parameter :: F00 = "('(cpl_mct_aVect_info) ',8a)"
   character(*),parameter :: F01 = "('(cpl_mct_aVect_info) ',a,i9)"
   character(*),parameter :: F02 = "('(cpl_mct_aVect_info) ',240a)"
   character(*),parameter :: F03 = "('(cpl_mct_aVect_info) ',a,2es11.3,i4,2x,a)"

!-------------------------------------------------------------------------------
! NOTE: has hard-coded knowledge/assumptions about mct aVect data type internals
!-------------------------------------------------------------------------------

   commOK = .false.
   peOK   = .false.

   if (present(pe)) then
     peOK = .true.
     pe_loc = pe
   endif
   if (present(comm)) then
     commOK = .true.
     comm_loc = comm
     if (.not.PEOK) then
       call shr_mpi_commrank(comm,pe_loc,subName)
       peOK = .true.
     endif
   endif

   nsize = cpl_mct_aVect_lsize(aVect)

   if (present(fld)) then
     nflds = 1
     ks = cpl_mct_aVect_indexRA(aVect,fld,perrWith=subName)
     ke = ks
   else
     nflds = cpl_mct_aVect_nRAttr(aVect)
     ks = 1
     ke = nflds
   endif

   if (flag >= 1) then
     if (present(istr)) write(6,*) trim(istr)
     write(6,F01) "local size =",nsize
     if (associated(aVect%iList%bf)) write(6,F02) "iList = ",aVect%iList%bf
     if (associated(aVect%rList%bf)) write(6,F02) "rList = ",aVect%rList%bf
   endif

   if (flag >= 2) then

     allocate(minl(nflds))
     allocate(maxl(nflds))

     do k=ks,ke
       minl(k) = minval(aVect%rAttr(k,:))
       maxl(k) = maxval(aVect%rAttr(k,:))
     enddo

     if (flag >= 4 .and. commOK) then
       allocate(ming(nflds))
       allocate(maxg(nflds))
       ming = 0.
       maxg = 0.
       call shr_mpi_min(minl,ming,comm,subName)
       call shr_mpi_max(maxl,maxg,comm,subName)
     endif

     do k=ks,ke
       call cpl_mct_aVect_getRList(item,k,aVect)
       itemc = cpl_mct_string_toChar(item)
       call cpl_mct_string_clean(item)
       write(6,F03) 'l min/max ',minl(k),maxl(k),k,trim(itemc)
       if (flag >= 3 .and. commOK) then
         if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
           write(6,F03) 'g min/max ',ming(k),maxg(k),k,trim(itemc)
         endif
       endif
       if (flag >= 4 .and. commOK) then
         if ((peOK .and. pe_loc == 0) .or. .not.peOK) then
           write(6,*) trim(subName),'g min/max ',ming(k),maxg(k),k,trim(itemc)
         endif
       endif
     enddo

      deallocate(minl)
      deallocate(maxl)
      if (flag >= 4 .and. commOK) then
         deallocate(ming)
         deallocate(maxg)
      endif

   endif

   call shr_sys_flush(6)

end subroutine cpl_mct_aVect_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_mct_aVect_getRAttr - get real F90 array data out of an aVect
!
! !DESCRIPTION:
!     Get real F90 array data (a single field) out of an MCT aVect.
!
! !REVISION HISTORY:
!     2002 Apr xx - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_mct_aVect_getRAttr(aVect,str,data,rcode)

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_mct_aVect),intent(in)  :: aVect    ! an Attribute vector
   character(*)       ,intent(in)  :: str      ! field name string
   real(R8)           ,intent(out) :: data(:)  ! an F90 array
   integer(IN)        ,intent(out) :: rcode    ! return code

!EOP

   !--- local ---
   integer(IN) :: k,n,m
   integer(IN) :: aVsize

   !--- formats ---
   character(*),parameter :: subName = "(cpl_mct_aVect_getRAttr) "
   character(*),parameter :: F00 = "('(cpl_mct_aVect_getRAttr) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rcode = 0

   n = cpl_mct_aVect_lsize(aVect)
   m = size(data)
   if (n /= m) then
      write(6,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
      data = cpl_const_spval
      rcode = 1
      return
   end if
   
   k = cpl_mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
   if ( k < 1) then
      write(6,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
      data = cpl_const_spval
      rcode = 2
      return
   end if

   data(:) = aVect%rAttr(k,:)

end subroutine cpl_mct_aVect_getRAttr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_mct_aVect_putRAttr - put real F90 array data into an aVect
!
! !DESCRIPTION:
!     Put real F90 array data (a single field) into an MCT aVect.
!
! !REVISION HISTORY:
!     2002 Apr xx - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_mct_aVect_putRAttr(aVect,str,data,rcode)

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_mct_aVect),intent(out) :: aVect ! Attribute vector
   character(*)       ,intent(in)  :: str
   real(R8)           ,intent(in)  :: data(:)
   integer(IN)        ,intent(out) :: rcode

!EOP

   !--- local ---
   integer(IN) :: k,n,m
   integer(IN) :: aVsize

   !--- formats ---
   character(*),parameter :: subName = "(cpl_mct_aVect_putRAttr) "
   character(*),parameter :: F00 = "('(cpl_mct_aVect_putRAttr) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rcode = 0

   n = cpl_mct_aVect_lsize(aVect)
   m = size(data)
   if (n /= m) then
      write(6,*) subName,"ERROR: size aV,data,attr = ",n,m,trim(str)
      rcode = 1
      return
   end if
   
   k = cpl_mct_aVect_indexRA(aVect,trim(str) ,perrWith=subName)
   if ( k < 1) then
      write(6,*) subName,"ERROR: attribute not found, var = ",trim(str),", k=",k
      rcode = 2
      return
   end if

   aVect%rAttr(k,:) = data(:) 

end subroutine cpl_mct_aVect_putRAttr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_mct_aVect_accum - accumulate attributes from one aVect to another
!
! !DESCRIPTION:
! This routine accumulates from input argment {\tt aVin} into the output 
! {\tt AttrVect} argument {\tt aVout} the real and integer attributes specified in 
! input {\tt CHARACTER} argument {\tt iList} and {\tt rList}. The attributes can
! be listed in any order.  If neither {\tt iList} nor {\tt rList} are provided, 
! all attributes shared between {\tt aVin} and {\tt aVout} will be copied.
!
! If any attributes in {\tt aVout} have different names but represent the
! the same quantity and should still be copied, you must provide a translation
! argument {\tt TrList} and/or {\tt TiList}.  The translation arguments should
! be identical to the {\tt rList} or {\tt iList} but with the correct {\tt aVout}
! name subsititued at the appropriate place.
!
! {\bf N.B.:}  This routine will fail if the {\tt aVout} is not initialized or
! if any of the specified attributes are not present in either {\tt aVout} or {\tt aVin}.
!
! !REVISION HISTORY:
!    2002 Sep 15 - ? - initial version.
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_mct_aVect_accum(aVin, rList, TrList, iList, TiList, aVout)


! !USES:

   use m_die ,          only : die
   use m_stdio ,        only : stderr
   use m_String ,       only : String_toChar => toChar
   use m_String ,       only : String
   use m_String ,       only : String_init
   use m_String ,       only : String_clean => clean
   use m_List ,         only : List
   use m_List,          only : List_get => get
   use m_List,          only : List_nullify => nullify
   use m_List,          only : List_clean => clean
   use m_List,          only : init,nitem
   use m_AttrVect,      only : AttrVect
   use m_AttrVect,      only : lsize
   use m_AttrVect,      only : SharedAttrIndexList

   implicit none

!INPUT/OUTPUT PARAMETERS

   type(AttrVect)        ,intent(in)    :: aVin
   character(*), optional,intent(in)    :: iList
   character(*), optional,intent(in)    :: rList
   character(*), optional,intent(in)    :: TiList
   character(*), optional,intent(in)    :: TrList
   type(AttrVect)        ,intent(inout) :: aVout

!EOP

   !--- local ---
   type(List)   :: rcpList       !  The list of real attributes to accum
   type(List)   :: icpList       !  The list of integer attributes to accum
   type(List)   :: TrcpList      !  Names of output attributes corresponding to input
   type(List)   :: TicpList      !  Names of output attributes corresponding to input
   type(String) :: attr          !  an individual attribute
   type(String) :: attr2         !  an individual attribute
   integer(IN)  :: i,j           ! generic indicies
   integer(IN)  :: rcode         ! return code
   integer(IN)  :: inx,outx
   integer(IN)  :: num_indices   ! Overlapping attribute index number

   !--- Overlapping attribute index storage arrays: ---
   integer(IN), dimension(:), pointer :: aVinindices, aVoutindices

   character(7) :: data_flag      ! character variable used as data type flag

   !--- formats ---
   character(*),parameter :: myname_='cpl_mct_accum'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call List_nullify(rcpList)
   call List_nullify(icpList)
   call List_nullify(TrcpList)
   call List_nullify(TicpList)

   if (lsize(aVin) .ne. lsize(aVout)) then
      write(stderr,'(2a)')myname_, &
      'MCTERROR:  Input aV and output aV do not have the same size'
      write(stderr,*)myname_, &
      'MCTERROR: ',lsize(aVin),lsize(aVout)
      call die(myname_,'lsize check',rcode)
   endif

   !----------------------------------------------------------------------------
   ! Accum the listed real attributes
   !----------------------------------------------------------------------------
   if ( present(rList) .and. (len_trim(rList)>0) ) then

      call init(rcpList,rList)    ! init.List()

      !--- check translation list ---
      if ( present(TrList) .and. (len_trim(TrList)>0) ) then
         call init(TrcpList,TrList)
         if ( nitem(rcpList) .ne. nitem(TrcpList)) then
            write(stderr,'(2a)')myname_, &
            'MCTERROR:  Input rList and TrList do not have the same size'
            call die(myname_,'nitem TrList check',rcode)
         endif
      endif

      if (nitem(rcpList) .ge. 1) then
         do i=1,lsize(aVin)
         do j=1,nitem(rcpList)
            call List_get(attr,j,rcpList)
            if (present(TrList)) then
               call List_get(attr2,j,TrcpList)
            else
               call String_init(attr2,attr)
            endif
            inx=cpl_mct_aVect_indexRA(aVin,String_toChar(attr),dieWith=myname_//'real aVin')
            outx=cpl_mct_aVect_indexRA(aVout,String_toChar(attr2),dieWith=myname_//'real aVout')
            aVout%rAttr(outx,i)=aVout%rAttr(outx,i)+aVin%rAttr(inx,i)
            call String_clean(attr)
            call String_clean(attr2)
         enddo
         enddo
      endif

     call List_clean(rcpList)
     if (present(TrList)) call List_clean(TrcpList)

   endif

   !----------------------------------------------------------------------------
   ! Accum the listed integer attributes 
   !----------------------------------------------------------------------------
   if ( present(iList) .and. (len_trim(iList)>0) ) then

      call init(icpList,iList)    ! init.List()

      !--- check translation list ---
      if ( present(TiList) .and. (len_trim(TiList)>0) ) then
         call init(TicpList,TiList)
         if ( nitem(icpList) .ne. nitem(TicpList)) then
             write(stderr,'(2a)')myname_, &
            'MCTERROR:  Input iList and TiList do not have the same size'
            call die(myname_,'nitem TiList check',rcode)
         endif
      endif

     if (nitem(icpList) .ge. 1) then
        do i=1,lsize(aVin)
        do j=1,nitem(icpList)
           call List_get(attr,j,icpList)
           if (present(TiList)) then
              call List_get(attr2,j,TicpList)
           else
              call String_init(attr2,attr)
           endif
           inx =cpl_mct_aVect_indexIA(aVin ,String_toChar(attr) ,dieWith=myname_//'int aVin')
           outx=cpl_mct_aVect_indexIA(aVout,String_toChar(attr2),dieWith=myname_//'int aVout')
           aVout%iAttr(outx,i)=aVout%iAttr(outx,i)+aVin%iAttr(inx,i)
           call String_clean(attr)
           call String_clean(attr2)
        enddo
        enddo
     endif

     call List_clean(icpList)
     if (present(TrList)) call List_clean(TicpList)

   endif

   !----------------------------------------------------------------------------
   ! if neither rList nor iList is present, accum shared attibutes from in to out
   !----------------------------------------------------------------------------
   if ( .not.present(rList) .and. .not.present(iList)) then

      data_flag = 'REAL'
      call SharedAttrIndexList(aVin, aVout, data_flag, num_indices, &
                                aVinindices, aVoutindices)
      if (num_indices .gt. 0) then
#ifdef CPP_VECTOR
         do j=1,num_indices
!CDIR SELECT(VECTOR)
         do i=1,lsize(aVin)
#else
         do i=1,lsize(aVin)
         do j=1,num_indices
#endif
            aVout%rAttr(aVoutindices(j),i)= &
            & aVout%rAttr(aVoutindices(j),i)+aVin%rAttr(aVinindices(j),i)
         enddo
         enddo
      endif
      deallocate(aVinindices, aVoutindices,stat=rcode)
      if (rcode /= 0) call die(myname_,'deallocate real(Vinindices...',rcode)

      data_flag = 'INTEGER'
      call SharedAttrIndexList(aVin, aVout, data_flag, num_indices, &
                                aVinindices, aVoutindices)
      if (num_indices .gt. 0) then
#ifdef CPP_VECTOR
         do j=1,num_indices
!CDIR SELECT(VECTOR)
         do i=1,lsize(aVin)
#else
         do i=1,lsize(aVin)
         do j=1,num_indices
#endif
            aVout%iAttr(aVoutindices(j),i)= &
            & aVout%iAttr(aVoutindices(j),i)+aVin%iAttr(aVinindices(j),i)
         enddo
         enddo
      endif
      deallocate(aVinindices, aVoutindices,stat=rcode)
      if (rcode /= 0) call die(myname_,'deallocate int(Vinindices...',rcode)

   endif

end subroutine cpl_mct_aVect_accum

!===============================================================================

end module cpl_mct_mod
