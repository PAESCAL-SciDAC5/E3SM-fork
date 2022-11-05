module cnst_cpl_utils
!------------------------------------------------------------------------------------------
! Utitily subroutines used by various process coupling schemes for gas and aerosol tracers
!------------------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none

contains

subroutine get_saved_dqdt( pbuf, fldname, state, pcols, pver, ptend )
!------------------------------------------------------------------------------------------
! This subroutine retrieves the tracer mass tendencies (pdel*dqdt) saved as fldname in pbuf,
! convert them to mixing ratio tendencies, and save the results to ptend.  
!
! History
!  - Hui Wan, 2022-11, initial version
!------------------------------------------------------------------------------------------

  use constituents,   only: pcnst, cnst_get_ind, cnst_get_type_byind
  use physics_types,  only: physics_ptend, physics_ptend_init, physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  ! Arguments

  type(physics_buffer_desc), pointer :: pbuf(:)  ! physics puffer
  character(len=*),intent(in) :: fldname         ! name of the pbuf field to retrieve info from 
  type(physics_state), intent(in) :: state       ! state variable containing rpdel(dry) and dimensin size info
  integer, intent(in) :: pcols, pver             ! # of grid columns and vertical layers in this chunk

  type(physics_ptend),intent(out) :: ptend       ! contains dqdt to be passed to the calling routine

  ! Local variables

  real(r8),pointer :: ptr2d(:,:)  ! pointer a slice of data in pbuf with the shape (pcols,pver)
  integer :: fldidx               ! pbuf field index
  integer :: im, ims              ! tracer loop index and its start value
  logical :: lq(pcnst)
  integer :: ncol

  !-----------
  ! Inquiries
  !-----------
  ncol = state%ncol
  call cnst_get_ind('NUMSNO', ims); ims = ims+1   ! Exclude water species from the list of tracers to be processed

  !------------------
  ! Initialize ptend
  !------------------
  lq(ims:pcnst) = .TRUE.
  call physics_ptend_init(ptend, state%psetcols, fldname, lq=lq)

  !-------------------------------------------------
  ! Retrieve pdel*dqdt from pbuf; convert to dqdt
  !-------------------------------------------------
  fldidx = pbuf_get_index( fldname )

  do im = ims,pcnst

     call pbuf_get_field( pbuf, fldidx, ptr2d,   &! in, in, out
                          start=(/1,1,im/),      &! in
                          kount=(/pcols,pver,1/) )! in

     if (cnst_get_type_byind(im).eq.'dry') then
        ptend%q(:ncol,:,im) = ptr2d(:ncol,:)*state%rpdeldry(:ncol,:)
     else
        ptend%q(:ncol,:,im) = ptr2d(:ncol,:)*state%rpdel(:ncol,:)
     end if

  end do

end subroutine get_saved_dqdt 


subroutine calculate_dqdt_and_save_to_pbuf( state_old, state_new, dtime, pbuf, fldname, pcols,pver )
!-------------------------------------------------------------------------------------------------
! This subroutine diagnoses tracer mixing ratio tendencies from two state snapshots and then
! convert them to mass tendencies and save to pbuf.
!
! History
!  - Hui Wan, 2022-11, initial version
!-------------------------------------------------------------------------------------------------

  use constituents,   only: pcnst, cnst_get_ind, cnst_get_type_byind
  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  ! Arguments

  type(physics_state),intent(in) :: state_old, state_new
  real(r8),           intent(in) :: dtime
  type(physics_buffer_desc), pointer :: pbuf(:)  ! physics puffer
  character(len=*),intent(in)        :: fldname  ! name of the pbuf field to copy values to 
  integer, intent(in) :: pcols, pver             ! # of grid columns and vertical layers in this chunk

  ! Local variables

  real(r8),pointer :: ptr2d(:,:)  ! pointer a slice of data in pbuf with the shape (pcols,pver)
  integer :: fldidx               ! pbuf field index
  integer :: im, ims              ! tracer loop index and its start value
  integer :: ncol

  !-----------
  ! Inquiries
  !-----------
  ncol = state_old%ncol
  call cnst_get_ind('NUMSNO', ims); ims = ims+1   ! Exclude water species from the list of tracers to be processed

  !-----------------------------------------------------
  ! Calculate tracer mass tendencies and save to pbuf
  !-----------------------------------------------------
  fldidx = pbuf_get_index( fldname )

  do im = ims,pcnst
     call pbuf_get_field( pbuf, fldidx, ptr2d,   &! in, in, out
                          start=(/1,1,im/),      &! in
                          kount=(/pcols,pver,1/) )! in

     if (cnst_get_type_byind(im).eq.'dry') then
        ptr2d(:ncol,:) = ( state_new%q(:ncol,:,im)*state_new%pdeldry(:ncol,:) &
                          -state_old%q(:ncol,:,im)*state_old%pdeldry(:ncol,:) )/dtime
     else
        ptr2d(:ncol,:) = ( state_new%q(:ncol,:,im)*state_new%pdel(:ncol,:) &
                          -state_old%q(:ncol,:,im)*state_old%pdel(:ncol,:) )/dtime
     end if

  end do

end subroutine calculate_dqdt_and_save_to_pbuf

end module cnst_cpl_utils
