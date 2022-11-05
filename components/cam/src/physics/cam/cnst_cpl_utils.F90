module cnst_cpl_utils
!---------------------------------------------------------------------------------
! Utitily subroutines used by various process coupling schemes 
!---------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none

contains

!------------------------------------------------------------------------------------------
subroutine copy_dqdt_from_pbuf_to_ptend( pbuf, ptend, pbuf_fldname, ims, ime, pcols, pver )

  use physics_types,  only: physics_ptend
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  ! Arguments

  type(physics_buffer_desc), pointer :: pbuf(:)  ! physics puffer
  type(physics_ptend) :: ptend                   ! indivdual parameterization tendencies
  character(len=*),intent(in) :: pbuf_fldname    ! name of the pbuf field to copy values to 
  integer, intent(in) :: ims, ime                ! start and end indices of tracers for which 
                                                 ! the tendencies should be copied to pbuf
  integer, intent(in) :: pcols, pver             ! # of grid columns and vertical layers
                                                 ! for which the tendencies should be copied 

  ! Local variables

  real(r8),pointer :: ptr2d(:,:)  ! (pcols,pver) pointer pointing to a slide of data in pbuf
  integer :: fldidx ! pbuf field index
  integer :: im     ! tracer loop index

  !---------------------------------------
  fldidx = pbuf_get_index( pbuf_fldname )

  do im = ims,ime
     call pbuf_get_field( pbuf, fldidx, ptr2d,   &! in, in, out
                          start=(/1,1,im/),      &! in
                          kount=(/pcols,pver,1/) )! in

     ptend%q(:,:,im) = ptr2d(:,:)
  end do

end subroutine copy_dqdt_from_pbuf_to_ptend


!-------------------------------------------------------------------------------------------------
subroutine calculate_dqdt_and_save_to_pbuf( state_old, state_new, dtime, pbuf, pbuf_fldname, &
                                            ims,ime, pcols,pver )

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  ! Arguments

  type(physics_state),intent(in) :: state_old, state_new
  real(r8),           intent(in) :: dtime
  type(physics_buffer_desc), pointer :: pbuf(:)  ! physics puffer
  character(len=*),intent(in)    :: pbuf_fldname ! name of the pbuf field to copy values to 
  integer, intent(in) :: ims, ime                ! start and end indices of tracers for which 
                                                 ! the tendencies should be copied to pbuf
  integer, intent(in) :: pcols, pver             ! # of grid columns and vertical layers
                                                 ! for which the tendencies should be copied 
  ! Local variables

  real(r8),pointer :: ptr2d(:,:)  ! (pcols,pver) pointer pointing to a slide of data in pbuf
  integer :: fldidx ! pbuf field index
  integer :: im     ! tracer loop index

  fldidx = pbuf_get_index( pbuf_fldname )

  !---------------------------------------
  do im = ims,ime
     call pbuf_get_field( pbuf, fldidx, ptr2d,   &! in, in, out
                          start=(/1,1,im/),      &! in
                          kount=(/pcols,pver,1/) )! in

     ptr2d(:,:) = ( state_new%q(:,:,im) - state_old%q(:,:,im) )/dtime
  end do

end subroutine calculate_dqdt_and_save_to_pbuf

end module cnst_cpl_utils
