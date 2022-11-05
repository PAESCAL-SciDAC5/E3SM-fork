module cnst_cpl_utils
!---------------------------------------------------------------------------------
! Utitily subroutines used by various process coupling schemes 
!---------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only : endrun

  implicit none

  private :: copy_dqdt_from_ptend_to_pbuf

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
subroutine copy_dqdt_from_ptend_to_pbuf( ptend, pbuf, pbuf_fldname, ims, ime, pcols, pver )

  use physics_types,  only: physics_ptend
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  type(physics_ptend) :: ptend                   ! indivdual parameterization tendencies
  type(physics_buffer_desc), pointer :: pbuf(:)  ! physics puffer
  character(len=*) :: pbuf_fldname               ! name of the pbuf field to copy values to 
  integer, intent(in) :: ims, ime                ! start and end indices of tracers for which 
                                                 ! the tendencies should be copied to pbuf
  integer, intent(in) :: pcols, pver             ! # of grid columns and vertical layers
                                                 ! for which the tendencies should be copied 

  real(r8),pointer :: ptr2d(:,:)  ! (pcols,pver) pointer pointing to a slide of data in pbuf
  integer :: fldidx ! pbuf field index
  integer :: im     ! tracer loop index

  fldidx = pbuf_get_index( pbuf_fldname )

  do im = ims,ime

     call pbuf_get_field( pbuf, fldidx, ptr2d, start=(/1,1,im/), kount=(/pcols,pver,1/) )

     if (ptend%lq(im)) then
        ptr2d(:,:) = ptend%q(:,:,im) 
     else
        ptr2d(:,:) = 0._r8 
     end if

  end do

end subroutine copy_dqdt_from_ptend_to_pbuf


subroutine copy_dqdt_from_pbuf_to_ptend( pbuf, ptend, pbuf_fldname, ims, ime, pcols, pver )

  use physics_types,  only: physics_ptend
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  type(physics_buffer_desc), pointer :: pbuf(:)  ! physics puffer
  type(physics_ptend) :: ptend                   ! indivdual parameterization tendencies
  character(len=*) :: pbuf_fldname               ! name of the pbuf field to copy values to 
  integer, intent(in) :: ims, ime                ! start and end indices of tracers for which 
                                                 ! the tendencies should be copied to pbuf
  integer, intent(in) :: pcols, pver             ! # of grid columns and vertical layers
                                                 ! for which the tendencies should be copied 

  real(r8),pointer :: ptr2d(:,:)  ! (pcols,pver) pointer pointing to a slide of data in pbuf
  integer :: fldidx ! pbuf field index
  integer :: im     ! tracer loop index


  fldidx = pbuf_get_index( pbuf_fldname )

  do im = ims,ime
     call pbuf_get_field( pbuf, fldidx, ptr2d, start=(/1,1,im/), kount=(/pcols,pver,1/) )
     ptend%q(:,:,im) = ptr2d(:,:)
  end do

end subroutine copy_dqdt_from_pbuf_to_ptend


!-------------------------------------------------------------------------------------------------
subroutine calculate_dqdt_and_save_to_pbuf( state_old, state_new, dtime, pbuf, pbuf_fldname, &
                                            ims,ime, pcols,pver )

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

  type(physics_state) :: state_old, state_new
  real(r8) :: dtime
  type(physics_buffer_desc), pointer :: pbuf(:)  ! physics puffer
  character(len=*) :: pbuf_fldname               ! name of the pbuf field to copy values to 
  integer, intent(in) :: ims, ime                ! start and end indices of tracers for which 
                                                 ! the tendencies should be copied to pbuf
  integer, intent(in) :: pcols, pver             ! # of grid columns and vertical layers
                                                 ! for which the tendencies should be copied 

  real(r8),pointer :: ptr2d(:,:)  ! (pcols,pver) pointer pointing to a slide of data in pbuf
  integer :: fldidx ! pbuf field index
  integer :: im     ! tracer loop index

  fldidx = pbuf_get_index( pbuf_fldname )

  do im = ims,ime
     call pbuf_get_field( pbuf, fldidx, ptr2d, start=(/1,1,im/), kount=(/pcols,pver,1/) )
     ptr2d(:,:) = ( state_new%q(:,:,im) - state_old%q(:,:,im) )/dtime
  end do

end subroutine calculate_dqdt_and_save_to_pbuf

end module cnst_cpl_utils
