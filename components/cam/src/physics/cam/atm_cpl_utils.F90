module atm_cpl_utils
!---------------------------------------------------------------------------------
! Utitily subroutines used by various process coupling schemes 
!---------------------------------------------------------------------------------

  use cam_abortutils, only : endrun

subroutine copy_dqdt_from_ptend_to_pbuf( ptend, pbuf, pbuf_fldname, ims, ime, pcols, pver )

  use physics_types,  only: physics_ptend
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field

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

  if (.not.ptend%lq) call endrun("copy_dqdt_from_ptend_to_pbuf: ptend%lq is .false.")

  fldidx = pbuf_get_index( pbuf_fldname )

  do im = ims,ime
     call pbuf_get_field( pbuf, fldidx, ptr2d, start=(/1,1,im/), kount=(/pcols,pver,1/) )
     ptr2d(:,:) = ptend%q(:,:,im) 
  end do

end subroutine copy_dqdt_from_ptend_to_pbuf


subroutine copy_dqdt_from_pbuf_to_ptend( pbuf, ptend, pbuf_fldname, ims, ime, pcols, pver )

  use physics_types,  only: physics_ptend
  use physics_buffer, only: physics_buffer_desc, pbuf_get_field

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

  if (.not.ptend%lq) call endrun("copy_dqdt_from_pbuf_to_ptend: ptend%lq is .false.")

  fldidx = pbuf_get_index( pbuf_fldname )

  do im = ims,ime
     call pbuf_get_field( pbuf, fldidx, ptr2d, start=(/1,1,im/), kount=(/pcols,pver,1/) )
     ptend%q(:,:,im) = ptr2d(:,:)
  end do

end subroutine copy_dqdt_from_pbuf_to_ptend

end module atm_cpl_utils
