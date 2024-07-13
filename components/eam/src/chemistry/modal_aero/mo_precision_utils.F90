module mo_precision_utils

use shr_kind_mod,   only: r8 => shr_kind_r8
use shr_kind_mod,   only: r4 => shr_kind_r4
use cam_abortutils, only: endrun

interface change_precision
  module procedure change_precision_for_1d_array
  module procedure change_precision_for_2d_array
  module procedure change_precision_for_3d_array
end interface

contains

!==============================================================================
!==============================================================================
subroutine change_precision_for_1d_array( array_in, precision_out, array_out )

  real(r8),intent(in)              :: array_in (:)
  real(r8),intent(out),allocatable :: array_out(:)
  integer, intent(in)              :: precision_out

  ! local variables

  integer :: shape1d(1)
  integer :: ierr

  character(len=*), parameter :: subname= "change_precision_for_1d_array"

  ! allocate memory

  shape1d = shape(array_in)

  allocate( array_out(shape1d(1)), stat=ierr )
  if (ierr/=0) call endrun(subname//": allocation error")

  ! convert the input values to a different precision and save to output array

  select case(precision_out)
  case(16)
    array_out = real(real(array_in,kind=r4),kind=r8)
  case(32)
    array_out = array_in
  case default
    call endrun(subname//": unknown choice of precision_out.")
  end select

end subroutine change_precision_for_1d_array

!==============================================================================
!==============================================================================
subroutine change_precision_for_2d_array( array_in, precision_out, array_out )

  real(r8),intent(in)              :: array_in (:,:)
  real(r8),intent(out),allocatable :: array_out(:,:)
  integer, intent(in)              :: precision_out

  ! local variables

  integer :: shape2d(2)
  integer :: ierr

  character(len=*), parameter :: subname= "change_precision_for_2d_array"

  ! allocate memory

  shape2d = shape(array_in)

  allocate( array_out(shape2d(1),shape2d(2)), stat=ierr )
  if (ierr/=0) call endrun(subname//": allocation error")

  ! convert the input values to a different precision and save to output array

  select case(precision_out)
  case(16)
    array_out = real(real(array_in,kind=r4),kind=r8)
  case(32)
    array_out = array_in
  case default
    call endrun(subname//": unknown choice of precision_out.")
  end select

end subroutine change_precision_for_2d_array

!==============================================================================
!==============================================================================
subroutine change_precision_for_3d_array( array_in, precision_out, array_out )

  real(r8),intent(in)              :: array_in (:,:,:)
  real(r8),intent(out),allocatable :: array_out(:,:,:)
  integer, intent(in)              :: precision_out

  ! local variables

  integer :: shape3d(3)
  integer :: ierr

  character(len=*), parameter :: subname= "change_precision_for_3d_array"

  ! allocate memory

  shape3d = shape(array_in)

  allocate( array_out(shape3d(1),shape3d(2),shape3d(3)), stat=ierr )
  if (ierr/=0) call endrun(subname//": allocation error")

  ! convert the input values to a different precision and save to output array

  select case(precision_out)
  case(16)
    array_out = real(real(array_in,kind=r4),kind=r8)
  case(32)
    array_out = array_in
  case default
    call endrun(subname//": unknown choice of precision_out.")
  end select

end subroutine change_precision_for_3d_array


end module mo_precision_utils
