module mo_precision_utils

use shr_kind_mod,   only: r8 => shr_kind_r8
use shr_kind_mod,   only: r4 => shr_kind_r4
use rp_emulator
use cam_abortutils, only: endrun

interface change_precision
  module procedure change_precision_for_1d_array
  module procedure change_precision_for_2d_array
  module procedure change_precision_for_3d_array
end interface

contains

!==============================================================================
subroutine change_precision_for_1d_array( array_in, precision_out, array_out )

  real(r8),intent(in)              :: array_in (:)
  real(r8),intent(out),allocatable :: array_out(:)
  integer, intent(in)              :: precision_out

  ! local variables

  type(rpe_var),allocatable :: array_rpe(:)

  integer :: shape1d(1)
  integer :: ierr

  character(len=*), parameter :: subname= "change_precision_for_1d_array"

  !--------------------------------------
  ! Allocate memory for the output array
  !--------------------------------------
  shape1d = shape(array_in)
  allocate( array_out(shape1d(1)), stat=ierr )
  if (ierr/=0) call endrun(subname//": allocation error for array_out")

  !-----------------------------------------------------------------------------------
  ! Convert the input values to a different precision and assign to the output array
  !-----------------------------------------------------------------------------------
  select case(precision_out)
  case(16) ! Single precision.

    array_out = real(real(array_in,kind=r4),kind=r8)

  case(32) ! Double precision
    array_out = array_in

  case default ! Emulated precision

    allocate( array_rpe(shape1d(1)), stat=ierr )
    if (ierr/=0) call endrun(subname//": allocation error occurred for array_rpe")

    array_rpe%sbits = precision_out
    array_rpe = array_in
    array_out = array_rpe

    deallocate(array_rpe)

  end select
  !----------

end subroutine change_precision_for_1d_array

!==============================================================================
subroutine change_precision_for_2d_array( array_in, precision_out, array_out )

  real(r8),intent(in)              :: array_in (:,:)
  real(r8),intent(out),allocatable :: array_out(:,:)
  integer, intent(in)              :: precision_out

  ! local variables

  type(rpe_var),allocatable :: array_rpe(:,:)

  integer :: shape2d(2)
  integer :: ierr

  character(len=*), parameter :: subname= "change_precision_for_2d_array"

  !--------------------------------------
  ! Allocate memory for the output array
  !--------------------------------------
  shape2d = shape(array_in)
  allocate( array_out(shape2d(1),shape2d(2)), stat=ierr )
  if (ierr/=0) call endrun(subname//": allocation error for array_out")

  !-----------------------------------------------------------------------------------
  ! Convert the input values to a different precision and assign to the output array
  !-----------------------------------------------------------------------------------
  select case(precision_out)
  case(16) ! Single precision.

    array_out = real(real(array_in,kind=r4),kind=r8)

  case(32) ! Double precision

    array_out = array_in

  case default ! Emulated precision

    allocate( array_rpe(shape2d(1),shape2d(2)), stat=ierr )
    if (ierr/=0) call endrun(subname//": allocation error occurred for array_rpe")

    array_rpe%sbits = precision_out
    array_rpe = array_in
    array_out = array_rpe

    deallocate(array_rpe)

  end select
  !----------

end subroutine change_precision_for_2d_array

!==============================================================================
subroutine change_precision_for_3d_array( array_in, precision_out, array_out )

  real(r8),intent(in)              :: array_in (:,:,:)
  real(r8),intent(out),allocatable :: array_out(:,:,:)
  integer, intent(in)              :: precision_out

  ! local variables

  type(rpe_var),allocatable :: array_rpe(:,:,:)

  integer :: shape3d(3)
  integer :: ierr

  character(len=*), parameter :: subname= "change_precision_for_3d_array"

  !--------------------------------------
  ! Allocate memory for the output array
  !--------------------------------------
  shape3d = shape(array_in)
  allocate( array_out(shape3d(1),shape3d(2),shape3d(3)), stat=ierr )
  if (ierr/=0) call endrun(subname//": allocation error for array_out")

  !-----------------------------------------------------------------------------------
  ! Convert the input values to a different precision and assigne to the output array
  !-----------------------------------------------------------------------------------
  select case(precision_out)
  case(16) ! Single precision.

    array_out = real(real(array_in,kind=r4),kind=r8)

  case(32) ! Double precision

    array_out = array_in

  case default ! Emulated precision

    allocate( array_rpe(shape3d(1),shape3d(2),shape3d(3)), stat=ierr )
    if (ierr/=0) call endrun(subname//": allocation error occurred for array_rpe")

    array_rpe%sbits = precision_out
    array_rpe = array_in
    array_out = array_rpe

    deallocate(array_rpe)

  end select
  !----------

end subroutine change_precision_for_3d_array


end module mo_precision_utils
