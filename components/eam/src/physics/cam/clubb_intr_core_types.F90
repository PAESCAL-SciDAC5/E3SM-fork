module clubb_intr_core_types

  use clubb_precision,  only: core_rknd
  use cam_abortutils,   only: endrun

  implicit none

  public 

  type core_state_t
    real(core_rknd), dimension(:), allocatable :: p_in_Pa
    real(core_rknd), dimension(:), allocatable :: exner
    real(core_rknd), dimension(:), allocatable :: rho_zt
    real(core_rknd), dimension(:), allocatable :: rho_zm
    real(core_rknd), dimension(:), allocatable :: rho_ds_zt
    real(core_rknd), dimension(:), allocatable :: rho_ds_zm
    real(core_rknd), dimension(:), allocatable :: invrs_rho_ds_zm
    real(core_rknd), dimension(:), allocatable :: invrs_rho_ds_zt
  end type core_state_t

  type core_forcing_t
    real(core_rknd), dimension(:), allocatable :: fld
  end type core_forcing_t

  type core_srcflx_t
    real(core_rknd), dimension(:), allocatable :: fld
  end type core_srcflx_t

contains

  !------------------------------------------------------ 
  subroutine clubb_core_state_alloc( core_state, pverp )

    type(core_state_t),intent(inout) :: core_state
    integer,           intent(in)    :: pverp

    character(len=100) :: routine='allocate_clubb_core_state'
    integer :: ierr

    allocate( core_state% p_in_Pa         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% exner           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rho_zt          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rho_zm          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rho_ds_zt       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rho_ds_zm       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% invrs_rho_ds_zm (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% invrs_rho_ds_zt (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

   !(pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
   !(pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
   !(pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

  end subroutine clubb_core_state_alloc

  !----------------------------------------------- 
  subroutine clubb_core_state_dealloc( core_state )

    type(core_state_t),intent(inout) :: core_state

    character(len=100) :: routine='clubb_core_state_dealloc'
    integer :: ierr

    deallocate( core_state% p_in_Pa        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% exner          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rho_zt         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rho_zm         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rho_ds_zt      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rho_ds_zm      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% invrs_rho_ds_zm, stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% invrs_rho_ds_zt, stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

  end subroutine clubb_core_state_dealloc

end module clubb_intr_core_types
