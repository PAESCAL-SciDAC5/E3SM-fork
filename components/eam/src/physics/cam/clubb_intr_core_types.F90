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

    real(core_rknd), dimension(:), allocatable :: thv_ds_zt
    real(core_rknd), dimension(:), allocatable :: thv_ds_zm

    real(core_rknd), dimension(:), allocatable :: wm_zt
    real(core_rknd), dimension(:), allocatable :: wm_zm

    real(core_rknd), dimension(:), allocatable :: um
    real(core_rknd), dimension(:), allocatable :: vm

    real(core_rknd), dimension(:), allocatable :: thlm
    real(core_rknd), dimension(:), allocatable :: rtm
    real(core_rknd), dimension(:), allocatable :: rcm

    real(core_rknd), dimension(:), allocatable :: up2
    real(core_rknd), dimension(:), allocatable :: vp2
    real(core_rknd), dimension(:), allocatable :: wp2
    real(core_rknd), dimension(:), allocatable :: wp3

    real(core_rknd), dimension(:), allocatable :: upwp
    real(core_rknd), dimension(:), allocatable :: vpwp

    real(core_rknd), dimension(:), allocatable :: rtp2
    real(core_rknd), dimension(:), allocatable :: thlp2
    real(core_rknd), dimension(:), allocatable :: wprtp
    real(core_rknd), dimension(:), allocatable :: wpthlp
    real(core_rknd), dimension(:), allocatable :: rtpthlp

    real(core_rknd), dimension(:), allocatable :: wp2thvp
    real(core_rknd), dimension(:), allocatable :: wpthvp
    real(core_rknd), dimension(:), allocatable :: rtpthvp
    real(core_rknd), dimension(:), allocatable :: thlpthvp

    real(core_rknd), dimension(:), allocatable :: wprcp
    real(core_rknd), dimension(:), allocatable :: thlprcp

    real(core_rknd), dimension(:), allocatable :: khzm
    real(core_rknd), dimension(:), allocatable :: khzt
    real(core_rknd), dimension(:), allocatable :: cloud_frac
    real(core_rknd), dimension(:), allocatable :: qclvar      ! cloud water variance [kg^2/kg^2]

    ! Passive scalars
    real(core_rknd), dimension(:,:), allocatable :: edsclr    ! eddy passive scalars

  end type core_state_t

  !------------------------------------
  type core_forcing_t

    real(core_rknd), dimension(:), allocatable :: um    ! u wind forcing (thermodynamic levels)         [m/s/s]
    real(core_rknd), dimension(:), allocatable :: vm    ! v wind forcing (thermodynamic levels)         [m/s/s]

    real(core_rknd), dimension(:), allocatable :: thlm  ! theta_l forcing (thermodynamic levels)        [K/s]
    real(core_rknd), dimension(:), allocatable :: rtm   ! r_t forcing (thermodynamic levels)            [(kg/kg)/s]

    real(core_rknd), dimension(:), allocatable :: rtp2
    real(core_rknd), dimension(:), allocatable :: thlp2
    real(core_rknd), dimension(:), allocatable :: rtpthlp 
    real(core_rknd), dimension(:), allocatable :: wprtp
    real(core_rknd), dimension(:), allocatable :: wpthlp

    ! Passive scalars
    real(core_rknd), dimension(:,:), allocatable :: edsclr    ! eddy passive scalars

  end type core_forcing_t

  !------------------------------------
  type core_sfc_t

    real(core_rknd) :: upwp
    real(core_rknd) :: vpwp

    real(core_rknd) :: wpthlp
    real(core_rknd) :: wprtp
 
    real(core_rknd), dimension(:), allocatable :: wpedsclrp

  end type core_sfc_t

contains

  !------------------------------------------------------ 
  subroutine clubb_core_fld_alloc( core_state, core_forcing, core_sfc, pverp, ntracer )

    type(core_state_t),  intent(inout) :: core_state
    type(core_forcing_t),intent(inout) :: core_forcing
    type(core_sfc_t),    intent(inout) :: core_sfc

    integer,           intent(in)    :: pverp
    integer,           intent(in)    :: ntracer

    character(len=100) :: routine='clubb_core_fld_alloc'
    integer :: ierr

    allocate( core_state% p_in_Pa         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% exner           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% rho_zt          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rho_zm          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rho_ds_zt       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rho_ds_zm       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% invrs_rho_ds_zm (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% invrs_rho_ds_zt (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% thv_ds_zt       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% thv_ds_zm       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% wm_zt           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% wm_zm           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% um              (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% vm              (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% thlm            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rtm             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rcm             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% up2             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% vp2             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% wp2             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% wp3             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% upwp            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% vpwp            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% rtp2            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% thlp2           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% wprtp           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% wpthlp          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rtpthlp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% wp2thvp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% wpthvp          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% rtpthvp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% thlpthvp        (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% wprcp           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% thlprcp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% khzm            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% khzt            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% cloud_frac      (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_state% qclvar          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_state% edsclr  (pverp,ntracer), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !-------------------------
    ! Forcing terms
    !-------------------------
    allocate( core_forcing% um            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_forcing% vm            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_forcing% thlm          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_forcing% rtm           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_forcing% rtp2          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_forcing% thlp2         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_forcing% wprtp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_forcing% wpthlp        (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_forcing% rtpthlp       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_forcing% edsclr(pverp,ntracer), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !-------------------------
    ! Sfc values
    !-------------------------
    ! Only need to allocate array(s) for passive scalar(s)
    allocate( core_sfc% wpedsclrp(ntracer), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

  end subroutine clubb_core_fld_alloc

  !----------------------------------------------- 
  subroutine clubb_core_fld_dealloc( core_state, core_forcing, core_sfc )

    type(core_state_t),  intent(inout) :: core_state
    type(core_forcing_t),intent(inout) :: core_forcing
    type(core_sfc_t),    intent(inout) :: core_sfc

    character(len=100) :: routine='clubb_core_fld_dealloc'
    integer :: ierr

    deallocate( core_state% p_in_Pa        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% exner          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% rho_zt         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rho_zm         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rho_ds_zt      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rho_ds_zm      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% invrs_rho_ds_zm, stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% invrs_rho_ds_zt, stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% thv_ds_zt      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% thv_ds_zm      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% wm_zt          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% wm_zm          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% um             , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% vm             , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% thlm           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rtm            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rcm            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% up2            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% vp2            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% wp2            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% wp3            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% upwp           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% vpwp           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% rtp2           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% thlp2          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% wprtp          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% wpthlp         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rtpthlp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% wp2thvp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% wpthvp         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% rtpthvp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% thlpthvp       , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% wprcp          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% thlprcp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% khzm           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% khzt           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% cloud_frac     , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_state% qclvar         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_state% edsclr         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !-------------------------
    ! Forcing terms
    !-------------------------
    deallocate( core_forcing% um           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_forcing% vm           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_forcing% thlm         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_forcing% rtm          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_forcing% rtp2         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_forcing% thlp2        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_forcing% wprtp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_forcing% wpthlp       , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_forcing% rtpthlp      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_forcing% edsclr       , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !-------------------------
    ! Sfc values
    !-------------------------
    ! Only need to deallocate array(s) for passive scalar(s)

    deallocate( core_sfc% wpedsclrp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

  end subroutine clubb_core_fld_dealloc

end module clubb_intr_core_types
