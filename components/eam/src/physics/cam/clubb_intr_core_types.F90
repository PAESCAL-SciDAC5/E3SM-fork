module clubb_intr_core_types

  use clubb_precision,  only: core_rknd
  use cam_abortutils,   only: endrun

  implicit none

  public 

  !--------------------------------------------------------------------
  ! Auxiliary variables describing the grid-box mean atmospheric state
  !--------------------------------------------------------------------
  type core_auxil_t

    integer         :: icol, lchnk   ! column and chunk indicies in the host model
    real(core_rknd) :: lat, lon      ! latitude and longitude [radian]

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

  end type core_auxil_t

  !------------------------------
  ! CLUBB's prognostic variables
  !------------------------------
  type core_prog_t

    real(core_rknd), dimension(:), allocatable :: um
    real(core_rknd), dimension(:), allocatable :: vm

    real(core_rknd), dimension(:), allocatable :: thlm
    real(core_rknd), dimension(:), allocatable :: rtm

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

    real(core_rknd), dimension(:,:), allocatable :: edsclr    ! eddy passive scalars

  end type core_prog_t

  !------------------------------
  ! CLUBB's diagnostic variables
  !------------------------------
  type core_diag_t

    real(core_rknd), dimension(:), allocatable :: wp2thvp
    real(core_rknd), dimension(:), allocatable :: wpthvp
    real(core_rknd), dimension(:), allocatable :: rtpthvp
    real(core_rknd), dimension(:), allocatable :: thlpthvp

    real(core_rknd), dimension(:), allocatable :: wprcp
    real(core_rknd), dimension(:), allocatable :: thlprcp

    real(core_rknd), dimension(:), allocatable :: khzm
    real(core_rknd), dimension(:), allocatable :: khzt

    real(core_rknd), dimension(:), allocatable :: rcm
    real(core_rknd), dimension(:), allocatable :: cloud_frac
    real(core_rknd), dimension(:), allocatable :: qclvar        ! cloud water variance [kg^2/kg^2]
    real(core_rknd), dimension(:), allocatable :: rcm_in_layer  ! CLUBB output of in-cloud liq. wat. mix. ratio [kg/kg] 
    real(core_rknd), dimension(:), allocatable :: cloud_cover   ! CLUBB output of in-cloud cloud fraction       [fraction] 

  end type core_diag_t

  !-------------------------
  ! Forcing terms
  !-------------------------
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

    real(core_rknd), dimension(:,:), allocatable :: edsclr    ! forcing of eddy passive scalars

  end type core_forcing_t

  !-------------------------
  ! Sfc values
  !-------------------------
  type core_sfc_t

    real(core_rknd) :: upwp
    real(core_rknd) :: vpwp

    real(core_rknd) :: wpthlp
    real(core_rknd) :: wprtp
 
    real(core_rknd), dimension(:), allocatable :: wpedsclrp

  end type core_sfc_t

  !-------------------------
  !-------------------------
  type clubb_misc_t

    real(core_rknd) :: varmu

    real(core_rknd), dimension(:), allocatable :: qrl_zt
    real(core_rknd), dimension(:), allocatable :: prer_evap
    real(core_rknd), dimension(:), allocatable :: rfrzm      ! Total ice-phase water mixing ratio        [kg/kg] 

  end type clubb_misc_t

contains

  !------------------------------------------------------ 
  subroutine clubb_core_fld_alloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc, pverp, ntracer )

    type(core_auxil_t),  intent(inout) :: core_auxil
    type(core_prog_t),   intent(inout) :: core_prog
    type(core_diag_t),   intent(inout) :: core_diag
    type(core_forcing_t),intent(inout) :: core_forcing
    type(core_sfc_t),    intent(inout) :: core_sfc
    type(clubb_misc_t),  intent(inout) :: clubb_misc

    integer,           intent(in)    :: pverp
    integer,           intent(in)    :: ntracer

    character(len=100) :: routine='clubb_core_fld_alloc'
    integer :: ierr

    !--------------------------------------------------------------------
    ! Auxiliary variables describing the grid-box mean atmospheric state
    !--------------------------------------------------------------------
    allocate( core_auxil% p_in_Pa         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% exner           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_auxil% rho_zt          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% rho_zm          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% rho_ds_zt       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% rho_ds_zm       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% invrs_rho_ds_zm (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% invrs_rho_ds_zt (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_auxil% thv_ds_zt       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% thv_ds_zm       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_auxil% wm_zt           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_auxil% wm_zm           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !------------------------------
    ! Misc.
    !------------------------------
    allocate( clubb_misc% qrl_zt          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( clubb_misc% prer_evap       (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( clubb_misc% rfrzm           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !------------------------------
    ! CLUBB's prognostic variables
    !------------------------------
    allocate( core_prog% um              (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% vm              (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_prog% thlm            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% rtm             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_prog% up2             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% vp2             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% wp2             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% wp3             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_prog% upwp            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% vpwp            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_prog% rtp2            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% thlp2           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% wprtp           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% wpthlp          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_prog% rtpthlp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_prog% edsclr  (pverp,ntracer), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !------------------------------
    ! CLUBB's diagnostic variables
    !------------------------------
    allocate( core_diag% wp2thvp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% wpthvp          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% rtpthvp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% thlpthvp        (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_diag% wprcp           (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% thlprcp         (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_diag% khzm            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% khzt            (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    allocate( core_diag% rcm             (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% cloud_frac      (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% qclvar          (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% rcm_in_layer    (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    allocate( core_diag% cloud_cover     (pverp), stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

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
  subroutine clubb_core_fld_dealloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc )

    type(core_auxil_t),  intent(inout) :: core_auxil
    type(core_prog_t),   intent(inout) :: core_prog
    type(core_diag_t),   intent(inout) :: core_diag
    type(core_forcing_t),intent(inout) :: core_forcing
    type(core_sfc_t),    intent(inout) :: core_sfc
    type(clubb_misc_t),  intent(inout) :: clubb_misc

    character(len=100) :: routine='clubb_core_fld_dealloc'
    integer :: ierr

    !--------------------------------------------------------------------
    ! Auxiliary variables describing the grid-box mean atmospheric state
    !--------------------------------------------------------------------
    deallocate( core_auxil% p_in_Pa        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% exner          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_auxil% rho_zt         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% rho_zm         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% rho_ds_zt      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% rho_ds_zm      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% invrs_rho_ds_zm, stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% invrs_rho_ds_zt, stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_auxil% thv_ds_zt      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% thv_ds_zm      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_auxil% wm_zt          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_auxil% wm_zm          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !------------------------------
    ! Misc.
    !------------------------------
    deallocate( clubb_misc% qrl_zt         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( clubb_misc% prer_evap      , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( clubb_misc% rfrzm          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !------------------------------
    ! CLUBB's prognostic variables
    !------------------------------
    deallocate( core_prog% um             , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% vm             , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_prog% thlm           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% rtm            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_prog% up2            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% vp2            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% wp2            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% wp3            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_prog% upwp           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% vpwp           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_prog% rtp2           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% thlp2          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% wprtp          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% wpthlp         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_prog% rtpthlp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_prog% edsclr         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    !------------------------------
    ! CLUBB's diagnostic variables
    !------------------------------
    deallocate( core_diag% wp2thvp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% wpthvp         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% rtpthvp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% thlpthvp       , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_diag% wprcp          , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% thlprcp        , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_diag% khzm           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% khzt           , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

    deallocate( core_diag% rcm            , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% cloud_frac     , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% qclvar         , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% rcm_in_layer   , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))
    deallocate( core_diag% cloud_cover    , stat=ierr ); if (ierr/=0) call endrun('error in '//trim(routine))

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
