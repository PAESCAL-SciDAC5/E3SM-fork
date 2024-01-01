module clubb_intr

  !----------------------------------------------------------------------------------------------------- !
  ! Module to interface CAM with Cloud Layers Unified by Bi-normals (CLUBB), developed                   !
  !    by the University of Wisconsin Milwaukee Group (UWM).                                             !
  !                                                                                                      !
  ! CLUBB replaces the exisiting turbulence, shallow convection, and macrophysics in CAM5                !
  !                                                                                                      !
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
  ! differencing the diffused and initial states.                                                        !
  !                                                                                                      !
  ! Calling sequence:                                                                                    !
  !                                                                                                      !
  !---------------------------Code history-------------------------------------------------------------- !
  ! Authors:  P. Bogenschutz, C. Craig, A. Gettelman                                                     !
  !                                                                                                      !
  !----------------------------------------------------------------------------------------------------- !
  !  2020-01  O. Guba Correct energy density function
  !-----------------------------------------------------------------------------------------------------
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use shr_log_mod ,  only: errMsg => shr_log_errMsg
  use ppgrid,        only: pcols, pver, pverp
  use phys_control,  only: phys_getopts
  use physconst,     only: rair, cpair, gravit, latvap, latice, zvir, rh2o, karman, &
                           tms_orocnst, tms_z0fac, pi
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc
  use constituents,  only: pcnst, cnst_add
  use pbl_utils,     only: calc_ustar, calc_obklen
  use perf_mod,      only: t_startf, t_stopf
  use mpishorthand
  use cam_abortutils,  only: endrun
  use cam_history_support, only: fillvalue
  use cam_history,   only: outfld
#ifdef CLUBB_SGS
  use clubb_api_module, only: pdf_parameter, clubb_fatal_error, fstderr
  use clubb_precision,  only: core_rknd
#else
  use shr_kind_mod,     only: core_rknd=>shr_kind_r8
#endif

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: clubb_ini_cam, clubb_register_cam, clubb_tend_cam, &
#ifdef CLUBB_SGS
            ! This utilizes CLUBB specific variables in its interface
            stats_init_clubb, &
#endif
            stats_end_timestep_clubb, &
            clubb_surface, &
            clubb_readnl, &
            clubb_init_cnst, &
            clubb_implements_cnst

#ifdef CLUBB_SGS
  ! Both of these utilize CLUBB specific variables in their interface
  private :: stats_zero, stats_avg
#endif

  logical, public :: do_cldcool

  ! ------------ !
  ! Private data !
  ! ------------ !

  integer, parameter :: &
      grid_type    = 3, &               ! The 2 option specifies stretched thermodynamic levels
      hydromet_dim = 0                  ! The hydromet array in SAM-CLUBB is currently 0 elements

  real(core_rknd), dimension(0) :: &
      sclr_tol = 1.e-8_core_rknd        ! Total water in kg/kg

  character(len=6), parameter :: &
      saturation_equation = "flatau"    ! Flatau polynomial approximation for SVP

  real(core_rknd), parameter :: &
      theta0   = 300._core_rknd, &      ! Reference temperature                     [K]
      ts_nudge = 86400._core_rknd, &    ! Time scale for u/v nudging (not used)     [s]
      p0_clubb = 100000._core_rknd

  integer, parameter :: &
    sclr_dim = 0                        ! Higher-order scalars, set to zero

  real(r8), parameter :: &
    wp3_const = 1._r8                   ! Constant to add to wp3 when moments are advected

  real(r8), parameter :: &
    wpthlp_const = 10.0_r8              ! Constant to add to wpthlp when moments are advected

  real(r8), parameter :: &
    wprtp_const = 0.01_r8               ! Constant to add to wprtp when moments are advected

  real(r8), parameter :: &
    rtpthlp_const = 0.01_r8             ! Constant to add to rtpthlp when moments are advected

  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  !PMA
  real(r8), parameter :: qsmall = 1.e-18_r8 ! qsmall used in MG

  real(r8), parameter ::  rad_to_deg = 180.0_r8/pi !converts radians to degrees

  real(r8) :: clubb_timestep = unset_r8  ! Default CLUBB timestep, unless overwriten by namelist
  real(r8) :: clubb_rnevap_effic = unset_r8

  !namelist variables
  real(r8) :: clubb_liq_deep = unset_r8
  real(r8) :: clubb_liq_sh   = unset_r8
  real(r8) :: clubb_ice_deep = unset_r8
  real(r8) :: clubb_ice_sh   = unset_r8
  real(r8) :: clubb_tk1      = unset_r8
  real(r8) :: clubb_tk2      = unset_r8

!  Constant parameters
  logical, parameter, private :: &
    l_uv_nudge       = .false.,       &  ! Use u/v nudging (not used)
    l_input_fields   = .false.,       &  ! Always false for EAM-CLUBB.
    l_implemented    = .true.,        &  ! Implemented in a host model (always true)
    l_host_applies_sfc_fluxes = .false.  ! Whether the host model applies the surface fluxes

  logical            :: do_tms
  logical            :: linearize_pbl_winds
  logical            :: lq(pcnst)
  logical            :: prog_modal_aero
  logical            :: do_rainturb
  logical            :: do_expldiff
  logical            :: clubb_do_adv
  logical            :: clubb_do_deep
  logical            :: micro_do_icesupersat
  logical            :: history_budget
  logical            :: use_sgv !PMA This flag controls tuning for tpert and gustiness
  integer            :: history_budget_histfile_num
  integer            :: edsclr_dim       ! Number of scalars to transport in CLUBB
  integer            :: offset

!  define physics buffer indicies here
  integer :: &
    wp2_idx, &          ! vertical velocity variances
    wp3_idx, &          ! third moment of vertical velocity
    wpthlp_idx, &       ! turbulent flux of thetal
    wprtp_idx, &        ! turbulent flux of total water
    rtpthlp_idx, &      ! covariance of thetal and rt
    rtp2_idx, &         ! variance of total water
    thlp2_idx, &        ! variance of thetal
    up2_idx, &          ! variance of east-west wind
    vp2_idx, &          ! variance of north-south wind
    upwp_idx, &         ! east-west momentum flux
    vpwp_idx, &         ! north-south momentum flux
    um_pert_idx, &      ! perturbed east-west momentum flux
    vm_pert_idx, &      ! perturbed north-south momentum flux
    upwp_pert_idx, &    ! perturbed east-west momentum flux
    vpwp_pert_idx, &    ! perturbed north-south momentum flux
    thlm_idx, &         ! mean thetal
    rtm_idx, &          ! mean total water mixing ratio
    um_idx, &           ! mean of east-west wind
    vm_idx, &           ! mean of north-south wind
    cld_idx, &          ! Cloud fraction
    concld_idx, &       ! Convective cloud fraction
    ast_idx, &          ! Stratiform cloud fraction
    alst_idx, &         ! Liquid stratiform cloud fraction
    aist_idx, &         ! Ice stratiform cloud fraction
    qlst_idx, &         ! Physical in-cloud LWC
    qist_idx, &         ! Physical in-cloud IWC
    dp_frac_idx, &      ! deep convection cloud fraction
    sh_frac_idx, &      ! shallow convection cloud fraction
    rel_idx, &          ! Rel
    kvh_idx, &          ! CLUBB eddy diffusivity on thermo levels
    kvm_idx, &          ! CLUBB eddy diffusivity on mom levels
    pblh_idx, &         ! PBL pbuf
    icwmrdp_idx, &      ! In cloud mixing ratio for deep convection
    tke_idx, &          ! turbulent kinetic energy
    tpert_idx, &        ! temperature perturbation from PBL
    fice_idx, &         ! fice_idx index in physics buffer
    cmeliq_idx, &       ! cmeliq_idx index in physics buffer
    relvar_idx, &       ! relative cloud water variance
    accre_enhan_idx, &  ! optional accretion enhancement factor for MG
    prer_evap_idx, &    ! rain evaporation rate
    qrl_idx, &          ! longwave cooling rate
    radf_idx

 integer :: &          ! newly added pbuf fields for CLUBB
    wpthvp_idx, &       ! < w'th_v' >
    wp2thvp_idx, &      ! < w'^2 th_v' >
    rtpthvp_idx, &      ! < r_t'th_v' >
    thlpthvp_idx, &     ! < th_l'th_v' >
    rcm_idx, &          ! Cloud water mixing ratio
    cloud_frac_idx !, & ! Cloud fraction

  integer :: &          ! added pbuf fields for clubb to have restart bfb when ipdf_call_placement=2
    pdf_zm_w_1_idx, &
    pdf_zm_w_2_idx, &
    pdf_zm_varnce_w_1_idx, &
    pdf_zm_varnce_w_2_idx, &
    pdf_zm_mixt_frac_idx

  integer :: &          !PMA adds pbuf fields for ZM gustiness
    prec_dp_idx, &
    snow_dp_idx, &
    vmag_gust_idx, &
    wsresp_idx, &
    tau_est_idx

  integer, public :: &
    ixthlp2 = 0, &
    ixwpthlp = 0, &
    ixwprtp = 0, &
    ixwp2 = 0, &
    ixwp3 = 0, &
    ixrtpthlp = 0, &
    ixrtp2 = 0, &
    ixup2 = 0, &
    ixvp2 = 0

  integer :: cmfmc_sh_idx = 0

  real(r8) :: dp1 !set in namelist; assigned in cloud_fraction.F90
  !  Output arrays for CLUBB statistics
  real(r8), allocatable, dimension(:,:,:) :: out_zt, out_zm, out_radzt, out_radzm, out_sfc

  character(len=16)  :: eddy_scheme      ! Default set in phys_control.F90
  character(len=16)  :: deep_scheme      ! Default set in phys_control.F90

  integer, parameter :: ncnst=9
  character(len=8)   :: cnst_names(ncnst)
  logical            :: do_cnst=.false.

#ifdef CLUBB_SGS
  type(pdf_parameter), target, allocatable :: pdf_params_chnk(:,:)    ! PDF parameters (thermo. levs.) [units vary]
  type(pdf_parameter), target, allocatable :: pdf_params_zm_chnk(:,:) ! PDF parameters on momentum levs. [units vary]
#endif

  logical :: liqcf_fix = .FALSE.  ! HW for liquid cloud fraction fix
  logical :: relvar_fix = .FALSE. !PMA for relvar fix

  real(r8) :: micro_mg_accre_enhan_fac = huge(1.0_r8) !Accretion enhancement factor from namelist

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_register_cam( )
!-------------------------------------------------------------------------------
! Description:
!   Register the constituents and fields in the physics buffer
! Author: P. Bogenschutz, C. Craig, A. Gettelman
!
!-------------------------------------------------------------------------------
#ifdef CLUBB_SGS

    !------------------------------------------------ !
    ! Register physics buffer fields and constituents !
    !------------------------------------------------ !

    !  Add CLUBB fields to pbuf
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dyn_time_lvls
    use ppgrid,          only: pver, pverp, pcols

    call phys_getopts( eddy_scheme_out                 = eddy_scheme, &
                       deep_scheme_out                 = deep_scheme, &
                       do_tms_out                      = do_tms,      &
                       linearize_pbl_winds_out         = linearize_pbl_winds,      &
                       history_budget_out              = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       micro_do_icesupersat_out        = micro_do_icesupersat, &
                       micro_mg_accre_enhan_fac_out    = micro_mg_accre_enhan_fac)

    if (clubb_do_adv) then
       cnst_names =(/'THLP2  ','RTP2   ','RTPTHLP','WPTHLP ','WPRTP  ','WP2    ','WP3    ','UP2    ','VP2    '/)
       do_cnst=.true.
       !  If CLUBB moments are advected, do not output them automatically which is typically done.  Some moments
       !    need a constant added to them before they are advected, thus this would corrupt the output.
       !    Users should refer to the "XXXX_CLUBB" (THLP2_CLUBB for instance) output variables for these moments
       call cnst_add(trim(cnst_names(1)),0._r8,0._r8,0._r8,ixthlp2,longname='second moment vertical velocity',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(2)),0._r8,0._r8,0._r8,ixrtp2,longname='second moment rtp',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(3)),0._r8,0._r8,-999999._r8,ixrtpthlp,longname='covariance rtp thlp',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(4)),0._r8,0._r8,-999999._r8,ixwpthlp,longname='CLUBB heat flux',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(5)),0._r8,0._r8,-999999._r8,ixwprtp,longname='CLUBB moisture flux',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(6)),0._r8,0._r8,0._r8,ixwp2,longname='CLUBB wp2',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(7)),0._r8,0._r8,-999999._r8,ixwp3,longname='CLUBB 3rd moment vert velocity',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(8)),0._r8,0._r8,0._r8,ixup2,longname='CLUBB 2nd moment u wind',cam_outfld=.false.)
       call cnst_add(trim(cnst_names(9)),0._r8,0._r8,0._r8,ixvp2,longname='CLUBB 2nd moment v wind',cam_outfld=.false.)
    end if

    !  put pbuf_add calls here (see macrop_driver.F90 for sample) use indicies defined at top
    call pbuf_add_field('pblh',       'global', dtype_r8, (/pcols/),                    pblh_idx)
    call pbuf_add_field('tke',        'global', dtype_r8, (/pcols, pverp/),             tke_idx)
    call pbuf_add_field('kvh',        'global', dtype_r8, (/pcols, pverp/),             kvh_idx)
    call pbuf_add_field('kvm',        'global', dtype_r8, (/pcols, pverp/),             kvm_idx)
    call pbuf_add_field('tpert',      'global', dtype_r8, (/pcols/),                    tpert_idx)
    call pbuf_add_field('AST',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    ast_idx)
    call pbuf_add_field('AIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    aist_idx)
    call pbuf_add_field('ALST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    alst_idx)
    call pbuf_add_field('QIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qist_idx)
    call pbuf_add_field('QLST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qlst_idx)
    call pbuf_add_field('CONCLD',     'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    concld_idx)
    call pbuf_add_field('CLD',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    cld_idx)
    call pbuf_add_field('FICE',       'physpkg',dtype_r8, (/pcols,pver/),               fice_idx)
    call pbuf_add_field('RAD_CLUBB',  'global', dtype_r8, (/pcols,pver/),               radf_idx)
    call pbuf_add_field('CMELIQ',     'physpkg',dtype_r8, (/pcols,pver/),                  cmeliq_idx)

    call pbuf_add_field('WP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp2_idx)
    call pbuf_add_field('WP3_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp3_idx)
    call pbuf_add_field('WPTHLP_nadv',     'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wpthlp_idx)
    call pbuf_add_field('WPRTP_nadv',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wprtp_idx)
    call pbuf_add_field('RTPTHLP_nadv',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtpthlp_idx)
    call pbuf_add_field('RTP2_nadv',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtp2_idx)
    call pbuf_add_field('THLP2_nadv',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlp2_idx)
    call pbuf_add_field('UP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), up2_idx)
    call pbuf_add_field('VP2_nadv',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vp2_idx)

    call pbuf_add_field('UPWP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), upwp_idx)
    call pbuf_add_field('VPWP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vpwp_idx)
    if (linearize_pbl_winds) then
       call pbuf_add_field('UM_PERT',    'physpkg', dtype_r8, (/pcols,pverp/), um_pert_idx)
       call pbuf_add_field('VM_PERT',    'physpkg', dtype_r8, (/pcols,pverp/), vm_pert_idx)
       call pbuf_add_field('UPWP_PERT',  'physpkg', dtype_r8, (/pcols,pverp/), upwp_pert_idx)
       call pbuf_add_field('VPWP_PERT',  'physpkg', dtype_r8, (/pcols,pverp/), vpwp_pert_idx)
    end if
    call pbuf_add_field('THLM',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlm_idx)
    call pbuf_add_field('RTM',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtm_idx)
    call pbuf_add_field('UM',         'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), um_idx)
    call pbuf_add_field('VM',         'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), vm_idx)

    call pbuf_add_field('WPTHVP',        'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wpthvp_idx)
    call pbuf_add_field('WP2THVP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wp2thvp_idx)
    call pbuf_add_field('RTPTHVP',       'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rtpthvp_idx)
    call pbuf_add_field('THLPTHVP',      'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), thlpthvp_idx)
    call pbuf_add_field('RCM',           'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), rcm_idx)
    call pbuf_add_field('CLOUD_FRAC',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), cloud_frac_idx)

    call pbuf_add_field('pdf_zm_w_1',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_w_1_idx)
    call pbuf_add_field('pdf_zm_w_2',    'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_w_2_idx)
    call pbuf_add_field('pdf_zm_var_w_1', 'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_varnce_w_1_idx)
    call pbuf_add_field('pdf_zm_var_w_2', 'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_varnce_w_2_idx)
    call pbuf_add_field('pdf_zm_mixt_frac',  'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), pdf_zm_mixt_frac_idx)

    call pbuf_add_field('vmag_gust',       'global', dtype_r8, (/pcols/),      vmag_gust_idx) !PMA total gustiness

    if (linearize_pbl_winds) then
       call pbuf_add_field('wsresp',          'global', dtype_r8, (/pcols/),      wsresp_idx)
       call pbuf_add_field('tau_est',         'global', dtype_r8, (/pcols/),      tau_est_idx)
    end if
#endif

  end subroutine clubb_register_cam
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

function clubb_implements_cnst(name)

  !----------------------------------------------------------------------------- !
  !                                                                              !
  ! Return true if specified constituent is implemented by this package          !
  !                                                                              !
  !----------------------------------------------------------------------------- !

   character(len=*), intent(in) :: name      ! constituent name
   logical :: clubb_implements_cnst     ! return value

   !-----------------------------------------------------------------------

   clubb_implements_cnst = (do_cnst .and. any(name == cnst_names))

end function clubb_implements_cnst


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

subroutine clubb_init_cnst(name, q, gcid)

#ifdef CLUBB_SGS
    use constants_clubb,        only: w_tol_sqd, rt_tol, thl_tol
#endif

   !----------------------------------------------------------------------- !
   !                                                                        !
   ! Initialize the state if clubb_do_adv                                   !
   !                                                                        !
   !----------------------------------------------------------------------- !

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

#ifdef CLUBB_SGS
   if (clubb_do_adv) then
      if (trim(name) == trim(cnst_names(1))) q = real(thl_tol**2, kind = r8)
      if (trim(name) == trim(cnst_names(2))) q = real(rt_tol**2, kind = r8)
      if (trim(name) == trim(cnst_names(3))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(4))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(5))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(6))) q = real(w_tol_sqd, kind = r8)
      if (trim(name) == trim(cnst_names(7))) q = 0.0_r8
      if (trim(name) == trim(cnst_names(8))) q = real(w_tol_sqd, kind = r8)
      if (trim(name) == trim(cnst_names(9))) q = real(w_tol_sqd, kind = r8)
   end if
#endif

end subroutine clubb_init_cnst


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_readnl(nlfile)

#ifdef CLUBB_SGS
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_abortutils,  only: endrun
    use stats_variables, only: l_stats, l_output_rad_files
    use mpishorthand
    use model_flags,     only: l_diffuse_rtm_and_thlm, l_stability_correct_Kh_N2_zm, &
                               l_vert_avg_closure, l_trapezoidal_rule_zt, &
                               l_trapezoidal_rule_zm, l_call_pdf_closure_twice,&
                               ipdf_call_placement

    use parameters_tunable, only: clubb_param_readnl
#endif

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

#ifdef CLUBB_SGS
    logical :: clubb_history, clubb_rad_history, clubb_cloudtop_cooling, clubb_rainevap_turb, &
               clubb_stabcorrect, clubb_expldiff ! Stats enabled (T/F)
    logical :: clubb_use_sgv !PMA This flag controls tuning for tpert and gustiness
    logical :: clubb_vert_avg_closure !XZheng This flag sets four clubb config flags for pdf_closure and the trapezoidal rule to  compute the varibles that are output from high order closure
    integer :: clubb_ipdf_call_placement  !XZheng This flag sets options for the placement of the call to CLUBB's PDF.

    integer :: iunit, read_status

    namelist /clubb_his_nl/ clubb_history, clubb_rad_history
    namelist /clubbpbl_diff_nl/ clubb_cloudtop_cooling, clubb_rainevap_turb, clubb_expldiff, &
                                clubb_do_adv, clubb_do_deep, clubb_timestep, clubb_stabcorrect, &
                                clubb_rnevap_effic, clubb_liq_deep, clubb_liq_sh, clubb_ice_deep, &
                                clubb_ice_sh, clubb_tk1, clubb_tk2, relvar_fix, clubb_use_sgv, &
                                clubb_vert_avg_closure, clubb_ipdf_call_placement


    !----- Begin Code -----

    !  Determine if we want clubb_history to be output
    clubb_history      = .false.   ! Initialize to false
    l_stats            = .false.   ! Initialize to false
    l_output_rad_files = .false.   ! Initialize to false
    do_cldcool         = .false.   ! Initialize to false
    do_rainturb        = .false.   ! Initialize to false
    do_expldiff        = .false.   ! Initialize to false
    relvar_fix         = .false.   ! Initialize to false
    clubb_do_adv       = .false.   ! Initialize to false
    clubb_do_deep      = .false.   ! Initialize to false
    use_sgv            = .false.
    clubb_vert_avg_closure = .true.
    clubb_ipdf_call_placement = -999


    !  Read namelist to determine if CLUBB history should be called
    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(nlfile), status='old' )

      call find_group_name(iunit, 'clubb_his_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_his_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_readnl:  error reading namelist'//errmsg(__FILE__,__LINE__))
         end if
      end if

      call find_group_name(iunit, 'clubbpbl_diff_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubbpbl_diff_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_readnl:  error reading namelist'//errmsg(__FILE__,__LINE__))
         end if
      end if

      close(unit=iunit)
      call freeunit(iunit)
    end if

#ifdef SPMD
! Broadcast namelist variables
      call mpibcast(clubb_history,            1,   mpilog,   0, mpicom)
      call mpibcast(clubb_rad_history,        1,   mpilog,   0, mpicom)
      call mpibcast(clubb_cloudtop_cooling,   1,   mpilog,   0, mpicom)
      call mpibcast(clubb_rainevap_turb,      1,   mpilog,   0, mpicom)
      call mpibcast(clubb_expldiff,           1,   mpilog,   0, mpicom)
      call mpibcast(clubb_do_adv,             1,   mpilog,   0, mpicom)
      call mpibcast(clubb_do_deep,            1,   mpilog,   0, mpicom)
      call mpibcast(clubb_timestep,           1,   mpir8,   0, mpicom)
      call mpibcast(clubb_stabcorrect,        1,   mpilog,   0, mpicom)
      call mpibcast(clubb_rnevap_effic,       1,   mpir8,   0, mpicom)
      call mpibcast(clubb_liq_deep,           1,   mpir8,   0, mpicom)
      call mpibcast(clubb_liq_sh,             1,   mpir8,   0, mpicom)
      call mpibcast(clubb_ice_deep,           1,   mpir8,   0, mpicom)
      call mpibcast(clubb_ice_sh,             1,   mpir8,   0, mpicom)
      call mpibcast(clubb_tk1,                1,   mpir8,   0, mpicom)
      call mpibcast(clubb_tk2,                1,   mpir8,   0, mpicom)
      call mpibcast(relvar_fix,               1,   mpilog,  0, mpicom)
      call mpibcast(clubb_use_sgv,            1,   mpilog,   0, mpicom)
      call mpibcast(clubb_vert_avg_closure,   1,   mpilog,   0, mpicom)
      call mpibcast(clubb_ipdf_call_placement,   1,   mpiint,   0, mpicom)
#endif

    !  Overwrite defaults if they are true
    if (clubb_ipdf_call_placement > 0) ipdf_call_placement = clubb_ipdf_call_placement
    if (clubb_history) l_stats = .true.
    if (clubb_rad_history) l_output_rad_files = .true.
    if (clubb_cloudtop_cooling) do_cldcool = .true.
    if (clubb_rainevap_turb) do_rainturb = .true.
    if (clubb_expldiff) do_expldiff = .true.
    if (clubb_use_sgv) use_sgv =.true.
    if (clubb_stabcorrect .and. clubb_expldiff)  then
      call endrun('clubb_readnl: clubb_stabcorrect and clubb_expldiff may not both be set to true at the same time'//errmsg(__FILE__,__LINE__))
    end if

    if (clubb_stabcorrect) then
      l_diffuse_rtm_and_thlm       = .true.   ! CLUBB flag set to true
      l_stability_correct_Kh_N2_zm = .true.   ! CLUBB flag set to true
    endif

    if (clubb_vert_avg_closure) then
      l_vert_avg_closure       = .true.   ! CLUBB flag set to true
      l_trapezoidal_rule_zt    = .true.   ! CLUBB flag set to true
      l_trapezoidal_rule_zm    = .true.   ! CLUBB flag set to true
      l_call_pdf_closure_twice = .true.   ! CLUBB flag set to true
    else
      l_vert_avg_closure       = .false.   ! CLUBB flag set to false
      l_trapezoidal_rule_zt    = .false.   ! CLUBB flag set to false
      l_trapezoidal_rule_zm    = .false.   ! CLUBB flag set to false
      l_call_pdf_closure_twice = .false.   ! CLUBB flag set to false
    endif

    ! read tunable parameters from namelist, handlings of masterproc vs others
    ! are done within clubb_param_readnl
    call clubb_param_readnl(nlfile)
#endif
  end subroutine clubb_readnl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_ini_cam(pbuf2d, dp1_in)
!-------------------------------------------------------------------------------
! Description:
!   Initialize UWM CLUBB.
! Author: Cheryl Craig March 2011
! Modifications: Pete Bogenschutz 2011 March and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------



#ifdef CLUBB_SGS

    !  From CAM libraries
    use physics_types,          only: physics_state, physics_ptend
    use cam_history,            only: addfld, horiz_only, add_default
    use ppgrid,                 only: pver, pverp, pcols, begchunk, endchunk
    use ref_pres,               only: pref_mid
    use hb_diff,                only: init_hb_diff
    use trb_mtn_stress,         only: init_tms
    use rad_constituents,       only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx

    !  From the CLUBB libraries

    !  From the CLUBB libraries
    use clubb_api_module, only: &
         setup_clubb_core_api, &
         time_precision, &
         core_rknd, &
         set_clubb_debug_level_api, &
         nparams, &
         read_parameters_api, &
         l_stats, &
         l_stats_samp, &
         l_grads, &
         stats_zt, &
         stats_zm, &
         stats_sfc, &
         stats_rad_zt, &
         stats_rad_zm, &
         w_tol_sqd, &
         rt_tol, &
         l_do_expldiff_rtm_thlm, &
         init_pdf_params_api
    use stats_variables,           only: l_output_rad_files

    use units,                     only: getunit, freeunit
    use error_messages,            only: handle_errmsg
    use time_manager,              only: is_first_step
    use constants_clubb,           only: thl_tol


    !  These are only needed if we're using a passive scalar
    use array_index,            only: iisclr_rt, iisclr_thl, iisclr_CO2, &    ! [kg/kg]/[K]/[1e6 mol/mol]
                                      iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2 ! "    "
    use constituents,           only: cnst_get_ind
    use phys_control,           only: phys_getopts

    use parameters_tunable, only: params_list

#endif

    use physics_buffer,         only: pbuf_get_index, pbuf_set_field, physics_buffer_desc
    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    real(r8) :: dp1_in

#ifdef CLUBB_SGS

    real(kind=time_precision) :: dum1, dum2, dum3

    real(core_rknd), dimension(nparams)  :: clubb_params    ! These adjustable CLUBB parameters (C1, C2 ...)

    logical :: clubb_history, clubb_rad_history, clubb_cloudtop_cooling, clubb_rainevap_turb, clubb_expldiff ! Stats enabled (T/F)

    ! The similar name to clubb_history is unfortunate...
    logical :: history_amwg, history_clubb

    character(len=128) :: errstring             ! error status for CLUBB init

    integer :: err_code, iunit                  ! Code for when CLUBB fails
    integer :: i, j, k, l, idx_chunk, idx_pcols ! Indices
    integer :: read_status                      ! Length of a string
    integer :: ntop_eddy                        ! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
    integer :: nbot_eddy                        ! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
    integer :: nmodes, nspec, pmam_ncnst, m
    integer :: ixnumliq
    integer :: lptr

    real(core_rknd)  :: zt_g(pverp)                        ! Height dummy array
    real(core_rknd)  :: zi_g(pverp)                        ! Height dummy array


    !----- Begin Code -----
    !$OMP PARALLEL
    l_do_expldiff_rtm_thlm = do_expldiff
    !$OMP END PARALLEL

    allocate( &
       pdf_params_chnk(pcols,begchunk:endchunk),   &
       pdf_params_zm_chnk(pcols,begchunk:endchunk) )

    do idx_chunk = begchunk, endchunk
        do idx_pcols = 1, pcols
            call init_pdf_params_api( pverp, pdf_params_chnk(idx_pcols,idx_chunk) )
            call init_pdf_params_api( pverp, pdf_params_zm_chnk(idx_pcols,idx_chunk) )
        end do
    end do

    ! ----------------------------------------------------------------- !
    ! Determine how many constituents CLUBB will transport.  Note that
    ! CLUBB does not transport aerosol consituents.  Therefore, need to
    ! determine how many aerosols constituents there are and subtract that
    ! off of pcnst (the total consituents)
    ! ----------------------------------------------------------------- !

    call phys_getopts(prog_modal_aero_out=prog_modal_aero, &
                      history_amwg_out=history_amwg, &
                      history_clubb_out=history_clubb,&
                      liqcf_fix_out   = liqcf_fix)

    !  Select variables to apply tendencies back to CAM

    ! Initialize all consituents to true to start
    lq(1:pcnst) = .true.
    edsclr_dim  = pcnst

    if (prog_modal_aero) then
       ! Turn off modal aerosols and decrement edsclr_dim accordingly
       call rad_cnst_get_info(0, nmodes=nmodes)

       do m = 1, nmodes
          call rad_cnst_get_mode_num_idx(m, lptr)
          lq(lptr)=.false.
          edsclr_dim = edsclr_dim-1

          call rad_cnst_get_info(0, m, nspec=nspec)
          do l = 1, nspec
             call rad_cnst_get_mam_mmr_idx(m, l, lptr)
             lq(lptr)=.false.
             edsclr_dim = edsclr_dim-1
          end do
       end do

       !  In addition, if running with MAM, droplet number is transported
       !  in dropmixnuc, therefore we do NOT want CLUBB to apply transport
       !  tendencies to avoid double counted.  Else, we apply tendencies.
       call cnst_get_ind('NUMLIQ',ixnumliq)
       lq(ixnumliq) = .false.
       edsclr_dim = edsclr_dim-1
    endif

    ! ----------------------------------------------------------------- !
    ! Set the debug level.  Level 2 has additional computational expense since
    ! it checks the array variables in CLUBB for invalid values.
    ! ----------------------------------------------------------------- !
    call set_clubb_debug_level_api( 0 )

    ! ----------------------------------------------------------------- !
    ! use pbuf_get_fld_idx to get existing physics buffer fields from other
    ! physics packages (e.g. tke)
    ! ----------------------------------------------------------------- !


    !  Defaults
    l_stats_samp = .false.
    l_grads = .false.

    !  Overwrite defaults if needbe
    if (l_stats) l_stats_samp = .true.

    !  Define physics buffers indexes
    cld_idx     = pbuf_get_index('CLD')         ! Cloud fraction
    concld_idx  = pbuf_get_index('CONCLD')      ! Convective cloud cover
    ast_idx     = pbuf_get_index('AST')         ! Stratiform cloud fraction
    alst_idx    = pbuf_get_index('ALST')        ! Liquid stratiform cloud fraction
    aist_idx    = pbuf_get_index('AIST')        ! Ice stratiform cloud fraction
    qlst_idx    = pbuf_get_index('QLST')        ! Physical in-stratus LWC
    qist_idx    = pbuf_get_index('QIST')        ! Physical in-stratus IWC
    dp_frac_idx = pbuf_get_index('DP_FRAC')     ! Deep convection cloud fraction
    icwmrdp_idx = pbuf_get_index('ICWMRDP')     ! In-cloud deep convective mixing ratio
    sh_frac_idx = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx  = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    accre_enhan_idx = pbuf_get_index('ACCRE_ENHAN') ! accretion enhancement for MG
    prer_evap_idx   = pbuf_get_index('PRER_EVAP')
    qrl_idx         = pbuf_get_index('QRL')
    cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')

    prec_dp_idx = pbuf_get_index('PREC_DP') !PMA ZM precip for gustiness
    snow_dp_idx = pbuf_get_index('SNOW_DP') !PMA ZM snow for gustiness
    vmag_gust_idx = pbuf_get_index('vmag_gust') !PMA ZM snow for gustiness

    iisclr_rt  = -1
    iisclr_thl = -1
    iisclr_CO2 = -1

    iiedsclr_rt  = -1
    iiedsclr_thl = -1
    iiedsclr_CO2 = -1

    ! ----------------------------------------------------------------- !
    ! Define number of tracers for CLUBB to diffuse
    ! ----------------------------------------------------------------- !

    if (do_expldiff) then
       offset = 2 ! diffuse temperature and moisture explicitly
       edsclr_dim = edsclr_dim + offset
    endif

    ! ----------------------------------------------------------------- !
    ! Setup CLUBB core
    ! ----------------------------------------------------------------- !

    !  Read in parameters for CLUBB.  Pack the default and updated (via nml)
    !  tunable parameters into clubb_params
!$OMP PARALLEL
    call read_parameters_api( -99, "", clubb_params )
!$OMP END PARALLEL

    ! Print the list of CLUBB parameters, if multi-threaded, it may print by each thread
    if (masterproc) then
       write(iulog,*)'CLUBB tunable parameters: total ',nparams
       write(iulog,*)'--------------------------------------------------'
       do i = 1, nparams
          write(iulog,*) params_list(i), " = ", clubb_params(i)
       enddo
    endif


    !  Fill in dummy arrays for height.  Note that these are overwrote
    !  at every CLUBB step to physical values.
    do k=1,pverp
       zt_g(k) = ((k-1)*1000._core_rknd)-500._core_rknd  !  this is dummy garbage
       zi_g(k) = (k-1)*1000._core_rknd                   !  this is dummy garbage
    enddo

    !  Set up CLUBB core.  Note that some of these inputs are overwrote
    !  when clubb_tend_cam is called.  The reason is that heights can change
    !  at each time step, which is why dummy arrays are read in here for heights
    !  as they are immediately overwrote.
!$OMP PARALLEL
    call setup_clubb_core_api     &
         ( pverp, theta0, ts_nudge, &                                 ! In
           hydromet_dim,  sclr_dim, &                                 ! In
           sclr_tol, edsclr_dim, clubb_params, &                      ! In
           l_host_applies_sfc_fluxes, &                               ! In
           l_uv_nudge, saturation_equation, l_input_fields,  &        ! In
           l_implemented, grid_type, zi_g(2), zi_g(1), zi_g(pverp), & ! In
           zi_g(1:pverp), zt_g(1:pverp), zi_g(1), &
           err_code )
!$OMP END PARALLEL

    ! ----------------------------------------------------------------- !
    ! Set-up HB diffusion.  Only initialized to diagnose PBL depth      !
    ! ----------------------------------------------------------------- !

    ! Initialize eddy diffusivity module

    ntop_eddy = 1    ! if >1, must be <= nbot_molec
    nbot_eddy = pver ! currently always pver

    call init_hb_diff( gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, karman, eddy_scheme )

    ! ----------------------------------------------------------------- !
    ! Initialize turbulent mountain stress module                       !
    ! ------------------------------------------------------------------!

    if ( do_tms) then
       call init_tms( r8, tms_orocnst, tms_z0fac, karman, gravit, rair, errstring)
       call handle_errmsg(errstring, subname="init_tms")

       call addfld( 'TAUTMSX' ,  horiz_only,  'A','N/m2',  'Zonal      turbulent mountain surface stress' )
       call addfld( 'TAUTMSY' ,  horiz_only,  'A','N/m2',  'Meridional turbulent mountain surface stress' )
       if (history_amwg) then
          call add_default( 'TAUTMSX ', 1, ' ' )
          call add_default( 'TAUTMSY ', 1, ' ' )
       end if
       if (masterproc) then
          write(iulog,*)'Using turbulent mountain stress module'
          write(iulog,*)'  tms_orocnst = ',tms_orocnst
          write(iulog,*)'  tms_z0fac = ',tms_z0fac
       end if
    endif

    ! ----------------------------------------------------------------- !
    ! Add output fields for the history files
    ! ----------------------------------------------------------------- !

    if (clubb_do_deep) then
       call addfld ('MU_CLUBB',horiz_only,'A','1/m','CLUBB value of entrainment')
    endif

    !  These are default CLUBB output.  Not the higher order history budgets
    call addfld ('RHO_CLUBB',    (/ 'ilev' /), 'A',        'kg/m3', 'Air Density')
    call addfld ('UP2_CLUBB',    (/ 'ilev' /), 'A',        'm2/s2', 'Zonal Velocity Variance')
    call addfld ('VP2_CLUBB',    (/ 'ilev' /), 'A',        'm2/s2', 'Meridional Velocity Variance')
    call addfld ('WP2_CLUBB',    (/ 'ilev' /), 'A',        'm2/s2', 'Vertical Velocity Variance')
    call addfld ('UPWP_CLUBB',    (/ 'ilev' /), 'A',       'm2/s2', 'Zonal Momentum Flux')
    call addfld ('VPWP_CLUBB',    (/ 'ilev' /), 'A',       'm2/s2', 'Meridional Momentum Flux')
    call addfld ('WP3_CLUBB',    (/ 'ilev' /), 'A',        'm3/s3', 'Third Moment Vertical Velocity')
    call addfld ('WPTHLP_CLUBB',     (/ 'ilev' /), 'A',     'W/m2', 'Heat Flux')
    call addfld ('WPRTP_CLUBB',     (/ 'ilev' /), 'A',      'W/m2', 'Moisture Flux')
    call addfld ('RTP2_CLUBB', (/ 'ilev' /), 'A',       'g^2/kg^2', 'Moisture Variance')
    call addfld ('THLP2_CLUBB',      (/ 'ilev' /), 'A',      'K^2', 'Temperature Variance')
    call addfld ('RTPTHLP_CLUBB',   (/ 'ilev' /), 'A',    'K g/kg', 'Temp. Moist. Covariance')
    call addfld ('RCM_CLUBB',     (/ 'ilev' /), 'A',        'g/kg', 'Cloud Water Mixing Ratio')
    call addfld ('WPRCP_CLUBB',     (/ 'ilev' /), 'A',      'W/m2', 'Liquid Water Flux')
    call addfld ('CLOUDFRAC_CLUBB', (/ 'lev' /),  'A',  '1', 'Cloud Fraction')
    call addfld ('RCMINLAYER_CLUBB',     (/ 'ilev' /), 'A', 'g/kg', 'Cloud Water in Layer')
    call addfld ('CLOUDCOVER_CLUBB', (/ 'ilev' /), 'A', '1', 'Cloud Cover')
    call addfld ('WPTHVP_CLUBB',     (/ 'lev' /),  'A',     'W/m2', 'Buoyancy Flux')
    call addfld ('RVMTEND_CLUBB',  (/ 'lev' /),  'A',    'g/kg /s', 'Water vapor tendency')
    call addfld ('TTEND_CLUBB',      (/ 'lev' /),  'A',      'k/s', 'Temperature tendency')
    call addfld ('RCMTEND_CLUBB',  (/ 'lev' /),  'A',    'g/kg /s', 'Cloud Liquid Water Tendency')
    call addfld ('RIMTEND_CLUBB',  (/ 'lev' /),  'A',    'g/kg /s', 'Cloud Ice Tendency')
    call addfld ('UTEND_CLUBB',   (/ 'lev' /),  'A',      'm/s /s', 'U-wind Tendency')
    call addfld ('VTEND_CLUBB',   (/ 'lev' /),  'A',      'm/s /s', 'V-wind Tendency')
    call addfld ('ZT_CLUBB',        (/ 'ilev' /), 'A',         'm', 'Thermodynamic Heights')
    call addfld ('ZM_CLUBB',        (/ 'ilev' /), 'A',         'm', 'Momentum Heights')
    call addfld ('UM_CLUBB',      (/ 'ilev' /), 'A',         'm/s', 'Zonal Wind')
    call addfld ('VM_CLUBB',      (/ 'ilev' /), 'A',         'm/s', 'Meridional Wind')
    call addfld ('THETAL',        (/ 'lev' /),  'A',           'K', 'Liquid Water Potential Temperature')
    call addfld ('PBLH',        horiz_only,     'A',             'm', 'PBL height')
    call addfld ('QT',    (/ 'lev' /),  'A',               'kg/kg', 'Total water mixing ratio')
    call addfld ('SL',     (/ 'lev' /),  'A',               'J/kg', 'Liquid water static energy')
    call addfld ('CLDST', (/ 'lev' /),  'A',            'fraction', 'Stratus cloud fraction')
    call addfld ('ZMDLF',  (/ 'lev' /),  'A',            'kg/kg/s', 'Detrained liquid water from ZM convection')
    call addfld ('TTENDICE',     (/ 'lev' /),  'A',         'K/s', 'T tendency from Ice Saturation Adjustment')
    call addfld ('QVTENDICE', (/ 'lev' /),  'A',        'kg/kg/s', 'Q tendency from Ice Saturation Adjustment')
    call addfld ('QITENDICE', (/ 'lev' /),  'A',        'kg/kg/s', 'CLDICE tendency from Ice Saturation Adjustment')
    call addfld ('NITENDICE', (/ 'lev' /),  'A',        'kg/kg/s', 'NUMICE tendency from Ice Saturation Adjustment')
    call addfld ('DPDLFLIQ', (/ 'lev' /),  'A',        'kg/kg/s', 'Detrained liquid water from deep convection')
    call addfld ('DPDLFICE', (/ 'lev' /),  'A',        'kg/kg/s', 'Detrained ice from deep convection')
    call addfld ('DPDLFT', (/ 'lev' /),  'A',        'K/s', 'T-tendency due to deep convective detrainment')
    call addfld ('RELVAR', (/ 'lev' /),  'A',        '-', 'Relative cloud water variance')
    call addfld ('RELVARC', (/ 'lev' /),  'A',        '-', 'Relative cloud water variance', flag_xyfill=.true.,fill_value=fillvalue)
    call addfld ('CONCLD', (/ 'lev' /),  'A',        'fraction', 'Convective cloud cover')
    call addfld ('CMELIQ', (/ 'lev' /),  'A',        'kg/kg/s', 'Rate of cond-evap of liq within the cloud')
!PMA gustiness output fields
    call addfld ('VMAGGUST',       horiz_only,     'A',             '-', 'Total gustiness enhancement')
    call addfld ('VMAGDP',        horiz_only,     'A',             '-', 'ZM gustiness enhancement')
    call addfld ('VMAGCL',        horiz_only,     'A',             '-', 'CLUBB gustiness enhancement')
    call addfld ('TPERTBLT',        horiz_only,     'A',             'K', 'perturbation temperature at PBL top')

    !  Initialize statistics, below are dummy variables
    dum1 = 300._r8
    dum2 = 1200._r8
    dum3 = 300._r8

    if (l_stats) then

       call stats_init_clubb( .true., dum1, dum2, &
                         pverp, pverp, pverp, dum3 )

       allocate(out_zt(pcols,pverp,stats_zt%num_output_fields))
       allocate(out_zm(pcols,pverp,stats_zm%num_output_fields))
       allocate(out_sfc(pcols,1,stats_sfc%num_output_fields))

       allocate(out_radzt(pcols,pverp,stats_rad_zt%num_output_fields))
       allocate(out_radzm(pcols,pverp,stats_rad_zm%num_output_fields))

    endif

    ! ----------------------------------------------------------------- !
    ! Make all of this output default, this is not CLUBB history
    ! ----------------------------------------------------------------- !
    if (clubb_do_adv .or. history_clubb) then
       call add_default('WP2_CLUBB',        1, ' ')
       call add_default('WP3_CLUBB',        1, ' ')
       call add_default('WPTHLP_CLUBB',     1, ' ')
       call add_default('WPRTP_CLUBB',      1, ' ')
       call add_default('RTP2_CLUBB',       1, ' ')
       call add_default('THLP2_CLUBB',      1, ' ')
       call add_default('RTPTHLP_CLUBB',    1, ' ')
       call add_default('UP2_CLUBB',        1, ' ')
       call add_default('VP2_CLUBB',        1, ' ')
    end if

    if (history_clubb) then

       if (clubb_do_deep) then
          call add_default('MU_CLUBB',         1, ' ')
       endif

       call add_default('RELVAR',           1, ' ')
       call add_default('RHO_CLUBB',        1, ' ')
       call add_default('UPWP_CLUBB',       1, ' ')
       call add_default('VPWP_CLUBB',       1, ' ')
       call add_default('RCM_CLUBB',        1, ' ')
       call add_default('WPRCP_CLUBB',      1, ' ')
       call add_default('CLOUDFRAC_CLUBB',  1, ' ')
       call add_default('RCMINLAYER_CLUBB', 1, ' ')
       call add_default('CLOUDCOVER_CLUBB', 1, ' ')
       call add_default('WPTHVP_CLUBB',     1, ' ')
       call add_default('RVMTEND_CLUBB',    1, ' ')
       call add_default('TTEND_CLUBB',      1, ' ')
       call add_default('RCMTEND_CLUBB',    1, ' ')
       call add_default('RIMTEND_CLUBB',    1, ' ')
       call add_default('UTEND_CLUBB',      1, ' ')
       call add_default('VTEND_CLUBB',      1, ' ')
       call add_default('ZT_CLUBB',         1, ' ')
       call add_default('ZM_CLUBB',         1, ' ')
       call add_default('UM_CLUBB',         1, ' ')
       call add_default('VM_CLUBB',         1, ' ')
       call add_default('SL',               1, ' ')
       call add_default('QT',               1, ' ')
       call add_default('CONCLD',           1, ' ')
    else
       call add_default('CLOUDFRAC_CLUBB',  1, ' ')
       call add_default('CONCLD',           1, ' ')
    end if

    if (history_amwg) then
       call add_default('PBLH',             1, ' ')
    end if

    if (history_budget) then
       call add_default('DPDLFLIQ',         history_budget_histfile_num, ' ')
       call add_default('DPDLFICE',         history_budget_histfile_num, ' ')
       call add_default('DPDLFT',           history_budget_histfile_num, ' ')
       call add_default('TTEND_CLUBB',      history_budget_histfile_num, ' ')
       call add_default('RCMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('RIMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('RVMTEND_CLUBB',    history_budget_histfile_num, ' ')
       call add_default('UTEND_CLUBB',      history_budget_histfile_num, ' ')
       call add_default('VTEND_CLUBB',      history_budget_histfile_num, ' ')
    endif


    ! --------------- !
    ! First step?     !
    ! Initialization  !
    ! --------------- !

    !  Is this the first time step?  If so then initialize CLUBB variables as follows
    if (is_first_step()) then

       call pbuf_set_field(pbuf2d, wp2_idx,     real(w_tol_sqd, kind = r8))
       call pbuf_set_field(pbuf2d, wp3_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, wpthlp_idx,  0.0_r8)
       call pbuf_set_field(pbuf2d, wprtp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthlp_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, rtp2_idx,    real(rt_tol**2, kind = r8))
       call pbuf_set_field(pbuf2d, thlp2_idx,   real(thl_tol**2, kind = r8))
       call pbuf_set_field(pbuf2d, up2_idx,     real(w_tol_sqd, kind = r8))
       call pbuf_set_field(pbuf2d, vp2_idx,     real(w_tol_sqd, kind = r8))

       call pbuf_set_field(pbuf2d, upwp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, vpwp_idx,    0.0_r8)
       if (linearize_pbl_winds) then
          call pbuf_set_field(pbuf2d, um_pert_idx, 0.0_r8)
          call pbuf_set_field(pbuf2d, vm_pert_idx, 0.0_r8)
          call pbuf_set_field(pbuf2d, upwp_pert_idx, 0.0_r8)
          call pbuf_set_field(pbuf2d, vpwp_pert_idx, 0.0_r8)
       end if
       call pbuf_set_field(pbuf2d, tke_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, kvh_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, fice_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, radf_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, qrl_idx,     0.0_r8)

       call pbuf_set_field(pbuf2d, wpthvp_idx,     0.0_r8)
       call pbuf_set_field(pbuf2d, wp2thvp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, rtpthvp_idx,    0.0_r8)
       call pbuf_set_field(pbuf2d, thlpthvp_idx,   0.0_r8)
       call pbuf_set_field(pbuf2d, rcm_idx,        0.0_r8)
       call pbuf_set_field(pbuf2d, cloud_frac_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d, pdf_zm_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_1_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_varnce_w_2_idx, 0.0_r8)
       call pbuf_set_field(pbuf2d, pdf_zm_mixt_frac_idx, 0.0_r8)

       call pbuf_set_field(pbuf2d, vmag_gust_idx,    1.0_r8)

       if (linearize_pbl_winds) then
          call pbuf_set_field(pbuf2d, wsresp_idx,    0.0_r8)
          call pbuf_set_field(pbuf2d, tau_est_idx,   0.0_r8)
       end if
    endif

    ! --------------- !
    ! End             !
    ! Initialization  !
    ! --------------- !

#endif
    dp1 = dp1_in !set via namelist, assigned in cloud_fraction.F90
    end subroutine clubb_ini_cam


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

   subroutine clubb_tend_cam( &
                              state,   ptend_all,   pbuf,     hdtime, &
                              cmfmc,   cam_in,   sgh30, &
                              macmic_it, cld_macmic_num_steps,dlf, det_s, det_ice, alst_o)

!-------------------------------------------------------------------------------
! Description: Provide tendencies of shallow convection, turbulence, and
!              macrophysics from CLUBB to CAM
!
! Author: Cheryl Craig, March 2011
! Modifications: Pete Bogenschutz, March 2011 and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------

   use physics_types,  only: physics_state, physics_ptend, &
                             physics_state_copy, physics_ptend_init, &
                             physics_ptend_sum !, set_dry_to_wet

   use physics_update_mod, only: physics_update

   use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                             pbuf_set_field, physics_buffer_desc

   use ppgrid,         only: pver, pverp, pcols
   use constituents,   only: cnst_get_ind, cnst_type
  !use co2_cycle,      only: co2_cycle_set_cnst_type
   use camsrfexch,     only: cam_in_t
   use ref_pres,       only: top_lev => trop_cloud_top_lev
   use time_manager,   only: is_first_step
   use cam_abortutils, only: endrun
   use wv_saturation,  only: qsat
   use micro_mg_cam,   only: micro_mg_version

#ifdef CLUBB_SGS

   use scamMOD,                   only: single_column,scm_clubb_iop_name
   use cldfrc2m,                  only: aist_vector
   use cam_history,               only: outfld
   use trb_mtn_stress,            only: compute_tms
   use macrop_driver_with_clubb,  only: ice_supersat_adj_tend
   use macrop_driver_with_clubb,  only: pblh_diag
   use macrop_driver_with_clubb,  only: gustiness
   use macrop_driver_with_clubb,  only: deepcu_detrainment_tend

   use clubb_api_module, only: &
        cleanup_clubb_core_api, &
        w_tol_sqd, &
        rt_tol, &
        thl_tol, &
        l_stats, &
        stats_tsamp, &
        stats_tout, &
        stats_zt, &
        stats_sfc, &
        stats_zm, &
        stats_rad_zt, &
        stats_rad_zm, &
        l_output_rad_files, &
        pdf_parameter, &
        stats_begin_timestep_api, &
        zt2zm_api, zm2zt_api

   use clubb_intr_types
   use clubb_intr_core_types, only: core_auxil_t, core_prog_t, core_diag_t, core_forcing_t, core_sfc_t, clubb_misc_t
   use clubb_intr_core_types, only: clubb_core_fld_alloc, clubb_core_fld_dealloc
#endif

   implicit none

   ! --------------- !
   ! Input Auguments !
   ! --------------- !

   type(physics_state), intent(in)    :: state                    ! Physics state variables                 [vary]
   type(cam_in_t),      intent(in)    :: cam_in
   real(r8),            intent(in)    :: hdtime                   ! Host model timestep                     [s]
   real(r8),            intent(in)    :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection [kg/ks/s]
   real(r8),            intent(in)    :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c           [kg/m2/s]
   real(r8),            intent(in)    :: sgh30(pcols)             ! std deviation of orography              [m]
   integer,             intent(in)    :: cld_macmic_num_steps     ! number of mac-mic iterations
   integer,             intent(in)    :: macmic_it                ! number of mac-mic iterations

   ! ---------------------- !
   ! Input-Output Auguments !
   ! ---------------------- !

   type(physics_buffer_desc), pointer :: pbuf(:)

   ! ---------------------- !
   ! Output Auguments !
   ! ---------------------- !

   type(physics_ptend), intent(out)   :: ptend_all                 ! package tendencies

   ! These two variables are needed for energy check
   real(r8),            intent(out)   :: det_s(pcols)              ! Integral of detrained static energy from ice
   real(r8),            intent(out)   :: det_ice(pcols)            ! Integral of detrained ice for energy check

   real(r8), intent(out) :: alst_o(pcols,pver)  ! H. Wang: for old liquid status fraction

   ! --------------- !
   ! Local Variables !
   ! --------------- !

#ifdef CLUBB_SGS

   type(physics_state),target :: state1                ! Local copy of state variable
   type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

   integer :: i, j, k, t, ixind
   integer :: ixcldice, ixcldliq, ixnumliq, ixnumice, ixq
   integer :: itim_old
   integer :: ncol, lchnk                       ! # of columns, and chunk identifier
   integer :: icnt, clubbtop


   integer :: n_clubb_core_step

  !===========================================================================================
  ! The variables defined as core_rknd is required by the advance_clubb_core_api() subroutine
  !===========================================================================================
   integer :: nz
   real(core_rknd) :: dtime                            ! CLUBB time step                              [s]

   !-----------------
   real(core_rknd) :: zt_bot                  ! height of themo level that is closest to the Earth's surface [m]
   real(core_rknd) :: core_rknd_landfrac
   real(core_rknd) :: core_rknd_rnevap_effic

   real(core_rknd) :: fcoriolis                        ! Coriolis forcing                              [s^-1]
   real(core_rknd) :: sfc_elevation                    ! Elevation of ground                           [m AMSL]

   type(core_auxil_t)   :: core_auxil
   type(core_prog_t)    :: core_prog
   type(core_diag_t)    :: core_diag
   type(core_forcing_t) :: core_forcing
   type(core_sfc_t)     :: core_sfc
   type(clubb_misc_t)   :: clubb_misc

   type(pdf_parameter), pointer :: pdf_params    ! PDF parameters (thermo. levs.) [units vary]
   type(pdf_parameter), pointer :: pdf_params_zm ! PDF parameters on momentum levs. [units vary]
   !-----------------

   real(core_rknd) :: C_10                             ! transfer coefficient                          [-]

   real(core_rknd) :: dum_core_rknd                    ! dummy variable  [units vary]
   real(core_rknd) :: hdtime_core_rknd                 ! host model's cloud macmic timestep in core_rknd

   real(core_rknd), pointer :: upwp_sfc_pert    ! u'w' at surface                               [m^2/s^2]
   real(core_rknd), pointer :: vpwp_sfc_pert    ! v'w' at surface                               [m^2/s^2]
   ! Pointers to temporary copies of particular columns of um_pert/upwp_pert fields
   real(core_rknd), pointer, dimension(:) :: um_pert_col ! Pointer to a particular column of um
   real(core_rknd), pointer, dimension(:) :: vm_pert_col
   real(core_rknd), pointer, dimension(:) :: upwp_pert_col
   real(core_rknd), pointer, dimension(:) :: vpwp_pert_col

   !===========================================================================================================================

   real(r8) :: apply_const
   real(r8) :: newfice(pcols,pver)              ! fraction of ice in cloud at CLUBB start       [-]
   real(r8) :: bflx22                           ! Variable for buoyancy flux for pbl            [K m/s]
   real(r8) :: invrs_hdtime                     ! Preculate 1/hdtime to reduce divide operations

  !real(r8) :: invrs_gravit                     ! Preculate 1/gravit to reduce divide operations

   real(r8) :: ubar                             ! surface wind                                  [m/s]
   real(r8) :: ustar                            ! surface stress                                [m/s]
   real(r8) :: z0                               ! roughness height                              [m]
   real(r8) :: zo                               ! roughness height                              [m]
   real(r8) :: dz_g(pver)                       ! thickness of layer                            [m]
   real(r8) :: minqn                            ! minimum total cloud liquid + ice threshold    [kg/kg]
   real(r8) :: tempqn                           ! temporary total cloud liquid + ice            [kg/kg]
   real(r8) :: relvarmax,relvarmin
   real(r8) :: qmin
   real(r8) :: varmu(pcols)
   real(r8) :: zt_out(pcols,pverp)              ! output for the thermo CLUBB grid              [m]
   real(r8) :: zi_out(pcols,pverp)              ! output for momentum CLUBB grid                [m]

   ! Variables below are needed to compute energy integrals for conservation
   real(r8) :: te_b
   real(r8) :: clubb_s(pver)

  !real(r8) :: exner_clubb(pcols,pverp)         ! Exner function consistent with CLUBB          [-]

   real(r8) :: wpthlp_output(pcols,pverp)       ! Heat flux output variable                     [W/m2]
   real(r8) :: wprtp_output(pcols,pverp)        ! Total water flux output variable              [W/m2]
   real(r8) :: wp3_output(pcols,pverp)          ! wp3 output                                    [m^3/s^3]
   real(r8) :: rtpthlp_output(pcols,pverp)      ! rtpthlp ouptut                                [K kg/kg]
   real(r8) :: qt_output(pcols,pver)            ! Total water mixing ratio for output           [kg/kg]
   real(r8) :: thetal_output(pcols,pver)        ! Liquid water potential temperature output     [K]
   real(r8) :: sl_output(pcols,pver)            ! Liquid water static energy                    [J/kg]
   real(r8) :: rho(pcols,pverp)                 ! Midpoint density in CAM                       [kg/m^3]
   real(r8) :: thv(pcols,pver)                        ! virtual potential temperature                 [K]

   real(r8) :: edsclr_out(pcols,pverp,edsclr_dim)     ! Scalars to be diffused through CLUBB          [units vary]

   real(r8) :: rcm_in_layer(pcols,pverp)        ! CLUBB in-cloud liquid water mixing ratio      [kg/kg]
   real(r8) :: cloud_cover(pcols,pverp)         ! CLUBB in-cloud cloud fraction                 [fraction]
   real(r8) :: wprcp(pcols,pverp)               ! CLUBB liquid water flux                       [m/s kg/kg]
   real(r8) :: wpthvp_diag(pcols,pverp)              ! CLUBB buoyancy flux                           [W/m^2]
   real(r8) :: eps                              ! Rv/Rd                                         [-]
   real(r8) :: dum1                             ! dummy variable                                [units vary]

   integer :: kkhost

   real(r8) :: ksrftms(pcols)                   ! Turbulent mountain stress surface drag        [kg/s/m2]
   real(r8) :: tautmsx(pcols)                   ! U component of turbulent mountain stress      [N/m2]
   real(r8) :: tautmsy(pcols)                   ! V component of turbulent mountain stress      [N/m2]

   real(r8) :: tmp_array(state%ncol,pverp)

   integer                               :: time_elapsed                ! time keep track of stats          [s]
   character(len=200)                    :: temp1, sub                  ! Strings needed for CLUBB output
   logical                               :: l_Lscale_plume_centered, l_use_ice_latent
  !character(len=3), dimension(pcnst)    :: cnst_type_loc               ! local override option for constituents cnst_type


   ! --------------- !
   ! Pointers        !
   ! --------------- !

   type(clubb_mean_2d_t) :: host_mean
   type(clubb_mnts_2d_t) :: host_mnts
   type(clubb_to_host_t) :: c2h

  !real(r8), pointer, dimension(:,:) :: wp2      ! vertical velocity variance                   [m^2/s^2]
  !real(r8), pointer, dimension(:,:) :: wp3      ! third moment of vertical velocity            [m^3/s^3]
  !real(r8), pointer, dimension(:,:) :: wpthlp   ! turbulent flux of thetal                     [m/s K]
  !real(r8), pointer, dimension(:,:) :: wprtp    ! turbulent flux of moisture                   [m/s kg/kg]
  !real(r8), pointer, dimension(:,:) :: rtpthlp  ! covariance of thetal and qt                  [kg/kg K]
  !real(r8), pointer, dimension(:,:) :: rtp2     ! moisture variance                            [kg^2/kg^2]
  !real(r8), pointer, dimension(:,:) :: thlp2    ! temperature variance                         [K^2]
  !real(r8), pointer, dimension(:,:) :: up2      ! east-west wind variance                      [m^2/s^2]
  !real(r8), pointer, dimension(:,:) :: vp2      ! north-south wind variance                    [m^2/s^2]

  !real(r8), pointer, dimension(:,:) :: wpthvp     ! < w'th_v' > (momentum levels)                [m/s K]
  !real(r8), pointer, dimension(:,:) :: wp2thvp    ! < w'^2 th_v' > (thermodynamic levels)        [m^2/s^2 K]
  !real(r8), pointer, dimension(:,:) :: rtpthvp    ! < r_t'th_v' > (momentum levels)              [kg/kg K]
  !real(r8), pointer, dimension(:,:) :: thlpthvp   ! < th_l'th_v' > (momentum levels)             [K^2]

  !real(r8), pointer, dimension(:,:) :: upwp     ! east-west momentum flux                      [m^2/s^2]
  !real(r8), pointer, dimension(:,:) :: vpwp     ! north-south momentum flux                    [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: um_pert  ! perturbed meridional wind                    [m/s]
   real(r8), pointer, dimension(:,:) :: vm_pert  ! perturbed zonal wind                         [m/s]
   real(r8), pointer, dimension(:,:) :: upwp_pert! perturbed meridional wind flux               [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vpwp_pert! perturbed zonal wind flux                    [m^2/s^2]

  !real(r8), pointer, dimension(:,:) :: um       ! mean east-west wind                          [m/s]
  !real(r8), pointer, dimension(:,:) :: vm       ! mean north-south wind                        [m/s]
  !real(r8), pointer, dimension(:,:) :: thlm     ! mean temperature                             [K]
  !real(r8), pointer, dimension(:,:) :: rtm      ! mean moisture mixing ratio                   [kg/kg]
  !real(r8), pointer, dimension(:,:) :: rcm        ! CLUBB cloud water mixing ratio               [kg/kg]

   real(r8) :: qclvar(pcols,pverp)              ! cloud water variance                          [kg^2/kg^2]
   real(r8), pointer, dimension(:,:) :: cloud_frac ! CLUBB cloud fraction                       [-]

   real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height                [m]
   real(r8), pointer, dimension(:,:) :: tke      ! turbulent kinetic energy                     [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: relvar   ! relative cloud water variance                [-]
   real(r8), pointer, dimension(:,:) :: accre_enhan ! accretion enhancement factor              [-]
   real(r8), pointer, dimension(:,:) :: cmeliq


  !real(r8), pointer, dimension(:,:) :: naai
   real(r8), pointer, dimension(:,:) :: prer_evap
   real(r8), pointer, dimension(:,:) :: qrl
   real(r8), pointer, dimension(:,:) :: radf_clubb
!PMA
   real(r8)  relvarc(pcols,pver)
  !logical :: lqice(pcnst)  ! det

   integer :: ixorg

   intrinsic :: selected_real_kind, max

  ! for linearize_pbl_winds 
   real(r8) :: sfc_v_diff_tau(pcols) ! Response to tau perturbation, m/s
   real(r8), parameter :: pert_tau = 0.1_r8 ! tau perturbation, Pa
   real(r8), pointer :: wsresp(:)
   real(r8), pointer :: tau_est(:)

!PMA adds gustiness and tpert
   real(r8), pointer :: prec_dp(:)                 ! total precipitation from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow precipitation from ZM convection
   real(r8), pointer :: vmag_gust(:)
   real(r8), pointer :: tpert(:)

#endif

#ifdef CLUBB_SGS

   !===================================
   ! Initialize derived-type variables 
   !===================================
   !  Copy the state to state1 array to use in this routine;
   !  Determine number of columns and which chunk computation is to be performed on

   call physics_state_copy(state,state1)

    ncol = state1%ncol
   lchnk = state1%lchnk

   ! Initialize physics tendency arrays for "all"

   call physics_ptend_init(ptend_all, state1%psetcols, 'macrop_driver_with_clubb')

   !  Determine time step of physics buffer

   itim_old = pbuf_old_tim_idx()

   ! Some tracer indices

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)
   call cnst_get_ind('CLDICE',ixcldice)
   call cnst_get_ind('NUMLIQ',ixnumliq)
   call cnst_get_ind('NUMICE',ixnumice)

   !===========================
   ! Ice Saturation Adjustment
   !===========================
   if (micro_do_icesupersat) then  ! This is .false. in default EAM
      call ice_supersat_adj_tend(state1,pbuf,hdtime,ptend_loc)  ! in, in, in, out
      call physics_ptend_sum(ptend_loc, ptend_all, ncol)        ! Add the ice tendency to the output tendency
      call physics_update(state1, ptend_loc, hdtime)            ! Update state1. ptend_loc is reset to zero by this call
   endif

   !===========================
   ! Turbulent mountain stress
   !===========================
    if ( do_tms) then
       call t_startf('compute_tms')
       call compute_tms( pcols,        pver,      ncol,     &! in
                         state1%u,     state1%v,  state1%t, &! in
                         state1%pmid,  state1%exner,        &! in 
                         state1%zm,    sgh30,               &! in
                         ksrftms,                           &! out. Used below for deriving input to CLUBB
                         tautmsx,   tautmsy,                &! out. Used below for outfld
                         cam_in%landfrac                    )! in
       call t_stopf('compute_tms')

      !Hui Wan note 2023-12: the following lines should be added? I don't see the outfld calls elsewhere.
      !call outfld( 'TAUTMSX'  , tautmsx, pcols, lchnk )
      !call outfld( 'TAUTMSY'  , tautmsy, pcols, lchnk )
    endif

   call pbuf_get_field(pbuf, cloud_frac_idx, cloud_frac, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   !==========================
   ! CLUBB
   !==========================
#include "clubb_tend_eam.inc"

   call physics_ptend_sum(ptend_loc,ptend_all,ncol)   !  Accumulate tendencies for output
   call physics_update(state1,ptend_loc,hdtime)       !  Update the tmp state - state1; ptend_loc is reset to zero after the call

   !===================================
   ! Detrained condensate from deep Cu
   !===================================
   call deepcu_detrainment_tend( state1, ixcldliq, ixcldice, ixnumliq, ixnumice, dlf,                              &! in
                                 clubb_tk1, clubb_tk2, clubb_liq_deep, clubb_liq_sh, clubb_ice_deep, clubb_ice_sh, &! in
                                 ptend_loc, det_s, det_ice                                                         )! out

   call physics_ptend_sum(ptend_loc,ptend_all,ncol)
   call physics_update(state1,ptend_loc,hdtime)

   ! ptend_all now has all accumulated tendencies.  Convert the tendencies for the
   ! dry constituents to dry air basis.
   do ixind = 1, pcnst
      if (lq(ixind) .and. cnst_type(ixind).eq.'dry') then
      ptend_all%q(:ncol,:pver,ixind) = ptend_all%q(:ncol,:pver,ixind)*state1%pdel(:ncol,:pver)/state1%pdeldry(1:ncol,:pver)
      end if
   end do

   !======================================================================================
   ! Some diagnostics from CLUBB
   ! NOTE: a few variables diagnosed here are multiplied by rho which is calculated with
   ! the updated T after detrainment, so moving the respective lines would result in nonBFB
   ! history output (although model integration should still be BFB).
   !======================================================================================
#include "clubb_misc_diag_and_outfld.inc"

   !===================================
   ! Diagnose various cloud fractions
   !===================================
   call cloudfrac_diags( state1, cam_in, cmfmc, itim_old, liqcf_fix, &! in
                         pbuf, cloud_frac,            &! inout
                         alst_o                       )! out 

   !======================
   ! Diagnose PBL height
   !======================
   call t_startf('pbl_depth_diag')
   call pbuf_get_field(pbuf, pblh_idx,    pblh)
   call pblh_diag( state1, cam_in, cloud_frac, dz_g(pver), use_sgv, pblh ) ! 5xin, 1xout
   call outfld('PBLH', pblh, pcols, lchnk)
   call t_stopf('pbl_depth_diag')

   !============
   ! Gustiness
   !============
   call pbuf_get_field(pbuf, prec_dp_idx,   prec_dp)
   call pbuf_get_field(pbuf, snow_dp_idx,   snow_dp)
   call pbuf_get_field(pbuf, vmag_gust_idx, vmag_gust)
   call pbuf_get_field(pbuf, tpert_idx,     tpert)

   call gustiness( use_sgv, state1, cam_in, host_mnts%up2(:,pver), host_mnts%vp2(:,pver), &! in
                   prec_dp, snow_dp, pblh, host_mnts%thlp2,                               &! in
                   vmag_gust, tpert                                                       )! out

   return
#endif

  end subroutine clubb_tend_cam

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

    subroutine clubb_surface (state, ptend, ztodt, cam_in, ustar, obklen)

!-------------------------------------------------------------------------------
! Description: Provide the obukov length and the surface friction velocity
!              for the dry deposition code in routine tphysac.  Since University
!              of Washington Moist Turbulence (UWMT) scheme is not called when
!              CLUBB is turned on the obukov length and ustar are never initialized
!              nor computed (sometimes never updated from NaN).  In addition, surface
!              fluxes are applied to the constituents.
!
! Author: Peter Bogenschutz, August 2011
! Origin: Based heavily on UWMT code (eddy_diff.F90)
! References:
!   None
!-------------------------------------------------------------------------------

    use physics_types,          only: physics_state, physics_ptend, &
                                      physics_ptend_init, &
                                      set_dry_to_wet, set_wet_to_dry
    use physconst,              only: gravit, zvir, latvap
    use ppgrid,                 only: pver, pcols
    use constituents,           only: pcnst, cnst_get_ind, cnst_type
    use co2_cycle,              only: co2_cycle_set_cnst_type
    use camsrfexch,             only: cam_in_t

    implicit none

    ! --------------- !
    ! Input Auguments !
    ! --------------- !

    type(physics_state), intent(inout)  :: state                ! Physics state variables
    type(cam_in_t),      intent(in)     :: cam_in

    real(r8),            intent(in)     :: ztodt                ! 2 delta-t        [ s ]

    ! ---------------- !
    ! Output Auguments !
    ! ---------------- !

    type(physics_ptend), intent(out)    :: ptend                ! Individual parameterization tendencies
    real(r8),            intent(out)    :: obklen(pcols)        ! Obukhov length [ m ]
    real(r8),            intent(out)    :: ustar(pcols)         ! Surface friction velocity [ m/s ]

#ifdef CLUBB_SGS

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer :: i                                                ! indicees
    integer :: ncol                                             ! # of atmospheric columns

    real(r8) :: th(pcols)                                       ! surface potential temperature
    real(r8) :: thv(pcols)                                      ! surface virtual potential temperature
    real(r8) :: kinheat                                         ! kinematic surface heat flux
    real(r8) :: kinwat                                          ! kinematic surface vapor flux
    real(r8) :: kbfs                                            ! kinematic surface buoyancy flux
    real(r8) :: tmp1(pcols)
    real(r8) :: rztodt                                          ! 1./ztodt
    integer  :: m
    integer  :: ixq,ixcldliq !PMA fix for thv
    real(r8) :: rrho                                            ! Inverse air density

    logical  :: lq(pcnst)

    character(len=3), dimension(pcnst) :: cnst_type_loc         ! local override option for constituents cnst_type

#endif
    obklen(pcols) = 0.0_r8
    ustar(pcols)  = 0.0_r8
#ifdef CLUBB_SGS

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! Assume 'wet' mixing ratios in surface diffusion code.
    ! don't convert co2 tracers to wet mixing ratios
    cnst_type_loc(:) = cnst_type(:)
    call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
    call set_dry_to_wet(state, cnst_type_loc)

    call cnst_get_ind('Q',ixq)
    if (use_sgv) then
       call cnst_get_ind('CLDLIQ',ixcldliq)
    endif

    lq(:) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'clubb_srf', lq=lq)

    ncol = state%ncol

    ! Compute the surface friction velocity and obukov length

    do i = 1, ncol
       th(i) = state%t(i,pver)*state%exner(i,pver)         ! diagnose potential temperature
       if (use_sgv) then
         thv(i) = th(i)*(1._r8+zvir*state%q(i,pver,ixq) & ! PMA corrects virtual potential temperature formula
                       - state%q(i,pver,ixcldliq))
       else
         thv(i) = th(i)*(1._r8+zvir*state%q(i,pver,ixq))  ! diagnose virtual potential temperature
       end if
    enddo

    do i = 1, ncol
       call calc_ustar( state%t(i,pver), state%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                        rrho, ustar(i) )
       call calc_obklen( th(i), thv(i), cam_in%cflx(i,1), cam_in%shf(i), rrho, ustar(i), &
                        kinheat, kinwat, kbfs, obklen(i) )
    enddo

    rztodt                 = 1._r8/ztodt
    ptend%q(:ncol,:pver,:) = state%q(:ncol,:pver,:)
    tmp1(:ncol)            = ztodt * gravit * state%rpdel(:ncol,pver)

    do m = 2, pcnst
      ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m) + tmp1(:ncol) * cam_in%cflx(:ncol,m)
    enddo

    ptend%q(:ncol,:pver,:) = (ptend%q(:ncol,:pver,:) - state%q(:ncol,:pver,:)) * rztodt

    ! Convert tendencies of dry constituents to dry basis.
    do m = 1,pcnst
       if (cnst_type(m).eq.'dry') then
          ptend%q(:ncol,:pver,m) = ptend%q(:ncol,:pver,m)*state%pdel(:ncol,:pver)/state%pdeldry(:ncol,:pver)
       endif
    end do
    ! convert wet mmr back to dry before conservation check
    ! avoid converting co2 tracers again
    cnst_type_loc(:) = cnst_type(:)
    call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
    call set_wet_to_dry(state, cnst_type_loc)

    return

#endif

    end subroutine clubb_surface

#ifdef CLUBB_SGS
! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!

real(r8) function diag_ustar( z, bflx, wnd, z0 )

use shr_const_mod, only : shr_const_karman, shr_const_pi, shr_const_g

implicit none

real(r8), parameter      :: am   =  4.8_r8   !   "          "         "
real(r8), parameter      :: bm   = 19.3_r8  !   "          "         "

real(r8), parameter      :: grav = shr_const_g
real(r8), parameter      :: vonk = shr_const_karman
real(r8), parameter      :: pi   = shr_const_pi

real(r8), intent (in)    :: z             ! height where u locates
real(r8), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
real(r8), intent (in)    :: wnd           ! wind speed at z
real(r8), intent (in)    :: z0            ! momentum roughness height


integer :: iterate
real(r8)    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

lnz   = log( z / z0 )
klnz  = vonk/lnz
c1    = pi / 2.0_r8 - 3.0_r8*log( 2.0_r8 )

ustar =  wnd*klnz
if (abs(bflx) > 1.e-6_r8) then
   do iterate=1,4

      if (ustar > 1.e-6_r8) then
         lmo   = -ustar**3 / ( vonk * bflx )
         zeta  = z/lmo
         if (zeta > 0._r8) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
         else
            x     = sqrt( sqrt( 1.0_r8 - bm*zeta ) )
            psi1  = 2._r8*log( 1.0_r8+x ) + log( 1.0_r8+x*x ) - 2._r8*atan( x ) + c1
            ustar = wnd*vonk/(lnz - psi1)
         end if

      endif

   end do
end if


diag_ustar = ustar

return


end function diag_ustar
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS

  subroutine stats_init_clubb( l_stats_in, stats_tsamp_in, stats_tout_in, &
                         nnzp, nnrad_zt,nnrad_zm, delt )
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of CAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with CAM output.

    !-----------------------------------------------------------------------


    use stats_variables, only: &
      stats_zt,      & ! Variables
      ztscr01, &
      ztscr02, &
      ztscr03, &
      ztscr04, &
      ztscr05, &
      ztscr06, &
      ztscr07, &
      ztscr08, &
      ztscr09, &
      ztscr10, &
      ztscr11, &
      ztscr12, &
      ztscr13, &
      ztscr14, &
      ztscr15, &
      ztscr16, &
      ztscr17, &
      ztscr18, &
      ztscr19, &
      ztscr20, &
      ztscr21

    use stats_variables, only: &
      stats_zm,      &
      zmscr01, &
      zmscr02, &
      zmscr03, &
      zmscr04, &
      zmscr05, &
      zmscr06, &
      zmscr07, &
      zmscr08, &
      zmscr09, &
      zmscr10, &
      zmscr11, &
      zmscr12, &
      zmscr13, &
      zmscr14, &
      zmscr15, &
      zmscr16, &
      zmscr17, &
      stats_rad_zt,  &
      stats_rad_zm,  &
      stats_sfc,     &
      l_stats, &
      l_output_rad_files, &
      stats_tsamp,   &
      stats_tout,    &
      l_stats_samp,  &
      l_stats_last, &
      fname_rad_zt, &
      fname_rad_zm, &
      fname_sfc, &
      l_netcdf, &
      l_grads

    use clubb_precision,        only: time_precision, core_rknd   !
    use stats_zm_module,        only: nvarmax_zm, stats_init_zm !
    use stats_zt_module,        only: nvarmax_zt, stats_init_zt !
    use stats_rad_zt_module,    only: nvarmax_rad_zt, stats_init_rad_zt !
    use stats_rad_zm_module,    only: nvarmax_rad_zm, stats_init_rad_zm !
    use stats_sfc_module,       only: nvarmax_sfc, stats_init_sfc !
    use error_code,             only: clubb_at_least_debug_level !
    use constants_clubb,        only: fstderr, var_length !
    use cam_history,            only: addfld, horiz_only
    use namelist_utils,         only: find_group_name
    use units,                  only: getunit, freeunit
    use cam_abortutils,         only: endrun

    implicit none

    ! Input Variables

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    real(kind=time_precision), intent(in) ::  &
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    integer, intent(in) :: nnzp     ! Grid points in the vertical [count]
    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]
    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real(kind=time_precision), intent(in) ::   delt         ! Timestep (dtmain in CLUBB)         [s]


    !  Local Variables

    !  Namelist Variables

    character(len=var_length), dimension(nvarmax_zt)     ::   clubb_vars_zt      ! Variables on the thermodynamic levels
    character(len=var_length), dimension(nvarmax_zm)     ::   clubb_vars_zm      ! Variables on the momentum levels
    character(len=var_length), dimension(nvarmax_rad_zt) ::   clubb_vars_rad_zt  ! Variables on the radiation levels
    character(len=var_length), dimension(nvarmax_rad_zm) ::   clubb_vars_rad_zm  ! Variables on the radiation levels
    character(len=var_length), dimension(nvarmax_sfc)    ::   clubb_vars_sfc     ! Variables at the model surface

    namelist /clubb_stats_nl/ &
      clubb_vars_zt, &
      clubb_vars_zm, &
      clubb_vars_rad_zt, &
      clubb_vars_rad_zm, &
      clubb_vars_sfc

    !  Local Variables

    logical :: l_error

    character(len=200) :: fname, temp1, sub

    integer :: i, ntot, read_status
    integer :: iunit

    !  Initialize
    l_error = .false.

    !  Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in

    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in

    if ( .not. l_stats ) then
       l_stats_samp  = .false.
       l_stats_last  = .false.
       return
    end if

    !  Initialize namelist variables

    clubb_vars_zt     = ''
    clubb_vars_zm     = ''
    clubb_vars_rad_zt = ''
    clubb_vars_rad_zm = ''
    clubb_vars_sfc    = ''

    !  Read variables to compute from the namelist
    if (masterproc) then
       iunit= getunit()
       open(unit=iunit,file="atm_in",status='old')
       call find_group_name(iunit, 'clubb_stats_nl', status=read_status)
       if (read_status == 0) then
          read(unit=iunit, nml=clubb_stats_nl, iostat=read_status)
          if (read_status /= 0) then
             call endrun('stats_init_clubb:  error reading namelist'//errmsg(__FILE__,__LINE__))
          end if
       end if
       close(unit=iunit)
       call freeunit(iunit)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(clubb_vars_zt,      var_length*nvarmax_zt,       mpichar,   0, mpicom)
    call mpibcast(clubb_vars_zm,      var_length*nvarmax_zm,       mpichar,   0, mpicom)
    call mpibcast(clubb_vars_rad_zt,  var_length*nvarmax_rad_zt,   mpichar,   0, mpicom)
    call mpibcast(clubb_vars_rad_zm,  var_length*nvarmax_rad_zm,   mpichar,   0, mpicom)
    call mpibcast(clubb_vars_sfc,     var_length*nvarmax_sfc,      mpichar,   0, mpicom)
#endif

    !  Hardcode these for use in CAM-CLUBB, don't want either
    l_netcdf = .false.
    l_grads  = .false.

    !  Check sampling and output frequencies

    !  The model time step length, delt (which is dtmain), should multiply
    !  evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - floor(stats_tsamp/delt) ) > 1.e-8_r8 ) then
       l_error = .true.  ! This will cause the run to stop.
       write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                        'delt (which is dtmain).  Check the appropriate ',  &
                        'model.in file.'
       write(fstderr,*) 'stats_tsamp = ', stats_tsamp
       write(fstderr,*) 'delt = ', delt
    endif

    !  Initialize zt (mass points)

    i = 1
    do while ( ichar(clubb_vars_zt(i)(1:1)) /= 0 .and. &
               len_trim(clubb_vars_zt(i))   /= 0 .and. &
               i <= nvarmax_zt )
       i = i + 1
    enddo
    ntot = i - 1
    if ( ntot == nvarmax_zt ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_zt than allowed for by nvarmax_zt."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_zt ",  &
                        "in the stats namelist, or change nvarmax_zt."
       write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
       call endrun ("stats_init_clubb:  number of zt statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
    endif

    stats_zt%num_output_fields = ntot
    stats_zt%kk = nnzp

    allocate( stats_zt%z( stats_zt%kk ) )

    allocate( stats_zt%accum_field_values( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%accum_num_samples( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    allocate( stats_zt%l_in_update( 1, 1, stats_zt%kk, stats_zt%num_output_fields ) )
    call stats_zero( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, &
                     stats_zt%accum_num_samples, stats_zt%l_in_update )

    allocate( stats_zt%file%var( stats_zt%num_output_fields ) )
    allocate( stats_zt%file%z( stats_zt%kk ) )

    !  Allocate scratch space

    allocate( ztscr01(stats_zt%kk) )
    allocate( ztscr02(stats_zt%kk) )
    allocate( ztscr03(stats_zt%kk) )
    allocate( ztscr04(stats_zt%kk) )
    allocate( ztscr05(stats_zt%kk) )
    allocate( ztscr06(stats_zt%kk) )
    allocate( ztscr07(stats_zt%kk) )
    allocate( ztscr08(stats_zt%kk) )
    allocate( ztscr09(stats_zt%kk) )
    allocate( ztscr10(stats_zt%kk) )
    allocate( ztscr11(stats_zt%kk) )
    allocate( ztscr12(stats_zt%kk) )
    allocate( ztscr13(stats_zt%kk) )
    allocate( ztscr14(stats_zt%kk) )
    allocate( ztscr15(stats_zt%kk) )
    allocate( ztscr16(stats_zt%kk) )
    allocate( ztscr17(stats_zt%kk) )
    allocate( ztscr18(stats_zt%kk) )
    allocate( ztscr19(stats_zt%kk) )
    allocate( ztscr20(stats_zt%kk) )
    allocate( ztscr21(stats_zt%kk) )

    ztscr01 = 0.0_core_rknd
    ztscr02 = 0.0_core_rknd
    ztscr03 = 0.0_core_rknd
    ztscr04 = 0.0_core_rknd
    ztscr05 = 0.0_core_rknd
    ztscr06 = 0.0_core_rknd
    ztscr07 = 0.0_core_rknd
    ztscr08 = 0.0_core_rknd
    ztscr09 = 0.0_core_rknd
    ztscr10 = 0.0_core_rknd
    ztscr11 = 0.0_core_rknd
    ztscr12 = 0.0_core_rknd
    ztscr13 = 0.0_core_rknd
    ztscr14 = 0.0_core_rknd
    ztscr15 = 0.0_core_rknd
    ztscr16 = 0.0_core_rknd
    ztscr17 = 0.0_core_rknd
    ztscr18 = 0.0_core_rknd
    ztscr19 = 0.0_core_rknd
    ztscr20 = 0.0_core_rknd
    ztscr21 = 0.0_core_rknd

    !  Default initialization for array indices for zt

    call stats_init_zt( clubb_vars_zt, l_error )

    !  Initialize zm (momentum points)

    i = 1
    do while ( ichar(clubb_vars_zm(i)(1:1)) /= 0  .and. &
               len_trim(clubb_vars_zm(i)) /= 0    .and. &
               i <= nvarmax_zm )
       i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zm ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_zm than allowed for by nvarmax_zm."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_zm ",  &
                        "in the stats namelist, or change nvarmax_zm."
       write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
       call endrun ("stats_init_clubb:  number of zm statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
    endif

    stats_zm%num_output_fields = ntot
    stats_zm%kk = nnzp

    allocate( stats_zm%z( stats_zm%kk ) )

    allocate( stats_zm%accum_field_values( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%accum_num_samples( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    allocate( stats_zm%l_in_update( 1, 1, stats_zm%kk, stats_zm%num_output_fields ) )
    call stats_zero( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, &
                     stats_zm%accum_num_samples, stats_zm%l_in_update )

    allocate( stats_zm%file%var( stats_zm%num_output_fields ) )
    allocate( stats_zm%file%z( stats_zm%kk ) )

    !  Allocate scratch space

    allocate( zmscr01(stats_zm%kk) )
    allocate( zmscr02(stats_zm%kk) )
    allocate( zmscr03(stats_zm%kk) )
    allocate( zmscr04(stats_zm%kk) )
    allocate( zmscr05(stats_zm%kk) )
    allocate( zmscr06(stats_zm%kk) )
    allocate( zmscr07(stats_zm%kk) )
    allocate( zmscr08(stats_zm%kk) )
    allocate( zmscr09(stats_zm%kk) )
    allocate( zmscr10(stats_zm%kk) )
    allocate( zmscr11(stats_zm%kk) )
    allocate( zmscr12(stats_zm%kk) )
    allocate( zmscr13(stats_zm%kk) )
    allocate( zmscr14(stats_zm%kk) )
    allocate( zmscr15(stats_zm%kk) )
    allocate( zmscr16(stats_zm%kk) )
    allocate( zmscr17(stats_zm%kk) )

    zmscr01 = 0.0_core_rknd
    zmscr02 = 0.0_core_rknd
    zmscr03 = 0.0_core_rknd
    zmscr04 = 0.0_core_rknd
    zmscr05 = 0.0_core_rknd
    zmscr06 = 0.0_core_rknd
    zmscr07 = 0.0_core_rknd
    zmscr08 = 0.0_core_rknd
    zmscr09 = 0.0_core_rknd
    zmscr10 = 0.0_core_rknd
    zmscr11 = 0.0_core_rknd
    zmscr12 = 0.0_core_rknd
    zmscr13 = 0.0_core_rknd
    zmscr14 = 0.0_core_rknd
    zmscr15 = 0.0_core_rknd
    zmscr16 = 0.0_core_rknd
    zmscr17 = 0.0_core_rknd

    call stats_init_zm( clubb_vars_zm, l_error )

    !  Initialize rad_zt (radiation points)

    if (l_output_rad_files) then

       i = 1
       do while ( ichar(clubb_vars_rad_zt(i)(1:1)) /= 0  .and. &
                  len_trim(clubb_vars_rad_zt(i))   /= 0  .and. &
                  i <= nvarmax_rad_zt )
          i = i + 1
       end do
       ntot = i - 1
       if ( ntot == nvarmax_rad_zt ) then
          write(fstderr,*) "There are more statistical variables listed in ",  &
                           "clubb_vars_rad_zt than allowed for by nvarmax_rad_zt."
          write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zt ",  &
                           "in the stats namelist, or change nvarmax_rad_zt."
          write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
          call endrun ("stats_init_clubb:  number of rad_zt statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
       endif

      stats_rad_zt%num_output_fields = ntot
      stats_rad_zt%kk = nnrad_zt

      allocate( stats_rad_zt%z( stats_rad_zt%kk ) )

      allocate( stats_rad_zt%accum_field_values( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%accum_num_samples( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%l_in_update( 1, 1, stats_rad_zt%kk, stats_rad_zt%num_output_fields ) )

      call stats_zero( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                     stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )

      allocate( stats_rad_zt%file%var( stats_rad_zt%num_output_fields ) )
      allocate( stats_rad_zt%file%z( stats_rad_zt%kk ) )

       fname = trim( fname_rad_zt )

       call stats_init_rad_zt( clubb_vars_rad_zt, l_error )

       !  Initialize rad_zm (radiation points)

       i = 1
       do while ( ichar(clubb_vars_rad_zm(i)(1:1)) /= 0 .and. &
                  len_trim(clubb_vars_rad_zm(i))   /= 0 .and. &
                  i <= nvarmax_rad_zm )
          i = i + 1
       end do
       ntot = i - 1
       if ( ntot == nvarmax_rad_zm ) then
          write(fstderr,*) "There are more statistical variables listed in ",  &
                           "clubb_vars_rad_zm than allowed for by nvarmax_rad_zm."
          write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zm ",  &
                           "in the stats namelist, or change nvarmax_rad_zm."
          write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
          call endrun ("stats_init_clubb:  number of rad_zm statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
       endif

       stats_rad_zm%num_output_fields = ntot
       stats_rad_zm%kk = nnrad_zm

       allocate( stats_rad_zm%z( stats_rad_zm%kk ) )

       allocate( stats_rad_zm%accum_field_values( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%accum_num_samples( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%l_in_update( 1, 1, stats_rad_zm%kk, stats_rad_zm%num_output_fields ) )

       call stats_zero( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                     stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update )

       allocate( stats_rad_zm%file%var( stats_rad_zm%num_output_fields ) )
       allocate( stats_rad_zm%file%z( stats_rad_zm%kk ) )

       fname = trim( fname_rad_zm )

       call stats_init_rad_zm( clubb_vars_rad_zm, l_error )
    end if ! l_output_rad_files


    !  Initialize sfc (surface point)

    i = 1
    do while ( ichar(clubb_vars_sfc(i)(1:1)) /= 0 .and. &
               len_trim(clubb_vars_sfc(i))   /= 0 .and. &
               i <= nvarmax_sfc )
       i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_sfc ) then
       write(fstderr,*) "There are more statistical variables listed in ",  &
                        "clubb_vars_sfc than allowed for by nvarmax_sfc."
       write(fstderr,*) "Check the number of variables listed for clubb_vars_sfc ",  &
                        "in the stats namelist, or change nvarmax_sfc."
       write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
       call endrun ("stats_init_clubb:  number of sfc statistical variables exceeds limit"//errmsg(__FILE__,__LINE__))
    endif

    stats_sfc%num_output_fields = ntot
    stats_sfc%kk = 1

    allocate( stats_sfc%z( stats_sfc%kk ) )

    allocate( stats_sfc%accum_field_values( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%accum_num_samples( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )
    allocate( stats_sfc%l_in_update( 1, 1, stats_sfc%kk, stats_sfc%num_output_fields ) )

    call stats_zero( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, &
                     stats_sfc%accum_num_samples, stats_sfc%l_in_update )

    allocate( stats_sfc%file%var( stats_sfc%num_output_fields ) )
    allocate( stats_sfc%file%z( stats_sfc%kk ) )

    fname = trim( fname_sfc )

    call stats_init_sfc( clubb_vars_sfc, l_error )

    ! Check for errors

    if ( l_error ) then
       call endrun ('stats_init:  errors found'//errmsg(__FILE__,__LINE__))
    endif

!   Now call add fields
    do i = 1, stats_zt%num_output_fields

      temp1 = trim(stats_zt%file%var(i)%name)
      sub   = temp1
      if (len(temp1) .gt. 16) sub = temp1(1:16)

       call addfld(trim(sub),(/ 'ilev' /),&
            'A',trim(stats_zt%file%var(i)%units),trim(stats_zt%file%var(i)%description))
    enddo

    do i = 1, stats_zm%num_output_fields

      temp1 = trim(stats_zm%file%var(i)%name)
      sub   = temp1
      if (len(temp1) .gt. 16) sub = temp1(1:16)

      call addfld(trim(sub),(/ 'ilev' /),&
           'A',trim(stats_zm%file%var(i)%units),trim(stats_zm%file%var(i)%description))
    enddo

    if (l_output_rad_files) then
      do i = 1, stats_rad_zt%num_output_fields
        call addfld(trim(stats_rad_zt%file%var(i)%name),(/ 'foobar' /),&
           'A',trim(stats_rad_zt%file%var(i)%units),trim(stats_rad_zt%file%var(i)%description))
      enddo

      do i = 1, stats_rad_zm%num_output_fields
        call addfld(trim(stats_rad_zm%file%var(i)%name),(/ 'foobar' /),&
           'A',trim(stats_rad_zm%file%var(i)%units),trim(stats_rad_zm%file%var(i)%description))
      enddo
    endif

    do i = 1, stats_sfc%num_output_fields
      call addfld(trim(stats_sfc%file%var(i)%name),horiz_only,&
           'A',trim(stats_sfc%file%var(i)%units),trim(stats_sfc%file%var(i)%description))
    enddo

    return


  end subroutine stats_init_clubb

#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !


    !-----------------------------------------------------------------------
  subroutine stats_end_timestep_clubb(lchnk,thecol,out_zt,out_zm,out_radzt,out_radzm,out_sfc)

    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------

#ifdef CLUBB_SGS

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use stats_variables, only: &
        stats_zt,  & ! Variable(s)
        stats_zm, &
        stats_rad_zt, &
        stats_rad_zm, &
        stats_sfc, &
        l_stats_last, &
        stats_tsamp, &
        stats_tout, &
        l_output_rad_files

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

    use cam_history, only: outfld

    use ppgrid,      only: pcols, pverp

    use cam_abortutils,  only: endrun

    implicit none


#endif

    integer :: lchnk
    integer :: thecol

    real(r8), intent(inout) :: out_zt(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_zm(:,:,:)     ! (pcols,pverp,zt%nn)
    real(r8), intent(inout) :: out_radzt(:,:,:)  ! (pcols,pverp,rad_zt%nn)
    real(r8), intent(inout) :: out_radzm(:,:,:)  ! (pcols,pverp,rad_zm%nn)
    real(r8), intent(inout) :: out_sfc(:,:,:)    ! (pcols,1,sfc%nn)

#ifdef CLUBB_SGS
    ! Local Variables

    integer :: i, k
    logical :: l_error

    !  Check if it is time to write to file

    if ( .not. l_stats_last ) return

    !  Initialize
    l_error = .false.

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zt statistics at each vertical level.
    do i = 1, stats_zt%num_output_fields
      do k = 1, stats_zt%kk

        if ( stats_zt%accum_num_samples(1,1,k,i) /= 0 .and.  &
             stats_zt%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(stats_zt%file%var(i)%name), ' in zt ',  &
                             'at k = ', k,  &
                             '; stats_zt%accum_num_samples(',k,',',i,') = ', stats_zt%accum_num_samples(1,1,k,i)
          endif

        endif

      enddo
    enddo

    !  Look for errors by checking the number of sampling points
    !  for each variable in the zm statistics at each vertical level.
    do i = 1, stats_zm%num_output_fields
      do k = 1, stats_zm%kk

        if ( stats_zm%accum_num_samples(1,1,k,i) /= 0 .and.  &
             stats_zm%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(stats_zm%file%var(i)%name), ' in zm ',  &
                             'at k = ', k,  &
                             '; stats_zm%accum_num_samples(',k,',',i,') = ', stats_zm%accum_num_samples(1,1,k,i)
          endif

        endif

      enddo
    enddo

    if (l_output_rad_files) then
      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zt statistics at each vertical level.
      do i = 1, stats_rad_zt%num_output_fields
        do k = 1, stats_rad_zt%kk

          if ( stats_rad_zt%accum_num_samples(1,1,k,i) /= 0 .and.  &
               stats_rad_zt%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(stats_rad_zt%file%var(i)%name), ' in rad_zt ',  &
                               'at k = ', k,  &
                               '; stats_rad_zt%accum_num_samples(',k,',',i,') = ', stats_rad_zt%accum_num_samples(1,1,k,i)
            endif

          endif

        enddo
      enddo

      !  Look for errors by checking the number of sampling points
      !  for each variable in the rad_zm statistics at each vertical level.
      do i = 1, stats_rad_zm%num_output_fields
        do k = 1, stats_rad_zm%kk

          if ( stats_rad_zm%accum_num_samples(1,1,k,i) /= 0 .and.  &
               stats_rad_zm%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(stats_rad_zm%file%var(i)%name), ' in rad_zm ',  &
                               'at k = ', k,  &
                               '; stats_rad_zm%accum_num_samples(',k,',',i,') = ', stats_rad_zm%accum_num_samples(1,1,k,i)
            endif

          endif

        enddo
      enddo
    end if ! l_output_rad_files

    !  Look for errors by checking the number of sampling points
    !  for each variable in the sfc statistics at each vertical level.
    do i = 1, stats_sfc%num_output_fields
      do k = 1, stats_sfc%kk

        if ( stats_sfc%accum_num_samples(1,1,k,i) /= 0 .and.  &
             stats_sfc%accum_num_samples(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(stats_sfc%file%var(i)%name), ' in sfc ',  &
                             'at k = ', k,  &
                             '; stats_sfc%accum_num_samples(',k,',',i,') = ', stats_sfc%accum_num_samples(1,1,k,i)
          endif

        endif

      enddo
    enddo

    !  Stop the run if errors are found.
    if ( l_error ) then
       write(fstderr,*) 'Possible statistical sampling error'
       write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                        'least 1 in the appropriate model.in file.'
       call endrun ('stats_end_timestep:  error(s) found'//errmsg(__FILE__,__LINE__))
    endif

    !  Compute averages
    call stats_avg( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, stats_zt%accum_num_samples )
    call stats_avg( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, stats_zm%accum_num_samples )
    if (l_output_rad_files) then
      call stats_avg( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                      stats_rad_zt%accum_num_samples )
      call stats_avg( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                      stats_rad_zm%accum_num_samples )
    end if
    call stats_avg( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, stats_sfc%accum_num_samples )

   !  Here we are not outputting the data, rather reading the stats into
   !  arrays which are conformable to CAM output.  Also, the data is "flipped"
   !  in the vertical level to be the same as CAM output.
    do i = 1, stats_zt%num_output_fields
      do k = 1, stats_zt%kk
         out_zt(thecol,k,i) = stats_zt%accum_field_values(1,1,stats_zt%kk-k+1,i)
         if(out_zt(thecol,k,i) .ne. out_zt(thecol,k,i)) out_zt(thecol,k,i) = 0.0_r8
      enddo
    enddo

    do i = 1, stats_zm%num_output_fields
      do k = 1, stats_zt%kk
         out_zm(thecol,k,i) = stats_zm%accum_field_values(1,1,stats_zt%kk-k+1,i)
         if(out_zm(thecol,k,i) .ne. out_zm(thecol,k,i)) out_zm(thecol,k,i) = 0.0_r8
      enddo
    enddo

    if (l_output_rad_files) then
      do i = 1, stats_rad_zt%num_output_fields
        do k = 1, stats_rad_zt%kk
          out_radzt(thecol,k,i) = stats_rad_zt%accum_field_values(1,1,stats_zt%kk-k+1,i)
          if(out_radzt(thecol,k,i) .ne. out_radzt(thecol,k,i)) out_radzt(thecol,k,i) = 0.0_r8
        enddo
      enddo

      do i = 1, stats_rad_zm%num_output_fields
        do k = 1, stats_rad_zm%kk
          out_radzm(thecol,k,i) = stats_rad_zm%accum_field_values(1,1,stats_zt%kk-k+1,i)
          if(out_radzm(thecol,k,i) .ne. out_radzm(thecol,k,i)) out_radzm(thecol,k,i) = 0.0_r8
        enddo
      enddo
    endif

    do i = 1, stats_sfc%num_output_fields
      out_sfc(thecol,1,i) = stats_sfc%accum_field_values(1,1,1,i)
      if(out_sfc(thecol,1,i) .ne. out_sfc(thecol,1,i)) out_sfc(thecol,1,i) = 0.0_r8
    enddo

    !  Reset sample fields
    call stats_zero( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, &
                     stats_zt%accum_num_samples, stats_zt%l_in_update )
    call stats_zero( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, &
                     stats_zm%accum_num_samples, stats_zm%l_in_update )
    if (l_output_rad_files) then
      call stats_zero( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                       stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )
      call stats_zero( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                       stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update )
    end if
    call stats_zero( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, &
                     stats_sfc%accum_num_samples, stats_sfc%l_in_update )

    return

#endif

  end subroutine stats_end_timestep_clubb


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS

    !-----------------------------------------------------------------------
  subroutine stats_zero( kk, nn, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd


    implicit none

    !  Input
    integer, intent(in) :: kk, nn

    !  Output
    real(kind=stat_rknd),    dimension(1,1,kk,nn), intent(out) :: x
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(out) :: n
    logical,                 dimension(1,1,kk,nn), intent(out) :: l_in_update

    !  Zero out arrays

    if ( nn > 0 ) then
       x(:,:,:,:) = 0.0_r8
       n(:,:,:,:) = 0
       l_in_update(:,:,:,:) = .false.
    end if

    return

  end subroutine stats_zero

#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !


#ifdef CLUBB_SGS
    !-----------------------------------------------------------------------
  subroutine stats_avg( kk, nn, x, n )

    !     Description:
    !     Compute the average of stats fields
    !-----------------------------------------------------------------------
    use clubb_precision, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    !  Input
    integer, intent(in) :: nn, kk
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(in) :: n

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(inout)  :: x

    !  Internal

    integer k,m

    !  Compute averages

    do m=1,nn
       do k=1,kk

          if ( n(1,1,k,m) > 0 ) then
             x(1,1,k,m) = x(1,1,k,m) / real( n(1,1,k,m) )
          end if

       end do
    end do

    return

  end subroutine stats_avg


  subroutine determine_clubb_dtime( hdtime, clubb_timestep_in, hdtime_core_rknd, dtime_core_rknd, n_clubb_core_step)

   real(r8),intent(in) :: hdtime              ! host model's cloud macmic timestep
   real(r8),intent(in) :: clubb_timestep_in   ! clubb timestep specified via namelist

   real(core_rknd),intent(out) :: hdtime_core_rknd   ! host model's cloud macmic timestep in core_rknd
   real(core_rknd),intent(out) :: dtime_core_rknd    ! clubb timestep after adjustments (if needed), in core_rknd
   integer,        intent(out) :: n_clubb_core_step  ! number of clubb timesteps per host model's cloud macmic timestep


   !  Determine CLUBB time step and make it sub-step friendly
   !  For now we want CLUBB time step to be 5 min since that is
   !  what has been scientifically validated.  However, there are certain
   !  instances when a 5 min time step will not be possible (based on
   !  host model time step or on macro-micro sub-stepping

    dtime_core_rknd = clubb_timestep_in
   hdtime_core_rknd = real(hdtime, kind = core_rknd)

   !  Now check to see if dtime is greater than the host model
   !    (or sub stepped) time step.  If it is, then simply
   !    set it equal to the host (or sub step) time step.
   !    This section is mostly to deal with small host model
   !    time steps (or small sub-steps)

   if (dtime_core_rknd .gt. hdtime_core_rknd) then
     dtime_core_rknd = hdtime_core_rknd
   endif

   !  Now check to see if CLUBB time step divides evenly into
   !    the host model time step.  If not, force it to divide evenly.
   !    We also want it to be 5 minutes or less.  This section is
   !    mainly for host model time steps that are not evenly divisible
   !    by 5 minutes

   if (mod(hdtime_core_rknd,dtime_core_rknd) .ne. 0) then
     dtime_core_rknd = hdtime_core_rknd/2._core_rknd
     do while (dtime_core_rknd .gt. 300._core_rknd)
       dtime_core_rknd = dtime_core_rknd/2._core_rknd
     end do
   endif

   !  If resulting host model time step and CLUBB time step do not divide evenly
   !    into each other, have model throw a fit.

   if (mod(hdtime_core_rknd,dtime_core_rknd) .ne. 0) then
     call endrun('clubb_tend_cam:  CLUBB time step and HOST time step NOT compatible'//errmsg(__FILE__,__LINE__))
   endif

   !  determine number of timesteps CLUBB core should be advanced,
   !  host time step divided by CLUBB time step
   n_clubb_core_step = max(hdtime_core_rknd/dtime_core_rknd,1._core_rknd)

 end subroutine determine_clubb_dtime

 subroutine column_total_energy( state1, cam_in, i, hdtime, te_b)

   use physics_types, only: physics_state
   use camsrfexch,    only: cam_in_t
   use constituents,  only: cnst_get_ind

   type(physics_state),intent(in) :: state1
   type(cam_in_t),     intent(in) :: cam_in

   integer, intent(in) :: i
   real(r8),intent(in) :: hdtime

   real(r8),intent(out) :: te_b

   integer :: ixq, ixcldliq
   integer :: k

   real(r8) :: invrs_gravit, se_b, ke_b, wv_b, wl_b

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Compute integrals of static energy, kinetic energy, water vapor, and liquid water
   ! for the computation of total energy before CLUBB is called.  This is for an
   ! effort to conserve energy since liquid water potential temperature (which CLUBB
   ! conserves) and static energy (which CAM conserves) are not exactly equal.

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)

   invrs_gravit = 1._r8 / gravit

   se_b = 0._r8  ! initialize vertical integrals
   ke_b = 0._r8
   wv_b = 0._r8
   wl_b = 0._r8

   do k=1,pver ! vertical integral
      ! use s=c_pT+g*z, total energy needs term c_pT but not gz
      se_b = se_b + (state1%s(i,k) - gravit*state1%zm(i,k) - state1%phis(i)) &
                   *  state1%pdel(i,k)*invrs_gravit
      ke_b = ke_b + 0.5_r8*(state1%u(i,k)**2+state1%v(i,k)**2)*state1%pdel(i,k)*invrs_gravit
      wv_b = wv_b + state1%q(i,k,ixq)*state1%pdel(i,k)*invrs_gravit
      wl_b = wl_b + state1%q(i,k,ixcldliq)*state1%pdel(i,k)*invrs_gravit
   enddo

   ! Total energy: sum up the components and also take into account the surface fluxes of heat and moisture
   te_b = se_b + ke_b + (latvap+latice)*wv_b + latice*wl_b
   te_b = te_b +(cam_in%shf(i)+(cam_in%cflx(i,1))*(latvap+latice))*hdtime
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end subroutine column_total_energy

  subroutine column_total_energy_fixer( te_b, thlm, rtm, rcm, exner_clubb, &! in
                                        um, vm, wp2, state1, i,            &! in
                                        clubb_s                            )! out

   use physics_types, only: physics_state
   type(physics_state),intent(in) :: state1
   integer, intent(in)  :: i
   real(r8),intent(in)  :: te_b
   real(r8),intent(in)  :: thlm(:,:), rtm(:,:), rcm(:,:)
   real(r8),intent(in)  :: exner_clubb(:,:), wp2(:,:)
   real(r8),intent(in)  :: um(:,:), vm(:,:) 
   real(r8),intent(out) :: clubb_s(pver)

   real(r8) :: invrs_gravit, se_a, ke_a, wv_a, wl_a, enthalpy, te_a, se_dis

   integer :: k, clubbtop
    

      invrs_gravit = 1._r8 / gravit

      ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
      ! after CLUBB is called.  This is for energy conservation purposes.


      se_a = 0._r8
      ke_a = 0._r8
      wv_a = 0._r8
      wl_a = 0._r8

      do k=1,pver
         enthalpy = cpair*((thlm(i,k)+(latvap/cpair)*rcm(i,k))/exner_clubb(i,k))
         clubb_s(k) = enthalpy + gravit*state1%zm(i,k)+state1%phis(i)
        !se_a = se_a + clubb_s(k)*state1%pdel(i,k)*invrs_gravit
         se_a = se_a + enthalpy * state1%pdel(i,k)*invrs_gravit
         ke_a = ke_a + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)*invrs_gravit
         wv_a = wv_a + (rtm(i,k)-rcm(i,k))*state1%pdel(i,k)*invrs_gravit
         wl_a = wl_a + (rcm(i,k))*state1%pdel(i,k)*invrs_gravit
      enddo

      ! Based on these integrals, compute the total energy before and after CLUBB call
      ! TE as in Williamson2015, E= \int_{whole domain} (K+c_p*T) +
      ! \int_{surface} p_s\phi_s (up to water forms), but we ignore surface term
      ! under assumption that CLUBB does not change surface pressure
      te_a = se_a + ke_a + (latvap+latice)*wv_a + latice*wl_a

      ! Limit the energy fixer to find highest layer where CLUBB is active
      ! Find first level where wp2 is higher than lowest threshold
      clubbtop = 1
      do while (wp2(i,clubbtop) .eq. w_tol_sqd .and. clubbtop .lt. pver-1)
         clubbtop = clubbtop + 1
      enddo

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis = (te_a - te_b)/(state1%pint(i,pverp)-state1%pint(i,clubbtop))

      ! Apply this fixer throughout the column evenly, but only at layers where
      ! CLUBB is active.
      do k=clubbtop,pver
         clubb_s(k) = clubb_s(k) - se_dis*gravit
      enddo

 end subroutine column_total_energy_fixer

#include "clubb_gather_host_fields.inc"
#include "clubb_setup_zgrid_1col.inc"
#include "advance_clubb_core_api_eam.inc"
#include "misc_clubb_intr_subs.inc"

#endif



#include "cloud_frac_diags.inc"


end module clubb_intr
