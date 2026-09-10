!c_doubli--------------------------------------------------------------
! SHOC parameterization
!   SHOC = Simplified Higher Order Closure
!   reference, Bogenschutz and Krueger 2013
!
! PDF-based parameterization for low clouds and turbulence
! email: bogenschutz1@llnl.gov
!--------------------------------------------------------------

! Include bit-for-bit math macros.
#include "bfb_math.inc"

module shoc

  use physics_utils,     only: rtype, rtype8, itype, btype
  use scream_abortutils, only: endscreamrun
  ! EAM host modules used by the MF code (condensation, cold-pool bookkeeping).
  ! The EAMxx CMake build of this file (libshoc: shoc_in_and_out driver, BFB
  ! tests) has no EAM host, so it uses equivalent stand-ins from
  ! eamxx/src/physics/shoc/shoc_eam_host_stubs.F90 instead.
#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_eam_host_stubs, only: qsat, no_ip_hltalt, is_first_step, is_first_restart_step, get_nstep
#else
  use wv_saturation,     only : qsat, no_ip_hltalt
  use time_manager,      only: is_first_step, is_first_restart_step, get_nstep
#endif

! Bit-for-bit math functions.
#ifdef SCREAM_CONFIG_IS_CMAKE
  use physics_share_f2c, only: scream_pow, scream_sqrt, scream_cbrt, scream_gamma, scream_log, &
                               scream_log10, scream_exp, scream_erf
#endif

implicit none
save ! for module variables

public  :: shoc_init, shoc_main, shoc_main_std

logical :: use_cxx = .true.

real(rtype), parameter, public :: largeneg = -99999999.99_rtype
real(rtype), parameter, public :: pi = 3.14159265358979323_rtype

!=========================================================
! Physical constants used in SHOC
!=========================================================

! These are set in initialization and should be set to
!  to the values used in whatever host model SHOC is
!  implemented in
real(rtype) :: ggr   ! gravity [m/s^2]
real(rtype) :: rgas  ! dry air gas constant [J/kg.K]
real(rtype) :: rv    ! water vapor gas constant [J/kg.K]
real(rtype) :: cp    ! specific heat of dry air [J/kg.K]
real(rtype) :: lcond ! latent heat of vaporization [J/kg]
real(rtype) :: lice  ! latent heat of fusion [J/kg]
real(rtype) :: eps   ! rh2o/rair - 1 [-]
real(rtype) :: vk    ! von karmann constant [-]
real(rtype) :: p0    ! Reference pressure, Pa

!=========================================================
! Tunable parameters used in SHOC
!=========================================================

! Set default values, if not overwritten by namelist.
!  All are unitless (unless units are stated)
real(rtype) :: thl2tune = 1.0_rtype ! Temperature variance tuning factor
real(rtype) :: qw2tune = 1.0_rtype ! Moisture variance tuning factor
real(rtype) :: qwthl2tune = 1.0_rtype ! Temperature moisture covariance
real(rtype) :: w2tune = 1.0_rtype ! Vertical velocity variance
real(rtype) :: length_fac = 0.5_rtype ! Length scale factor
real(rtype) :: c_diag_3rd_mom = 7.0_rtype ! w3 factor
real(rtype) :: lambda_low = 0.001_rtype ! lowest value for stability correction
real(rtype) :: lambda_high = 0.04_rtype ! highest value for stability correction
real(rtype) :: lambda_slope = 2.65_rtype ! stability correction slope
real(rtype) :: lambda_thresh = 0.02_rtype ! value to apply stability correction
real(rtype) :: Ckh = 0.1_rtype ! Eddy diffusivity coefficient for heat
real(rtype) :: Ckm = 0.1_rtype ! Eddy diffusivity coefficient for momentum
real(rtype) :: Ckh_s = 0.1_rtype ! Stable PBL diffusivity for heat
real(rtype) :: Ckm_s = 0.1_rtype ! Stable PBL diffusivity for momentum

! MJC: extra tunable parameters
real(rtype) :: l_inf_const = 150.0_rtype  ! [m] Asymptotic value of length scale L
real(rtype) :: tscale_const = 400.0_rtype ! [s] Eddy turnover timescale 
real(rtype) :: Cee_const = 1.0_rtype      ! Turbulent constant of TKE dissipation

! For EDMF:
real(rtype) :: mf_L0   = 50._rtype   ! Default in namelist_defaults_eam.xml: 50 m
real(rtype) :: mf_ent0 = 0.22_rtype  ! Default in namelist_defaults_eam.xml: 0.22
integer     :: mf_nup  = 10          ! Default in namelist_defaults_eam.xml: 10 plumes
real(rtype) :: mf_a = 1._rtype       ! Default in namelist_defaults_eam.xml: 1
real(rtype) :: mf_b = 0.5_rtype      ! Default in namelist_defaults_eam.xml: 0.5
real(rtype) :: mf_c = 0.5_rtype      ! Default in namelist_defaults_eam.xml: 0.5
real(rtype) :: mf_a_wcp = 1._rtype   ! Default in namelist_defaults_eam.xml: 1.0
real(rtype) :: mf_tau_wcp = 4._rtype ! Default in namelist_defaults_eam.xml: 4 h

logical :: do_edmf = .false.
logical :: do_condensation = .false.
logical :: do_precip = .false.
logical :: do_mf_diag = .false.
logical :: do_wthv_mf = .false.
logical :: do_dynamic_L = .false.
logical :: do_entr_tke = .false.
logical :: do_explicit = .false.
logical :: do_integral = .false.
logical :: do_implicit = .false.

! Persistent MF state for shoc_main_std (stands in for EAM's pbuf-carried
! intent(inout) MF fields; see shoc_main_std below)
real(rtype), allocatable, save, private :: mfp_mf_dry_a(:,:)
real(rtype), allocatable, save, private :: mfp_mf_moist_a(:,:)
real(rtype), allocatable, save, private :: mfp_mf_dry_w(:,:)
real(rtype), allocatable, save, private :: mfp_mf_moist_w(:,:)
real(rtype), allocatable, save, private :: mfp_mf_dry_qt(:,:)
real(rtype), allocatable, save, private :: mfp_mf_moist_qt(:,:)
real(rtype), allocatable, save, private :: mfp_mf_dry_thl(:,:)
real(rtype), allocatable, save, private :: mfp_mf_moist_thl(:,:)
real(rtype), allocatable, save, private :: mfp_mf_dry_u(:,:)
real(rtype), allocatable, save, private :: mfp_mf_moist_u(:,:)
real(rtype), allocatable, save, private :: mfp_mf_dry_v(:,:)
real(rtype), allocatable, save, private :: mfp_mf_moist_v(:,:)
real(rtype), allocatable, save, private :: mfp_mf_moist_qc(:,:)
real(rtype), allocatable, save, private :: mfp_mf_thlflx(:,:)
real(rtype), allocatable, save, private :: mfp_mf_qtflx(:,:)
real(rtype), allocatable, save, private :: mfp_mf_thvflx(:,:)
real(rtype), allocatable, save, private :: mfp_mf_qlflx(:,:)
real(rtype), allocatable, save, private :: mfp_mf_ae(:,:)
real(rtype), allocatable, save, private :: mfp_mf_aw(:,:)
real(rtype), allocatable, save, private :: mfp_mf_awthv(:,:)
real(rtype), allocatable, save, private :: mfp_mf_awthl(:,:)
real(rtype), allocatable, save, private :: mfp_mf_awqt(:,:)
real(rtype), allocatable, save, private :: mfp_mf_awql(:,:)
real(rtype), allocatable, save, private :: mfp_mf_awqi(:,:)
real(rtype), allocatable, save, private :: mfp_mf_awu(:,:)
real(rtype), allocatable, save, private :: mfp_mf_awv(:,:)
real(rtype), allocatable, save, private :: mfp_cfl_mf(:,:)
real(rtype), allocatable, save, private :: mfp_mf_auSthl(:,:)
real(rtype), allocatable, save, private :: mfp_mf_auSqt(:,:)
real(rtype), allocatable, save, private :: mfp_mf_auRR(:,:)
real(rtype), allocatable, save, private :: mfp_mf_thvflx_zt(:,:)
real(rtype), allocatable, save, private :: mfp_mf_qlflx_zt(:,:)
real(rtype), allocatable, save, private :: mfp_mf_dry_freq(:)
real(rtype), allocatable, save, private :: mfp_mf_moist_freq(:)
real(rtype), allocatable, save, private :: mfp_plumeheight(:)
real(rtype), allocatable, save, private :: mfp_plume_dry_height(:)
real(rtype), allocatable, save, private :: mfp_mf_w_cp(:)

!=========================================================
! Private module parameters
!=========================================================

! =========
! Below are options to activate certain features in SHOC

! Allow temperature skewness to be independent of moisture
!  variance
logical(btype), parameter :: dothetal_skew = .false.

! 
! ========
! Below define some parameters for SHOC

! Reference temperature [K]
real(rtype), parameter :: basetemp = 300._rtype
! Reference pressure [Pa]
real(rtype), parameter :: basepres = 100000._rtype
! Lower troposphere pressure [Pa]
real(rtype), parameter :: troppres = 80000._rtype
! Minimum surface friction velocity
real(rtype), parameter :: ustar_min = 0.01_rtype
! PBL max depth in pressure units
real(rtype), parameter :: pblmaxp = 4.e4_rtype
! third moment of vertical velocity clipping factor
real(rtype), parameter :: w3clip=1.2_rtype

! ========
! Set upper limits for certain SHOC quantities
! Note that these upper limits are quite high
! and they are rarely reached in a stable simulation

! Mixing length [m]
real(rtype), parameter :: maxlen = 20000.0_rtype
! Minimum Mixing length [m]
real(rtype), parameter :: minlen = 20.0_rtype
! Maximum TKE [m2/s2]
real(rtype), parameter :: maxtke = 50.0_rtype
! Minimum TKE [m2/s2]
real(rtype), parameter :: mintke = 0.0004_rtype

!===================
! const parameter for Diagnosis of PBL depth
real(rtype), parameter :: tinyw = 1.e-36_rtype    ! lower bound for wind magnitude
real(rtype), parameter :: fac   = 100._rtype      ! ustar parameter in height diagnosis
real(rtype), parameter :: ricr  =  0.3_rtype      ! Critical richardson number

! Maximum number of levels in pbl from surface
integer :: npbl

!==============================================================
! Begin SHOC parameterization code!
contains
!==============================================================

subroutine shoc_init( &
         nlev, gravit, rair, rh2o, cpair, &
         zvir, latvap, latice, karman, p0_shoc, &
         pref_mid, nbot_shoc, ntop_shoc, &
         thl2tune_in, qw2tune_in, qwthl2tune_in, &
         w2tune_in, length_fac_in, c_diag_3rd_mom_in, &
         lambda_low_in, lambda_high_in, lambda_slope_in, &
         lambda_thresh_in, Ckh_in, Ckm_in, Ckh_s_in, Ckm_s_in, &
         l_inf_const_in, tscale_const_in, Cee_const_in, &
         mf_L0_in, mf_ent0_in, mf_nup_in, mf_a_in, mf_b_in, mf_c_in, &
         mf_a_wcp_in, mf_tau_wcp_in, &
         do_edmf_in, do_condensation_in, do_precip_in, do_mf_diag_in, &
         do_wthv_mf_in, do_dynamic_L_in, do_entr_tke_in, &
         do_explicit_in, do_integral_in, do_implicit_in )



  implicit none

  ! Purpose:  Initialize constants for SHOC
  ! These should be set to the same values used in
  ! whatever host model SHOC is implemented in

  integer, intent(in)   :: nlev ! number of levels

  real(rtype), intent(in)  :: gravit ! gravity
  real(rtype), intent(in)  :: rair   ! dry air gas constant
  real(rtype), intent(in)  :: rh2o   ! water vapor gas constant
  real(rtype), intent(in)  :: cpair  ! specific heat of dry air
  real(rtype), intent(in)  :: zvir   ! rh2o/rair - 1
  real(rtype), intent(in)  :: latvap ! latent heat of vaporization
  real(rtype), intent(in)  :: latice ! latent heat of fusion
  real(rtype), intent(in)  :: karman ! Von Karman's constant
  real(rtype), intent(in)  :: p0_shoc! Reference pressure, Pa

  real(rtype), intent(in) :: pref_mid(nlev) ! reference pressures at midpoints

  integer, intent(in)   :: nbot_shoc ! Bottom level to which SHOC is applied
  integer, intent(in)   :: ntop_shoc ! Top level to which SHOC is applied

  ! Tunable parameter factors, all unitless
  real(rtype), intent(in), optional :: thl2tune_in ! temperature variance tuning factor
  real(rtype), intent(in), optional :: qw2tune_in ! moisture variance
  real(rtype), intent(in), optional :: qwthl2tune_in ! temperature moisture covariance
  real(rtype), intent(in), optional :: w2tune_in ! vertical velocity variance
  real(rtype), intent(in), optional :: length_fac_in ! length scale
  real(rtype), intent(in), optional :: c_diag_3rd_mom_in ! third moment vertical velocity
  real(rtype), intent(in), optional :: lambda_low_in ! stability correction, lowest value
  real(rtype), intent(in), optional :: lambda_high_in ! stability correciton, highest value
  real(rtype), intent(in), optional :: lambda_slope_in ! stability correction, slope term
  real(rtype), intent(in), optional :: lambda_thresh_in ! value to apply stability correction
  real(rtype), intent(in), optional :: Ckh_in ! eddy diffusivity coefficient for heat
  real(rtype), intent(in), optional :: Ckm_in ! eddy diffusivity coefficient for momentum
  real(rtype), intent(in), optional :: Ckh_s_in ! Stable PBL diffusivity for heat
  real(rtype), intent(in), optional :: Ckm_s_in ! Stable PBL diffusivity for momentum

  ! MJC: extra tunable parameters
  real(rtype), intent(in), optional :: l_inf_const_in  ! Asymptotic value of length scale
  real(rtype), intent(in), optional :: tscale_const_in ! Eddy turnover timescale 
  real(rtype), intent(in), optional :: Cee_const_in    ! Turbulent const of TKE dissipation 

  ! MJC: EDMF parameters and tunable constants
  real(rtype), intent(in), optional :: mf_L0_in 
  real(rtype), intent(in), optional :: mf_ent0_in 
  integer,     intent(in), optional :: mf_nup_in 
  real(rtype), intent(in), optional :: mf_a_in 
  real(rtype), intent(in), optional :: mf_b_in 
  real(rtype), intent(in), optional :: mf_c_in
  real(rtype), intent(in), optional :: mf_a_wcp_in
  real(rtype), intent(in), optional :: mf_tau_wcp_in
   
  logical, intent(in), optional :: do_edmf_in
  logical, intent(in), optional :: do_condensation_in
  logical, intent(in), optional :: do_precip_in
  logical, intent(in), optional :: do_mf_diag_in
  logical, intent(in), optional :: do_wthv_mf_in
  logical, intent(in), optional :: do_dynamic_L_in
  logical, intent(in), optional :: do_entr_tke_in
  logical, intent(in), optional :: do_explicit_in
  logical, intent(in), optional :: do_integral_in
  logical, intent(in), optional :: do_implicit_in

  integer :: k

  ggr = gravit   ! [m/s2]
  rgas = rair    ! [J/kg.K]
  rv = rh2o      ! [J/kg.K]
  cp = cpair     ! [J/kg.K]
  eps = zvir     ! [-]
  lcond = latvap ! [J/kg]
  lice = latice  ! [J/kg]
  vk = karman    ! [-]
  p0 = p0_shoc   ! [Pa]

  ! Tunable parameters, all unitless
  !  override default values if value is present
  if (present(thl2tune_in)) thl2tune = thl2tune_in
  if (present(qw2tune_in)) qw2tune = qw2tune_in
  if (present(qwthl2tune_in)) qwthl2tune = qwthl2tune_in
  if (present(w2tune_in)) w2tune=w2tune_in
  if (present(length_fac_in)) length_fac=length_fac_in
  if (present(qwthl2tune_in)) c_diag_3rd_mom=c_diag_3rd_mom_in
  if (present(lambda_low_in)) lambda_low=lambda_low_in
  if (present(lambda_high_in)) lambda_high=lambda_high_in
  if (present(lambda_slope_in)) lambda_slope=lambda_slope_in
  if (present(lambda_thresh_in)) lambda_thresh=lambda_thresh_in
  if (present(Ckh_in)) Ckh=Ckh_in
  if (present(Ckm_in)) Ckm=Ckm_in
  if (present(Ckh_s_in)) Ckh_s=Ckh_s_in
  if (present(Ckm_s_in)) Ckm_s=Ckm_s_in

  if (present(l_inf_const_in)) l_inf_const = l_inf_const_in
  if (present(tscale_const_in)) tscale_const = tscale_const_in
  if (present(Cee_const_in)) Cee_const = Cee_const_in

  if (present(mf_L0_in)) mf_L0 = mf_L0_in
  if (present(mf_ent0_in)) mf_ent0 = mf_ent0_in
  if (present(mf_nup_in)) mf_nup = mf_nup_in
  if (present(mf_a_in)) mf_a = mf_a_in
  if (present(mf_b_in)) mf_b = mf_b_in
  if (present(mf_c_in)) mf_c = mf_c_in
  if (present(mf_a_wcp_in)) mf_a_wcp = mf_a_wcp_in
  if (present(mf_tau_wcp_in)) mf_tau_wcp = mf_tau_wcp_in
  
  if (present(do_edmf_in)) do_edmf = do_edmf_in
  if (present(do_condensation_in)) do_condensation = do_condensation_in
  if (present(do_precip_in)) do_precip = do_precip_in
  if (present(do_mf_diag_in)) do_mf_diag = do_mf_diag_in
  if (present(do_wthv_mf_in)) do_wthv_mf = do_wthv_mf_in
  if (present(do_dynamic_L_in)) do_dynamic_L = do_dynamic_L_in
  if (present(do_entr_tke_in)) do_entr_tke = do_entr_tke_in
  if (present(do_explicit_in)) do_explicit = do_explicit_in
  if (present(do_integral_in)) do_integral = do_integral_in
  if (present(do_implicit_in)) do_implicit = do_implicit_in

  
   ! Limit pbl height to regions below 400 mb
   ! npbl = max number of levels (from bottom) in pbl

   npbl = 0
   do k=nbot_shoc,ntop_shoc,-1
      if (pref_mid(k) >= pblmaxp) then
         npbl = npbl + 1
      end if
   end do
   npbl = max(npbl,1)

   return

end subroutine shoc_init

!==============================================================
! Main driver for the SHOC scheme
! Host models should call the following routine to call SHOC

subroutine shoc_main ( &
     shcol, nlev, nlevi, dtime, nadv, &   ! Input
     host_dx, host_dy,thv, &              ! Input
     zt_grid,zi_grid,pres,presi,pdel,&    ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, & ! Input
     wtracer_sfc,num_qtracers,w_field, &  ! Input
     inv_exner,phis, &                    ! Input
     host_dse, tke, thetal, qw, &         ! Input/Output
     u_wind, v_wind,qtracers,&            ! Input/Output
     wthv_sec,tkh,tk,&                    ! Input/Output
     shoc_ql,shoc_cldfrac,&               ! Input/Output
     shoc_ql_orig,shoc_cldfrac_orig,mf_ql,& ! Output
     mf_moist_a_zt,mf_moist_qc_zt,&       ! Output (EDMF diagnostic)
     pblh,pblh_wthl,pblh_tke,&            ! Output
     shoc_mix, l_inf, isotropy,&          ! Output (diagnostic)
     a_diss, a_prod_bu, a_prod_sh,&       ! Output (diagnostic)
     w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
     wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
     uw_sec, vw_sec, w3,&                 ! Output (diagnostic)
     wqls_sec, brunt, shoc_ql2,&          ! Output (diagnostic)
     a1_out,C1_out,C2_out,ql1_out,ql2_out,& ! Output
     wthv_sec_ed, &                       ! Output (EDMF diagnostic: ED-only PDF buoyancy flux)
     mf_dry_a,   mf_moist_a,&             ! Output (EDMF diagnostic)
     mf_dry_w,   mf_moist_w,&             ! Output (EDMF diagnostic)
     mf_dry_qt,  mf_moist_qt,&            ! Output (EDMF diagnostic)
     mf_dry_thl, mf_moist_thl,&           ! Output (EDMF diagnostic)
     mf_dry_u,   mf_moist_u,&             ! Output (EDMF diagnostic)
     mf_dry_v,   mf_moist_v,&             ! Output (EDMF diagnostic)
                 mf_moist_qc,&            ! Output (EDMF diagnostic)
     mf_thlflx,  mf_qtflx,&               ! Output (EDMF diagnostic)
     mf_thvflx,  mf_thvflx_zt, &          ! Output (EDMF diagnostic)
     mf_qlflx,   mf_qlflx_zt, &           ! Output (EDMF diagnostic)
     mf_ae, mf_aw, &                      ! Output (EDMF diagnostic)
     mf_awthv, mf_awthl, mf_awqt,&        ! Output (EDMF diagnostic)
     mf_awql, mf_awqi, &                  ! Output (EDMF diagnostic)
     mf_awu, mf_awv, cfl_mf, &            ! Output (EDMF diagnostic)
     mf_auSthl,  mf_auSqt,  mf_auRR, &    ! Output (EDMF micrphysics source terms)
     mf_dry_freq, mf_moist_freq, &        ! Output (EDMF diagnostic)
     plumeheight, plume_dry_height, &     ! Output (EDMF diagnostic)
     ztop, dynamic_L0, ent_ensemble_mean,&! Output (EDMF diagnostic) 
     wstar, qstar, thstar, mf_w_cp, &
     wthl_sec_ed, wthl_sec_mf, &          ! Output (EDMF diagnostic) 
     wqw_sec_ed, wqw_sec_mf &             ! Output - EDMF)              
#ifdef SCREAM_CONFIG_IS_CMAKE
     , elapsed_s &
#endif
     )

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_main_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns in the array
  integer, intent(in) :: shcol
  ! number of levels [-]
  integer, intent(in) :: nlev
  ! number of levels on interface grid [-]
  integer, intent(in) :: nlevi
  ! number of tracers [-]
  integer, intent(in) :: num_qtracers
  ! number of times to loop SHOC
  integer, intent(in) :: nadv

  ! SHOC timestep [s]
  real(rtype), intent(in) :: dtime
  ! grid spacing of host model in x direction [m]
  real(rtype), intent(in) :: host_dx(shcol)
  ! grid spacing of host model in y direction [m]
  real(rtype), intent(in) :: host_dy(shcol)
  ! heights, for thermo grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights, for interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! pressure levels on thermo grid [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! pressure levels on interface grid [Pa]
  real(rtype), intent(in) :: presi(shcol,nlevi)
  ! Differences in pressure levels [Pa]
  real(rtype), intent(in) :: pdel(shcol,nlev)
  ! virtual potential temperature [K]
  real(rtype), intent(in) :: thv(shcol,nlev)
  ! large scale vertical velocity [m/s]
  real(rtype), intent(in) :: w_field(shcol,nlev)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! Surface flux for tracers [varies]
  real(rtype), intent(in) :: wtracer_sfc(shcol,num_qtracers)
  ! Inverse of the exner function [-]
  real(rtype), intent(in) :: inv_exner(shcol,nlev)
  ! Host model surface geopotential height
  real(rtype), intent(in) :: phis(shcol)

! INPUT/OUTPUT VARIABLES
  ! prognostic temp variable of host model
  ! dry static energy [J/kg]
  ! dse = Cp*T + g*z + phis
  real(rtype), intent(inout) :: host_dse(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)
  ! liquid water potential temperature [K]
  real(rtype), intent(inout) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(inout) :: qw(shcol,nlev)
  ! u wind component [m/s]
  real(rtype), intent(inout) :: u_wind(shcol,nlev)
  ! v wind component [m/s]
  real(rtype), intent(inout) :: v_wind(shcol,nlev)
  ! buoyancy flux [K m/s] (pbuf-carried; holds the total ED+MF flux that drives
  ! TKE when do_edmf .and. do_wthv_mf, otherwise the ED-only PDF flux)
  real(rtype), intent(inout) :: wthv_sec(shcol,nlev)
  ! ED-only (SHOC PDF) buoyancy flux, snapshot before wthv_sec is overwritten
  ! with the total; diagnostic only [K m/s]
  real(rtype), intent(out) :: wthv_sec_ed(shcol,nlev)
  ! tracers [varies]
  real(rtype), intent(inout) :: qtracers(shcol,nlev,num_qtracers)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(inout) :: tk(shcol,nlev)
  ! eddy coefficent for heat [m2/s]
  real(rtype), intent(inout) :: tkh(shcol,nlev)
   ! [MJC [10/30/24]: If do_edmf = true, shoc_cldfrac and shoc_ql include MF contributions (see update_cldfrac_ql)
  ! Cloud fraction [-]
  real(rtype), intent(inout) :: shoc_cldfrac(shcol,nlev)
  ! Grid-mean cloud liquid mixing ratio [kg/kg]
  real(rtype), intent(inout) :: shoc_ql(shcol,nlev)
  ! [MJC [10/30/24]:
  ! SHOC's cloud fraction for postprocessing purposes [-]
  real(rtype), intent(out) :: shoc_cldfrac_orig(shcol,nlev) 
  ! SHOC's grid-mean cloud liquid mixing ratio for postprocessing purposes [kg/kg]
  real(rtype), intent(out) :: shoc_ql_orig(shcol,nlev) 
  ! [MJC [5/20/25]: MF's cloud liquid
  real(rtype), intent(out) :: mf_ql(shcol,nlev)
  ! MF moist updraft area fraction interpolated to ZT grid [fraction]
  real(rtype), intent(out) :: mf_moist_a_zt(shcol,nlev)
  ! MF moist updraft condensate interpolated to ZT grid [kg/kg]
  real(rtype), intent(out) :: mf_moist_qc_zt(shcol,nlev)

  ! OUTPUT VARIABLES

  ! planetary boundary layer depth [m]
  real(rtype), intent(out) :: pblh(shcol)
  real(rtype), intent(out) :: pblh_wthl(shcol)  ! MJC: Calculated in edmf_pblh function
  real(rtype), intent(out) :: pblh_tke(shcol)   ! MJC: Calculated in edmf_pblh function
    
  ! cloud liquid mixing ratio variance [kg^2/kg^2]
  real(rtype), intent(out) :: shoc_ql2(shcol, nlev)
  ! MJC [11/03/24]: PDF-related output variables for diagnostic purposes
  real(rtype), intent(out) :: a1_out(shcol, nlev)
  real(rtype), intent(out) :: C1_out(shcol, nlev)
  real(rtype), intent(out) :: C2_out(shcol, nlev)
  real(rtype), intent(out) :: ql1_out(shcol, nlev)
  real(rtype), intent(out) :: ql2_out(shcol, nlev)

  ! also output variables, but part of the SHOC diagnostics
  !  to be output to history file by host model (if desired)

  ! Turbulent length scale [m]
  real(rtype), intent(out) :: shoc_mix(shcol,nlev)
  ! vertical velocity variance [m2/s2]
  real(rtype), intent(out) :: w_sec(shcol,nlev)
  ! temperature variance [K^2]
  real(rtype), intent(out) :: thl_sec(shcol,nlevi)
  ! moisture variance [kg2/kg2]
  real(rtype), intent(out) :: qw_sec(shcol,nlevi)
  ! temp moisture covariance [K kg/kg]
  real(rtype), intent(out) :: qwthl_sec(shcol,nlevi) 
  ! vertical heat flux [K m/s]
  real(rtype), intent(out) :: wthl_sec(shcol,nlevi)
  ! MJC: Extra output variables
  real(rtype), intent(out) :: wthl_sec_ed(shcol,nlevi)
  real(rtype), intent(out) :: wthl_sec_mf(shcol,nlevi)
  
  ! vertical moisture flux [K m/s]
  real(rtype), intent(out) :: wqw_sec(shcol,nlevi)
  ! MJC: Extra output variables
  real(rtype), intent(out) :: wqw_sec_ed(shcol,nlevi)
  real(rtype), intent(out) :: wqw_sec_mf(shcol,nlevi)
  
  ! vertical tke flux [m3/s3]
  real(rtype), intent(out) :: wtke_sec(shcol,nlevi)
  ! vertical zonal momentum flux [m2/s2]
  real(rtype), intent(out) :: uw_sec(shcol,nlevi)
  ! vertical meridional momentum flux [m2/s2]
  real(rtype), intent(out) :: vw_sec(shcol,nlevi)
  ! third moment vertical velocity [m3/s3]
  real(rtype), intent(out) :: w3(shcol,nlevi)
  ! liquid water flux [kg/kg m/s]
  real(rtype), intent(out) :: wqls_sec(shcol,nlev)
  ! brunt vaisala frequency [s-1]
  real(rtype), intent(out) :: brunt(shcol,nlev)
  ! return to isotropic timescale [s]
  real(rtype), intent(out) :: isotropy(shcol,nlev)
  ! MJC: Extra output variables
  real(rtype), intent(out) :: a_diss(shcol,nlev)
  real(rtype), intent(out) :: a_prod_bu(shcol,nlev)
  real(rtype), intent(out) :: a_prod_sh(shcol,nlev)
  real(rtype), intent(out) :: l_inf(shcol)
  
  !! MJC: EDMF VARIABLES
  ! mf_* are diagnostic variables
  real(rtype), intent(inout), dimension(shcol,nlevi) :: &
      mf_dry_a,   mf_moist_a,   & ! dry and moist plume area, respectively [-]
      mf_dry_w,   mf_moist_w,   & ! dry and moist plume mean vertical velocity, respectively [m/s]
      mf_dry_qt,  mf_moist_qt,  & ! dry and moist plume mean qt, respectively [kg/kg]
      mf_dry_thl, mf_moist_thl, & ! dry and moist plume mean thl, respectively [K]
      mf_dry_u,   mf_moist_u,   & ! dry and moist plume mean u, respectively [m/s]
      mf_dry_v,   mf_moist_v,   & ! dry and moist plume mean v, respectively [m/s]
                  mf_moist_qc,  & ! moist plume mean qc [kg/kg]
      mf_thlflx,  mf_qtflx,     & ! total MF turbulent flux of theta_l [K m/s], q_t [(kg/kg) m/s]
      mf_thvflx,  mf_qlflx,     & ! vertical flux of buoyancy (wthv) from mass flux plumes [K m/s]
      mf_ae,    & ! environmental area; by default ae=1 (set constant for now) [-]
      mf_aw,    & ! sum(a_i * w_i) [m/s]
      mf_awthv, & ! sum(a_i * w_i * thv_i) [K m/s]
      mf_awthl, & ! sum(a_i * w_i * thl_i) [K m/s]
      mf_awqt,  & ! sum(a_i * w_i * qt_i)  [(kg/kg) m/s]
      mf_awql,  & ! sum(a_i * w_i * ql_i) [(kg/kg) m/s]
      mf_awqi,  & ! sum(a_i * w_i * qi_i) [(kg/kg) m/s]
      mf_awu,   & ! sum(a_i * w_i * u_i) [m^2/s^2]
      mf_awv,   & ! sum(a_i * w_i * v_i) [m^2/s^2]
      cfl_mf,   &
      mf_auSthl,  mf_auSqt, mf_auRR 

  ! MF buoyancy (thv) flux interpolated to the zt grid [K m/s]
  real(rtype), intent(inout) :: mf_thvflx_zt(shcol,nlev)
  real(rtype), intent(inout) :: mf_qlflx_zt(shcol,nlev)

  ! 2D statistics of plume activation frequency, one for dry and one for moist plumes
  real(rtype), intent(inout) :: mf_dry_freq(shcol), mf_moist_freq(shcol)
  real(rtype), intent(inout) :: plumeheight(shcol), plume_dry_height(shcol)
  real(rtype), intent(inout) :: mf_w_cp(shcol)

  real(rtype), intent(out) :: wstar(shcol), qstar(shcol), thstar(shcol)

  real(rtype), intent(out) :: ztop(shcol), dynamic_L0(shcol)
  real(rtype), intent(out) :: ent_ensemble_mean(shcol,nlev)  
  !! MJC: end of EDMF VARIABLES


#ifdef SCREAM_CONFIG_IS_CMAKE
  real(rtype), optional, intent(out) :: elapsed_s ! duration of main loop in seconds
#endif


  !============================================================================
! LOCAL VARIABLES

  ! time counter
  integer :: t

  ! air density on thermo grid [kg/m3]
  real(rtype) :: rho_zt(shcol,nlev)
  ! SHOC water vapor [kg/kg]
  real(rtype) :: shoc_qv(shcol,nlev)
  ! SHOC temperature [K]
  real(rtype) :: shoc_tabs(shcol,nlev)

  ! Grid difference centereted on thermo grid [m]
  real(rtype) :: dz_zt(shcol,nlev)
  ! Grid difference centereted on interface grid [m]
  real(rtype) :: dz_zi(shcol,nlevi)

  ! Virtual potential temperature recomputed each nadv iteration from the
  ! current substep state (thv-staleness fix, PAESCAL E3SM-fork maint-3.0-scm
  ! commits faf05abdc9 / a80a2140b9). The intent(in) `thv` is computed once per
  ! host call and held constant; inside the nadv loop thetal/qw evolve and
  ! shoc_ql is re-diagnosed, so `thv` goes stale. shoc_thv is rebuilt from
  ! SHOC's own diagnosed temperature (shoc_tabs) and vapor (shoc_qv) with the
  ! same formula as shoc_intr's host-side thv, and is used in place of `thv`
  ! by shoc_length (-> compute_brunt_shoc_length) and by the MF plume model
  ! (integrate_mf and its thv_zi interpolation). With nadv = 1 per host call
  ! (the shoc_in_and_out driver) this is what keeps the plume environment
  ! current over a long integration; in EAM (host refreshes thv each step)
  ! it only matters for nadv > 1.
  real(rtype) :: shoc_thv(shcol,nlev)

  ! Surface friction velocity [m/s]
  real(rtype) :: ustar(shcol)
  ! Monin Obukhov length [m]
  real(rtype) :: obklen(shcol)
  ! Kinematic surface buoyancy flux [m^2/s^3]
  real(rtype) :: kbfs(shcol)

  ! Variables related to energy conservation
  real(rtype) :: se_b(shcol),ke_b(shcol),&
              wv_b(shcol),wl_b(shcol),&
              se_a(shcol),ke_a(shcol),&
              wv_a(shcol),wl_a(shcol)
              
  ! MJC: Extra local variables
  ! air density on interface grid [kg/m3]
  real(rtype) :: rho_zi_mf(shcol,nlevi)
  real(rtype) :: thv_zi(shcol,nlevi)
  real(rtype) :: shoc_ql_zi(shcol,nlevi)
  real(rtype) :: mf_auSthl_zt(shcol,nlev), mf_auSqt_zt(shcol,nlev)

  real(rtype) :: nstep_mf, nstep_help 

#ifdef SCREAM_CONFIG_IS_CMAKE
  integer :: clock_count1, clock_count_rate, clock_count_max, clock_count2, clock_count_diff
#endif

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call shoc_main_f(shcol, nlev, nlevi, dtime, nadv, npbl,& ! Input
                     host_dx, host_dy,thv, &                 ! Input
                     zt_grid,zi_grid,pres,presi,pdel,&       ! Input
                     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &    ! Input
                     wtracer_sfc,num_qtracers,w_field, &     ! Input
                     inv_exner,phis, &                           ! Input
                     host_dse, tke, thetal, qw, &            ! Input/Output
                     u_wind, v_wind,qtracers,&               ! Input/Output
                     wthv_sec,tkh,tk,&                       ! Input/Output
                     shoc_ql,shoc_cldfrac,&                  ! Input/Output
                     pblh,&                                  ! Output
                     shoc_mix, isotropy,&                    ! Output (diagnostic)
                     w_sec, thl_sec, qw_sec, qwthl_sec,&     ! Output (diagnostic)
                     wthl_sec, wqw_sec, wtke_sec,&           ! Output (diagnostic)
                     uw_sec, vw_sec, w3,&                    ! Output (diagnostic)
                     wqls_sec, brunt, shoc_ql2)              ! Output (diagnostic)
     return
  endif
#endif

#ifdef SCREAM_CONFIG_IS_CMAKE
    call system_clock(clock_count1, clock_count_rate, clock_count_max)
#endif

  ! Compute integrals of static energy, kinetic energy, water vapor, and liquid water
  ! for the computation of total energy before SHOC is called.  This is for an
  ! effort to conserve energy since liquid water potential temperature (which SHOC
  ! conserves) and static energy (which E3SM conserves) are not exactly equal.
  call shoc_energy_integrals(&
     shcol,nlev,host_dse,pdel,&             ! Input
     qw,shoc_ql,u_wind,v_wind,&             ! Input
     se_b,ke_b,wv_b,wl_b)                   ! Input/Output

  ! wthv_sec (pbuf-carried buoyancy flux) is consumed by shoc_tke at the top of
  ! the loop and overwritten with the total ED+MF flux after shoc_assumed_pdf
  ! (loop end). On the first sub-iteration shoc_tke uses the value carried over
  ! from the previous timestep, the same convention as tk/tkh.

  do t=1,nadv

    ! Check TKE to make sure values lie within acceptable
    !  bounds after host model performs horizontal advection
    call check_tke(shcol,nlev,&                 ! Input
           tke)                                 ! Input/Output

    ! Define vertical grid arrays needed for
    !   vertical derivatives in SHOC, also
    !   define air density
    call shoc_grid( &
       shcol,nlev,nlevi,&                   ! Input
       zt_grid,zi_grid,pdel,&               ! Input
       dz_zt,dz_zi,rho_zt)                  ! Output

    ! Compute the planetary boundary layer height, which is an
    !   input needed for the length scale calculation.

    ! Update SHOC water vapor, to be used by the next two routines
    call compute_shoc_vapor(&
       shcol,nlev,qw,shoc_ql,&              ! Input
       shoc_qv)                             ! Output

    ! Diagnose absolute temperature
    call compute_shoc_temperature(&
       shcol,nlev,thetal,shoc_ql,inv_exner,& ! Input
       shoc_tabs)                            ! Output

    ! Refresh virtual potential temperature for this substep:
    !   thv = T * inv_exner * (1 + eps*qv - qc)   (as shoc_intr's host-side thv)
    shoc_thv(:,:) = shoc_tabs(:,:) * inv_exner(:,:) &
                    * ( 1.0_rtype + eps * shoc_qv(:,:) - shoc_ql(:,:) )

    call shoc_diag_obklen(&
       shcol,uw_sfc,vw_sfc,&                          ! Input
       wthl_sfc,wqw_sfc,thetal(:shcol,nlev),&         ! Input
       shoc_ql(:shcol,nlev),shoc_qv(:shcol,nlev),&    ! Input
       ustar,kbfs,obklen)                             ! Output

    call pblintd(&
       shcol,nlev,nlevi,&                   ! Input
       zt_grid,zi_grid,thetal,shoc_ql,&     ! Input
       shoc_qv,u_wind,v_wind,&              ! Input
       ustar,obklen,kbfs,shoc_cldfrac,&     ! Input
       pblh)                                ! Output

    ! MJC: PBL height for integrate_mf subroutine
    call edmf_pblh(&
       shcol,nlev,nlevi,&           ! Input
       zt_grid,zi_grid,&            ! Input
       wthl_sec,tke,&       ! Input
       pblh_wthl,pblh_tke)          ! Output  
        
    ! MJC: If using EDMF plumes, diagnose plume properties here
    call linear_interp(zt_grid,zi_grid,shoc_thv,thv_zi,nlev,nlevi,shcol,0._rtype)
    call linear_interp(zt_grid,zi_grid,rho_zt,rho_zi_mf,nlev,nlevi,shcol,0._rtype)
    
    !! MJC: 
    nstep_mf = get_nstep()
    nstep_help = nstep_mf
    if ( nstep_help .eq. 0 ) then
      mf_w_cp(:) = 0._rtype
    endif

    if (do_edmf) then
       call integrate_mf(&
               shcol, nlev, nlevi, dtime,&               ! Input
               rho_zt, rho_zi_mf,&					          	 ! Input
               zt_grid, zi_grid, dz_zt, presi,thv_zi,&   ! Input
               u_wind, v_wind, thetal, shoc_thv, qw,&    ! Input (shoc_thv: refreshed each substep)
               ustar, wthl_sfc, wqw_sfc, shoc_ql, &      ! Input
               pblh, pblh_tke, pblh_wthl, tke, &       	 ! Input
               mf_dry_a,   mf_moist_a, &                 ! Output - updraft diagnostics
               mf_dry_w,   mf_moist_w, &                 ! Output - updraft diagnostics
               mf_dry_qt,  mf_moist_qt, &                ! Output - updraft diagnostics
               mf_dry_thl, mf_moist_thl, &               ! Output - updraft diagnostics
               mf_dry_u,   mf_moist_u,  &                ! Output - updraft diagnostics
               mf_dry_v,   mf_moist_v, &                 ! Output - updraft diagnostics
                           mf_moist_qc, &                ! Output - updraft diagnostics
               mf_ae,      mf_aw, &                      ! Output - for diffusion solver
               mf_awthv, &                               ! Output for total wthv
               mf_awthl,   mf_awqt, &                    ! Output - for diffusion solver
               mf_awql,    mf_awqi, &                    ! Output - for diffusion solver/PDF closure but not coupled yet
               mf_awu,     mf_awv,  &                    ! Output - for diffusion solver/PDF closure but not coupled yet
               mf_w_cp, &
               mf_auSthl,  mf_auSqt, mf_auRR, &          ! Output - source terms from microphysics
               mf_dry_freq,mf_moist_freq, plumeheight, & ! Output
               plume_dry_height, cfl_mf,               & ! Output
               ent_ensemble_mean, ztop, dynamic_L0, &
               wstar,     qstar,   thstar     ) ! Output - 2D statistics of plume activation frequency
    else
       mf_dry_a = 0._rtype
       mf_dry_w = 0._rtype
       mf_dry_qt = 0._rtype
       mf_dry_thl = 0._rtype
       mf_dry_u = 0._rtype
       mf_dry_v = 0._rtype

       mf_moist_a = 0._rtype
       mf_moist_w = 0._rtype
       mf_moist_qt = 0._rtype
       mf_moist_thl = 0._rtype
       mf_moist_u = 0._rtype
       mf_moist_v = 0._rtype
       mf_moist_qc = 0._rtype

       mf_ae = 1._rtype
       mf_aw = 0._rtype
       mf_awthv = 0._rtype
       mf_awthl = 0._rtype
       mf_awqt = 0._rtype
       mf_awql = 0._rtype
       mf_awqi = 0._rtype
       mf_awu = 0._rtype
       mf_awv = 0._rtype
       mf_auSthl = 0._rtype
       mf_auSqt  = 0._rtype
       mf_auRR  = 0._rtype
       mf_dry_freq = 0._rtype
       mf_moist_freq = 0._rtype
    endif
    mf_ae = 1._rtype
     
    ! MF buoyancy (thv) flux diagnostic: kinematic MF flux on the interface grid
    ! (calc_mf_vertflux), then interpolated to the zt grid.
    if (do_edmf) then
       call calc_mf_vertflux(shcol,nlev,nlevi,mf_aw,mf_awthv,thv,thv_zi,mf_thvflx)
       call linear_interp(zi_grid,zt_grid,mf_thvflx,mf_thvflx_zt,nlevi,nlev,shcol,0._rtype)
    else
       mf_thvflx    = 0._rtype
       mf_thvflx_zt = 0._rtype
    endif
    
    !!! Calculate mf_qlflx_zm
    call linear_interp(zt_grid,zi_grid,shoc_ql,shoc_ql_zi,nlev,nlevi,shcol,0._rtype)    
    call calc_mf_vertflux(shcol,nlev,nlevi,mf_aw,mf_awql,shoc_ql,shoc_ql_zi,mf_qlflx)
    call linear_interp(zi_grid,zt_grid,mf_qlflx,mf_qlflx_zt,nlevi,nlev,shcol,0._rtype)

    ! MJC: The buoyancy flux that drives TKE is wthv_sec. When
    ! do_edmf .and. do_wthv_mf, wthv_sec is overwritten after shoc_assumed_pdf
    ! (end of this loop) with the total (ED+MF) eq-7.15 flux; otherwise it stays
    ! the ED-only PDF flux. shoc_tke below consumes the value from the previous
    ! sub-iteration / timestep (pbuf-carried, same convention as tk/tkh).

    ! Update the turbulent length scale
    call shoc_length(&
       	shcol,nlev,nlevi,&                ! Input
       	host_dx,host_dy,&                 ! Input
       	zt_grid,zi_grid,dz_zt,&           ! Input
       	tke,shoc_thv,&                         ! Input
       	brunt,l_inf,shoc_mix)             ! Output

    ! Advance the SGS TKE equation using wthv_sec (total ED+MF buoyancy flux
    ! when do_edmf .and. do_wthv_mf, else the ED-only PDF flux).
    call shoc_tke(&
       shcol,nlev,nlevi,dtime,&                ! Input
       wthv_sec,shoc_mix,&                     ! Input
       dz_zi,dz_zt,pres,shoc_tabs,&            ! Input
       u_wind,v_wind,brunt,&                   ! Input
       zt_grid,zi_grid,pblh,&                  ! Input
       tke,tk,tkh,&                            ! Input/Output
       isotropy, a_diss, a_prod_bu, a_prod_sh) ! Output

    ! Update SHOC prognostic variables here
    !   via implicit diffusion solver
    call update_prognostics_implicit(&         ! Input
       shcol,nlev,nlevi,num_qtracers,&         ! Input
       dtime,dz_zt,dz_zi,rho_zt,&              ! Input
       zt_grid,zi_grid,tk,tkh,&                ! Input
       uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,&        ! Input
       wtracer_sfc,&                           ! Input
       do_mf_diag,mf_ae,mf_aw,mf_awu,mf_awv,&  ! EDMF Input
       mf_awthl,mf_awqt,&                      ! EDMF Input       
       thetal,qw,qtracers,tke,&                ! Input/Output
       u_wind,v_wind)                          ! Input/Output
     
    ! Update thetal and qw by adding the contribution from the 
    ! MF microphysics source terms (eq 5, 20 and 21 of Suselj 2019)
    if (do_edmf) then
      ! Interpolate source terms to midpoints grid
      call linear_interp(zi_grid,zt_grid,mf_auSthl,mf_auSthl_zt,nlevi,nlev,shcol,0._rtype)
      call linear_interp(zi_grid,zt_grid,mf_auSqt,mf_auSqt_zt,nlevi,nlev,shcol,0._rtype)
  
      call update_thermody_mf_source(&
      shcol,nlev,mf_auSthl_zt,mf_auSqt_zt,&  ! Input
      thetal,qw)                             ! Input/Output 
    endif

    ! Diagnose the second order moments
    call diag_second_shoc_moments(&
       shcol,nlev,nlevi, &                    ! Input
       thetal,qw,u_wind,v_wind,tke, &         ! Input
       isotropy,tkh,tk,&                      ! Input
       dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
       wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &   ! Input
       do_mf_diag, mf_ae,  mf_aw, &           ! EDMF Input
       mf_awthl,    mf_awqt, &                ! EDMF Input       
       thl_sec, qw_sec,wthl_sec,wqw_sec,&     ! Output
       qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Output
       w_sec, &                               ! Output
	     mf_thlflx, mf_qtflx, &			        	  ! Output - EDMF)
       wthl_sec_ed,wthl_sec_mf, &             ! Output - EDMF)
       wqw_sec_ed, wqw_sec_mf)                ! Output - EDMF)                  
       
       
    ! Diagnose the third moment of vertical velocity,
    !  needed for the PDF closure
    call diag_third_shoc_moments(&
       shcol,nlev,nlevi,&                   ! Input
       w_sec,thl_sec,&                      ! Input
       wthl_sec,isotropy,brunt,&            ! Input
       thetal,tke,&                         ! Input
       dz_zt,dz_zi,zt_grid,zi_grid,&        ! Input
       w3)                                  ! Output

    ! Call the PDF to close on SGS cloud and turbulence
 	if (do_edmf) then
       call shoc_assumed_pdf(&
       shcol,nlev,nlevi,&                   ! Input
       thetal,qw,w_field,thl_sec,qw_sec,&   ! Input
       !wthl_sec,w_sec,&                    ! Input
       wthl_sec_ed,w_sec,&                  ! Input
       !wqw_sec,qwthl_sec,w3,pres,&         ! Input
       wqw_sec_ed,qwthl_sec,w3,pres,&       ! Input
       mf_qlflx_zt, &					            	! Input ??
       zt_grid,zi_grid,&                    ! Input
       shoc_cldfrac,shoc_ql,&               ! Output
       wqls_sec,wthv_sec,shoc_ql2,&          ! Output
       a1_out,C1_out,C2_out,ql1_out,ql2_out)              ! Output
    else
       call shoc_assumed_pdf(&
       shcol,nlev,nlevi,&                   ! Input
       thetal,qw,w_field,thl_sec,qw_sec,&   ! Input
       wthl_sec,w_sec,&                     ! Input
       wqw_sec,qwthl_sec,w3,pres,&          ! Input
       mf_qlflx_zt, &					            	! Input ??
       zt_grid,zi_grid,&                    ! Input
       shoc_cldfrac,shoc_ql,&               ! Output
       wqls_sec,wthv_sec,shoc_ql2,&          ! Output
       a1_out,C1_out,C2_out,ql1_out,ql2_out)              ! Output
    endif

    ! Snapshot the ED-only PDF buoyancy flux for diagnostics (mirrors
    ! wthl_sec_ed / wqw_sec_ed), before wthv_sec is optionally overwritten below.
    wthv_sec_ed = wthv_sec

    ! MJC: Overwrite wthv_sec with the total (ED+MF) buoyancy flux via eq 7.15.
    ! Placed here, after shoc_assumed_pdf, so all inputs are from the current
    ! sub-iteration: total wthl_sec/wqw_sec, the PDF liquid flux wqls_sec, and
    ! mf_qlflx_zt. wthv_sec is pbuf-carried, so shoc_tke consumes this value at
    ! the top of the next sub-iteration / timestep (one-step lag, same convention
    ! as tk/tkh). Requires do_edmf: without it the MF fluxes are all zero.
    if (do_edmf .and. do_wthv_mf) then
       call buoyancy_total_fluxes(&
           shcol,nlev,nlevi,&                ! Input
           zt_grid,zi_grid,&                 ! Input
           pres,&                            ! Input
           wthl_sec,wqw_sec,&                ! Input
           wqls_sec,mf_qlflx_zt,&            ! Input
           wthv_sec)                         ! Output (overwrites wthv_sec with the total)
    endif

    ! Check TKE to make sure values lie within acceptable
    !  bounds after vertical advection, etc.
    call check_tke(shcol,nlev,tke)
   
    ! MJC [08/02/24]: Update shoc_cldfrac and shoc_ql to include MF contribution
    ! and save the SHOC cloud fraction and SHOC liquid water in shoc_cldfrac_orig and shoc_ql_orig
    if (do_edmf) then
      call update_cldfrac_ql(&
       shcol,nlev,nlevi,&                    ! Input
       zt_grid,zi_grid,&                     ! Input
       mf_moist_a,mf_moist_qc,&              ! Input
       shoc_cldfrac,shoc_ql,&                ! Input/Output
       shoc_cldfrac_orig,shoc_ql_orig,mf_ql,& ! Output
       mf_moist_a_zt,mf_moist_qc_zt)         ! Output
    else
    ! Lets keep a backup of shoc_cldfrac and shoc_ql for postprocessing analysis
      shoc_cldfrac_orig = shoc_cldfrac
      shoc_ql_orig = shoc_ql
      mf_ql = 0._rtype
      mf_moist_a_zt = 0._rtype
      mf_moist_qc_zt = 0._rtype
    endif

  enddo ! end time loop

  ! End SHOC parameterization

  ! Use SHOC outputs to update the host model
  !  temperature
  call update_host_dse(&
     shcol,nlev,thetal,&                   ! Input
     shoc_ql,inv_exner,zt_grid,phis,&          ! Input
     host_dse)                             ! Output

  call shoc_energy_integrals(&             ! Input
     shcol,nlev,host_dse,pdel,&            ! Input
     qw,shoc_ql,u_wind,v_wind,&            ! Input
     se_a,ke_a,wv_a,wl_a)                  ! Output

  call shoc_energy_fixer(&
     shcol,nlev,nlevi,dtime,nadv,&         ! Input
     zt_grid,zi_grid,&                     ! Input
     se_b,ke_b,wv_b,wl_b,&                 ! Input
     se_a,ke_a,wv_a,wl_a,&                 ! Input
     wthl_sfc,wqw_sfc,&                    ! Input
     rho_zt,tke,presi,&                    ! Input
     host_dse)                             ! Input/Output

  ! Remaining code is to diagnose certain quantities
  !  related to PBL.  No answer changing subroutines
  !  should be placed at this point onward.

  ! Update PBLH, as other routines outside of SHOC
  !  may require this variable.

  ! Update SHOC water vapor, to be used by the next two routines
  call compute_shoc_vapor(&
     shcol,nlev,qw,shoc_ql,&              ! Input
     shoc_qv)                             ! Output

  call shoc_diag_obklen(&
     shcol,uw_sfc,vw_sfc,&                          ! Input
     wthl_sfc,wqw_sfc,thetal(:shcol,nlev),&         ! Input
     shoc_ql(:shcol,nlev),shoc_qv(:shcol,nlev),&    ! Input
     ustar,kbfs,obklen)                             ! Output

  call pblintd(&
     shcol,nlev,nlevi,&                   ! Input
     zt_grid,zi_grid,thetal,shoc_ql,&     ! Input
     shoc_qv,u_wind,v_wind,&              ! Input
     ustar,obklen,kbfs,shoc_cldfrac,&     ! Input
     pblh)                                ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
  call system_clock(clock_count2, clock_count_rate, clock_count_max)
  clock_count_diff = clock_count2 - clock_count1
  if (present(elapsed_s)) then
    elapsed_s = real(clock_count_diff) / real(clock_count_rate)
  endif
#endif

  return

end subroutine shoc_main

!==============================================================================
! shoc_main_std: SHOC+MF entry point with standard SHOC's argument list.
!
! The SHOC+MF shoc_main takes ~65 extra, non-optional MF arguments (in EAM they
! are pbuf/history fields owned by shoc_intr). Hosts that only know standard
! SHOC's interface -- the EAMxx C bridge shoc_iso_c::shoc_main_c used by the
! shoc_in_and_out driver and the Fortran BFB tests -- call this wrapper instead.
!  * intent(out) MF/PDF diagnostics of shoc_main go to per-call locals.
!  * intent(inout) MF fields go to module-level arrays (mfp_*) that persist
!    between calls, mimicking pbuf persistence: plume_dry_height is read at the
!    top of the next step to size the plume ensemble (integrate_mf).
!  * Cold-pool velocity scale mfp_mf_w_cp starts at 0 and is only advanced by
!    shoc_main when do_precip = .true. (off for the BOMEX standalone runs).
! All MF switches/constants come from the module variables set by shoc_init or
! shoc_iso_c::shoc_set_param_c; with do_edmf = .false. this is standard SHOC.
!==============================================================================
subroutine shoc_main_std ( &
     shcol, nlev, nlevi, dtime, nadv, &   ! Input
     host_dx, host_dy,thv, &              ! Input
     zt_grid,zi_grid,pres,presi,pdel,&    ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, & ! Input
     wtracer_sfc,num_qtracers,w_field, &  ! Input
     inv_exner,phis, &                    ! Input
     host_dse, tke, thetal, qw, &         ! Input/Output
     u_wind, v_wind,qtracers,&            ! Input/Output
     wthv_sec,tkh,tk,&                    ! Input/Output
     shoc_ql,shoc_cldfrac,&               ! Input/Output
     pblh,&                               ! Output
     shoc_mix, isotropy,&                 ! Output (diagnostic)
     w_sec, thl_sec, qw_sec, qwthl_sec,&  ! Output (diagnostic)
     wthl_sec, wqw_sec, wtke_sec,&        ! Output (diagnostic)
     uw_sec, vw_sec, w3,&                 ! Output (diagnostic)
     wqls_sec, brunt, shoc_ql2 &          ! Output (diagnostic)
#ifdef SCREAM_CONFIG_IS_CMAKE
     , elapsed_s &
#endif
     )

  implicit none

  integer, intent(in) :: shcol, nlev, nlevi, num_qtracers, nadv
  real(rtype), intent(in) :: dtime
  real(rtype), intent(in) :: host_dx(shcol), host_dy(shcol)
  real(rtype), intent(in) :: zt_grid(shcol,nlev), zi_grid(shcol,nlevi)
  real(rtype), intent(in) :: pres(shcol,nlev), presi(shcol,nlevi), pdel(shcol,nlev)
  real(rtype), intent(in) :: thv(shcol,nlev), w_field(shcol,nlev)
  real(rtype), intent(in) :: wthl_sfc(shcol), wqw_sfc(shcol), uw_sfc(shcol), vw_sfc(shcol)
  real(rtype), intent(in) :: wtracer_sfc(shcol,num_qtracers)
  real(rtype), intent(in) :: inv_exner(shcol,nlev), phis(shcol)
  real(rtype), intent(inout) :: host_dse(shcol,nlev), tke(shcol,nlev), thetal(shcol,nlev), qw(shcol,nlev)
  real(rtype), intent(inout) :: u_wind(shcol,nlev), v_wind(shcol,nlev), wthv_sec(shcol,nlev)
  real(rtype), intent(inout) :: qtracers(shcol,nlev,num_qtracers)
  real(rtype), intent(inout) :: tk(shcol,nlev), tkh(shcol,nlev)
  real(rtype), intent(inout) :: shoc_cldfrac(shcol,nlev), shoc_ql(shcol,nlev)
  real(rtype), intent(out) :: pblh(shcol)
  real(rtype), intent(out) :: shoc_ql2(shcol,nlev), shoc_mix(shcol,nlev), w_sec(shcol,nlev)
  real(rtype), intent(out) :: thl_sec(shcol,nlevi), qw_sec(shcol,nlevi), qwthl_sec(shcol,nlevi)
  real(rtype), intent(out) :: wthl_sec(shcol,nlevi), wqw_sec(shcol,nlevi), wtke_sec(shcol,nlevi)
  real(rtype), intent(out) :: uw_sec(shcol,nlevi), vw_sec(shcol,nlevi), w3(shcol,nlevi)
  real(rtype), intent(out) :: wqls_sec(shcol,nlev), brunt(shcol,nlev), isotropy(shcol,nlev)
#ifdef SCREAM_CONFIG_IS_CMAKE
  real(rtype), optional, intent(out) :: elapsed_s ! duration of main loop in seconds
#endif

  ! Per-call intent(out) diagnostics of shoc_main (discarded)
  real(rtype), dimension(shcol,nlev) :: shoc_ql_orig
  real(rtype), dimension(shcol,nlev) :: shoc_cldfrac_orig
  real(rtype), dimension(shcol,nlev) :: mf_ql
  real(rtype), dimension(shcol,nlev) :: mf_moist_a_zt
  real(rtype), dimension(shcol,nlev) :: mf_moist_qc_zt
  real(rtype), dimension(shcol,nlev) :: a_diss
  real(rtype), dimension(shcol,nlev) :: a_prod_bu
  real(rtype), dimension(shcol,nlev) :: a_prod_sh
  real(rtype), dimension(shcol,nlev) :: a1_out
  real(rtype), dimension(shcol,nlev) :: C1_out
  real(rtype), dimension(shcol,nlev) :: C2_out
  real(rtype), dimension(shcol,nlev) :: ql1_out
  real(rtype), dimension(shcol,nlev) :: ql2_out
  real(rtype), dimension(shcol,nlev) :: wthv_sec_ed
  real(rtype), dimension(shcol,nlev) :: ent_ensemble_mean
  real(rtype), dimension(shcol,nlevi) :: wthl_sec_ed
  real(rtype), dimension(shcol,nlevi) :: wthl_sec_mf
  real(rtype), dimension(shcol,nlevi) :: wqw_sec_ed
  real(rtype), dimension(shcol,nlevi) :: wqw_sec_mf
  real(rtype), dimension(shcol) :: pblh_wthl
  real(rtype), dimension(shcol) :: pblh_tke
  real(rtype), dimension(shcol) :: l_inf
  real(rtype), dimension(shcol) :: wstar
  real(rtype), dimension(shcol) :: qstar
  real(rtype), dimension(shcol) :: thstar
  real(rtype), dimension(shcol) :: ztop
  real(rtype), dimension(shcol) :: dynamic_L0

  ! (Re)allocate the persistent MF state on the first call or on a size change
  if (.not. allocated(mfp_mf_dry_a)) then
    call shoc_main_std_alloc(shcol, nlev, nlevi)
  else if (size(mfp_mf_dry_a,1) /= shcol .or. size(mfp_mf_dry_a,2) /= nlevi) then
    call shoc_main_std_alloc(shcol, nlev, nlevi)
  end if

  call shoc_main( &
     shcol=shcol, &
     nlev=nlev, &
     nlevi=nlevi, &
     dtime=dtime, &
     nadv=nadv, &
     host_dx=host_dx, &
     host_dy=host_dy, &
     thv=thv, &
     zt_grid=zt_grid, &
     zi_grid=zi_grid, &
     pres=pres, &
     presi=presi, &
     pdel=pdel, &
     wthl_sfc=wthl_sfc, &
     wqw_sfc=wqw_sfc, &
     uw_sfc=uw_sfc, &
     vw_sfc=vw_sfc, &
     wtracer_sfc=wtracer_sfc, &
     num_qtracers=num_qtracers, &
     w_field=w_field, &
     inv_exner=inv_exner, &
     phis=phis, &
     host_dse=host_dse, &
     tke=tke, &
     thetal=thetal, &
     qw=qw, &
     u_wind=u_wind, &
     v_wind=v_wind, &
     qtracers=qtracers, &
     wthv_sec=wthv_sec, &
     tkh=tkh, &
     tk=tk, &
     shoc_ql=shoc_ql, &
     shoc_cldfrac=shoc_cldfrac, &
     pblh=pblh, &
     shoc_mix=shoc_mix, &
     isotropy=isotropy, &
     w_sec=w_sec, &
     thl_sec=thl_sec, &
     qw_sec=qw_sec, &
     qwthl_sec=qwthl_sec, &
     wthl_sec=wthl_sec, &
     wqw_sec=wqw_sec, &
     wtke_sec=wtke_sec, &
     uw_sec=uw_sec, &
     vw_sec=vw_sec, &
     w3=w3, &
     wqls_sec=wqls_sec, &
     brunt=brunt, &
     shoc_ql2=shoc_ql2, &
     shoc_ql_orig=shoc_ql_orig, &
     shoc_cldfrac_orig=shoc_cldfrac_orig, &
     mf_ql=mf_ql, &
     mf_moist_a_zt=mf_moist_a_zt, &
     mf_moist_qc_zt=mf_moist_qc_zt, &
     a_diss=a_diss, &
     a_prod_bu=a_prod_bu, &
     a_prod_sh=a_prod_sh, &
     a1_out=a1_out, &
     C1_out=C1_out, &
     C2_out=C2_out, &
     ql1_out=ql1_out, &
     ql2_out=ql2_out, &
     wthv_sec_ed=wthv_sec_ed, &
     ent_ensemble_mean=ent_ensemble_mean, &
     wthl_sec_ed=wthl_sec_ed, &
     wthl_sec_mf=wthl_sec_mf, &
     wqw_sec_ed=wqw_sec_ed, &
     wqw_sec_mf=wqw_sec_mf, &
     pblh_wthl=pblh_wthl, &
     pblh_tke=pblh_tke, &
     l_inf=l_inf, &
     wstar=wstar, &
     qstar=qstar, &
     thstar=thstar, &
     ztop=ztop, &
     dynamic_L0=dynamic_L0, &
     mf_dry_a=mfp_mf_dry_a, &
     mf_moist_a=mfp_mf_moist_a, &
     mf_dry_w=mfp_mf_dry_w, &
     mf_moist_w=mfp_mf_moist_w, &
     mf_dry_qt=mfp_mf_dry_qt, &
     mf_moist_qt=mfp_mf_moist_qt, &
     mf_dry_thl=mfp_mf_dry_thl, &
     mf_moist_thl=mfp_mf_moist_thl, &
     mf_dry_u=mfp_mf_dry_u, &
     mf_moist_u=mfp_mf_moist_u, &
     mf_dry_v=mfp_mf_dry_v, &
     mf_moist_v=mfp_mf_moist_v, &
     mf_moist_qc=mfp_mf_moist_qc, &
     mf_thlflx=mfp_mf_thlflx, &
     mf_qtflx=mfp_mf_qtflx, &
     mf_thvflx=mfp_mf_thvflx, &
     mf_qlflx=mfp_mf_qlflx, &
     mf_ae=mfp_mf_ae, &
     mf_aw=mfp_mf_aw, &
     mf_awthv=mfp_mf_awthv, &
     mf_awthl=mfp_mf_awthl, &
     mf_awqt=mfp_mf_awqt, &
     mf_awql=mfp_mf_awql, &
     mf_awqi=mfp_mf_awqi, &
     mf_awu=mfp_mf_awu, &
     mf_awv=mfp_mf_awv, &
     cfl_mf=mfp_cfl_mf, &
     mf_auSthl=mfp_mf_auSthl, &
     mf_auSqt=mfp_mf_auSqt, &
     mf_auRR=mfp_mf_auRR, &
     mf_thvflx_zt=mfp_mf_thvflx_zt, &
     mf_qlflx_zt=mfp_mf_qlflx_zt, &
     mf_dry_freq=mfp_mf_dry_freq, &
     mf_moist_freq=mfp_mf_moist_freq, &
     plumeheight=mfp_plumeheight, &
     plume_dry_height=mfp_plume_dry_height, &
     mf_w_cp=mfp_mf_w_cp &
#ifdef SCREAM_CONFIG_IS_CMAKE
     , elapsed_s=elapsed_s &
#endif
     )

end subroutine shoc_main_std

subroutine shoc_main_std_alloc(shcol, nlev, nlevi)
  ! Allocate and zero the persistent MF state used by shoc_main_std
  implicit none
  integer, intent(in) :: shcol, nlev, nlevi
    if (allocated(mfp_mf_dry_a)) deallocate(mfp_mf_dry_a)
    allocate(mfp_mf_dry_a(shcol,nlevi)); mfp_mf_dry_a = 0._rtype
    if (allocated(mfp_mf_moist_a)) deallocate(mfp_mf_moist_a)
    allocate(mfp_mf_moist_a(shcol,nlevi)); mfp_mf_moist_a = 0._rtype
    if (allocated(mfp_mf_dry_w)) deallocate(mfp_mf_dry_w)
    allocate(mfp_mf_dry_w(shcol,nlevi)); mfp_mf_dry_w = 0._rtype
    if (allocated(mfp_mf_moist_w)) deallocate(mfp_mf_moist_w)
    allocate(mfp_mf_moist_w(shcol,nlevi)); mfp_mf_moist_w = 0._rtype
    if (allocated(mfp_mf_dry_qt)) deallocate(mfp_mf_dry_qt)
    allocate(mfp_mf_dry_qt(shcol,nlevi)); mfp_mf_dry_qt = 0._rtype
    if (allocated(mfp_mf_moist_qt)) deallocate(mfp_mf_moist_qt)
    allocate(mfp_mf_moist_qt(shcol,nlevi)); mfp_mf_moist_qt = 0._rtype
    if (allocated(mfp_mf_dry_thl)) deallocate(mfp_mf_dry_thl)
    allocate(mfp_mf_dry_thl(shcol,nlevi)); mfp_mf_dry_thl = 0._rtype
    if (allocated(mfp_mf_moist_thl)) deallocate(mfp_mf_moist_thl)
    allocate(mfp_mf_moist_thl(shcol,nlevi)); mfp_mf_moist_thl = 0._rtype
    if (allocated(mfp_mf_dry_u)) deallocate(mfp_mf_dry_u)
    allocate(mfp_mf_dry_u(shcol,nlevi)); mfp_mf_dry_u = 0._rtype
    if (allocated(mfp_mf_moist_u)) deallocate(mfp_mf_moist_u)
    allocate(mfp_mf_moist_u(shcol,nlevi)); mfp_mf_moist_u = 0._rtype
    if (allocated(mfp_mf_dry_v)) deallocate(mfp_mf_dry_v)
    allocate(mfp_mf_dry_v(shcol,nlevi)); mfp_mf_dry_v = 0._rtype
    if (allocated(mfp_mf_moist_v)) deallocate(mfp_mf_moist_v)
    allocate(mfp_mf_moist_v(shcol,nlevi)); mfp_mf_moist_v = 0._rtype
    if (allocated(mfp_mf_moist_qc)) deallocate(mfp_mf_moist_qc)
    allocate(mfp_mf_moist_qc(shcol,nlevi)); mfp_mf_moist_qc = 0._rtype
    if (allocated(mfp_mf_thlflx)) deallocate(mfp_mf_thlflx)
    allocate(mfp_mf_thlflx(shcol,nlevi)); mfp_mf_thlflx = 0._rtype
    if (allocated(mfp_mf_qtflx)) deallocate(mfp_mf_qtflx)
    allocate(mfp_mf_qtflx(shcol,nlevi)); mfp_mf_qtflx = 0._rtype
    if (allocated(mfp_mf_thvflx)) deallocate(mfp_mf_thvflx)
    allocate(mfp_mf_thvflx(shcol,nlevi)); mfp_mf_thvflx = 0._rtype
    if (allocated(mfp_mf_qlflx)) deallocate(mfp_mf_qlflx)
    allocate(mfp_mf_qlflx(shcol,nlevi)); mfp_mf_qlflx = 0._rtype
    if (allocated(mfp_mf_ae)) deallocate(mfp_mf_ae)
    allocate(mfp_mf_ae(shcol,nlevi)); mfp_mf_ae = 0._rtype
    if (allocated(mfp_mf_aw)) deallocate(mfp_mf_aw)
    allocate(mfp_mf_aw(shcol,nlevi)); mfp_mf_aw = 0._rtype
    if (allocated(mfp_mf_awthv)) deallocate(mfp_mf_awthv)
    allocate(mfp_mf_awthv(shcol,nlevi)); mfp_mf_awthv = 0._rtype
    if (allocated(mfp_mf_awthl)) deallocate(mfp_mf_awthl)
    allocate(mfp_mf_awthl(shcol,nlevi)); mfp_mf_awthl = 0._rtype
    if (allocated(mfp_mf_awqt)) deallocate(mfp_mf_awqt)
    allocate(mfp_mf_awqt(shcol,nlevi)); mfp_mf_awqt = 0._rtype
    if (allocated(mfp_mf_awql)) deallocate(mfp_mf_awql)
    allocate(mfp_mf_awql(shcol,nlevi)); mfp_mf_awql = 0._rtype
    if (allocated(mfp_mf_awqi)) deallocate(mfp_mf_awqi)
    allocate(mfp_mf_awqi(shcol,nlevi)); mfp_mf_awqi = 0._rtype
    if (allocated(mfp_mf_awu)) deallocate(mfp_mf_awu)
    allocate(mfp_mf_awu(shcol,nlevi)); mfp_mf_awu = 0._rtype
    if (allocated(mfp_mf_awv)) deallocate(mfp_mf_awv)
    allocate(mfp_mf_awv(shcol,nlevi)); mfp_mf_awv = 0._rtype
    if (allocated(mfp_cfl_mf)) deallocate(mfp_cfl_mf)
    allocate(mfp_cfl_mf(shcol,nlevi)); mfp_cfl_mf = 0._rtype
    if (allocated(mfp_mf_auSthl)) deallocate(mfp_mf_auSthl)
    allocate(mfp_mf_auSthl(shcol,nlevi)); mfp_mf_auSthl = 0._rtype
    if (allocated(mfp_mf_auSqt)) deallocate(mfp_mf_auSqt)
    allocate(mfp_mf_auSqt(shcol,nlevi)); mfp_mf_auSqt = 0._rtype
    if (allocated(mfp_mf_auRR)) deallocate(mfp_mf_auRR)
    allocate(mfp_mf_auRR(shcol,nlevi)); mfp_mf_auRR = 0._rtype
    if (allocated(mfp_mf_thvflx_zt)) deallocate(mfp_mf_thvflx_zt)
    allocate(mfp_mf_thvflx_zt(shcol,nlev)); mfp_mf_thvflx_zt = 0._rtype
    if (allocated(mfp_mf_qlflx_zt)) deallocate(mfp_mf_qlflx_zt)
    allocate(mfp_mf_qlflx_zt(shcol,nlev)); mfp_mf_qlflx_zt = 0._rtype
    if (allocated(mfp_mf_dry_freq)) deallocate(mfp_mf_dry_freq)
    allocate(mfp_mf_dry_freq(shcol)); mfp_mf_dry_freq = 0._rtype
    if (allocated(mfp_mf_moist_freq)) deallocate(mfp_mf_moist_freq)
    allocate(mfp_mf_moist_freq(shcol)); mfp_mf_moist_freq = 0._rtype
    if (allocated(mfp_plumeheight)) deallocate(mfp_plumeheight)
    allocate(mfp_plumeheight(shcol)); mfp_plumeheight = 0._rtype
    if (allocated(mfp_plume_dry_height)) deallocate(mfp_plume_dry_height)
    allocate(mfp_plume_dry_height(shcol)); mfp_plume_dry_height = 0._rtype
    if (allocated(mfp_mf_w_cp)) deallocate(mfp_mf_w_cp)
    allocate(mfp_mf_w_cp(shcol)); mfp_mf_w_cp = 0._rtype
end subroutine shoc_main_std_alloc


!==============================================================
! MJC [08/02/24]: Update shoc_ql and shoc_cldfrac to include MF contribution
subroutine update_cldfrac_ql(&
        shcol,nlev,nlevi,&                     ! Input
        zt_grid,zi_grid, &                     ! Input
        mf_moist_a,mf_moist_qc,&               ! Input
        shoc_cldfrac,shoc_ql,&                 ! Input/Output
        shoc_cldfrac_orig,shoc_ql_orig,mf_ql,& ! Output
        mf_moist_a_zt,mf_moist_qc_zt)          ! Output
        
        implicit none
   
        ! INPUT VARIABLES
         ! number of columns [-]
         integer, intent(in) :: shcol
         ! number of mid-point levels [-]
         integer, intent(in) :: nlev
         ! number of interface levels [-]
         integer, intent(in) :: nlevi
         ! mid-point grid heights [m]
         real(rtype), intent(in) :: zt_grid(shcol,nlev)
         ! mid-point grid heights [m]
         real(rtype), intent(in) :: zi_grid(shcol,nlevi)
         ! MF cloud fraction interpolated to zt grid
         real(rtype), intent(in) :: mf_moist_a(shcol,nlevi)
         ! MF cloud liquid interpolated to zt grid
         real(rtype), intent(in) :: mf_moist_qc(shcol,nlevi)
         ! SHOC cloud fraction [-]
         real(rtype), intent(inout) :: shoc_cldfrac(shcol,nlev)
         ! SHOC cloud liquid 
         real(rtype), intent(inout) :: shoc_ql(shcol,nlev)
         ! SHOC+MF cloud fraction [-]
         real(rtype), intent(out) :: shoc_cldfrac_orig(shcol,nlev)
         ! SHOC+MF cloud liquid 
         real(rtype), intent(out) :: shoc_ql_orig(shcol,nlev)
         ! MF cloud liquid
         real(rtype), intent(out) :: mf_ql(shcol,nlev)
         
         ! MF moist updraft area fraction interpolated to ZT grid [fraction]
         real(rtype), intent(out) :: mf_moist_a_zt(shcol,nlev)
         ! MF moist updraft condensate interpolated to ZT grid [kg/kg]
         real(rtype), intent(out) :: mf_moist_qc_zt(shcol,nlev)

         ! Local variables
         integer :: i, k
    
         mf_ql = 0._rtype
         
         call linear_interp(zi_grid,zt_grid,mf_moist_a,mf_moist_a_zt,nlevi,nlev,shcol,largeneg)
         call linear_interp(zi_grid,zt_grid,mf_moist_qc,mf_moist_qc_zt,nlevi,nlev,shcol,largeneg)
         
         ! Lets keep a backup of shoc_cldfrac and shoc_ql for postprocessing analysis
         do k=1,nlev
          do i=1,shcol
            shoc_cldfrac_orig(i,k) = shoc_cldfrac(i,k)
            shoc_ql_orig(i,k) = shoc_ql(i,k)
          enddo
         enddo
         
         ! Here we update shoc_cldfrac and shoc_ql to include the MF contributions
         do k=1,nlev
          do i=1,shcol
            shoc_cldfrac(i,k) = min(1._rtype, shoc_cldfrac(i,k) + mf_moist_a_zt(i,k) )
            shoc_ql(i,k) = shoc_ql(i,k) + mf_moist_a_zt(i,k)*mf_moist_qc_zt(i,k)
            mf_ql(i,k) = mf_moist_a_zt(i,k)*mf_moist_qc_zt(i,k)
          enddo
         enddo
    
end subroutine update_cldfrac_ql
    

!==============================================================
! MJC [08/01/24]: Calculates the buoyancy flux from the total (ED+MF) fluxes,
! using the SAME linearization as shoc_assumed_pdf_compute_buoyancy_flux so the
! only difference from the PDF's ED-only wthv_sec is the use of total fluxes:
!   wthv = wthl + ((1-eps_term)/eps_term)*T0*wqt
!        + ((Lv/cp)*Pi^-1 - T0/eps_term) * (wqls + mf_qlflx)
! where Pi^-1 = (p0/p)^(R/cp) is the inverse Exner function and T0 = basetemp.
! The liquid water flux term follows eq 7.15 of the SHOC tech doc and includes
! both the PDF-derived liquid water flux (wqls_sec) and the MF liquid water
! flux (mf_qlflx_zt).

subroutine buoyancy_total_fluxes(&
        shcol,nlev,nlevi, &       ! Input
        zt_grid,zi_grid, &        ! Input
        pres,&                    ! Input
        wthl_sec,wqw_sec,&        ! Input
        wqls_sec,mf_qlflx_zt,&   ! Input
        wthv_out)                 ! Output

     implicit none

    ! INPUT VARIABLES
     ! number of columns [-]
     integer, intent(in) :: shcol
     ! number of mid-point levels [-]
     integer, intent(in) :: nlev
     ! number of interface levels [-]
     integer, intent(in) :: nlevi
     ! mid-point grid heights [m]
     real(rtype), intent(in) :: zt_grid(shcol,nlev)
     ! interface grid heights [m]
     real(rtype), intent(in) :: zi_grid(shcol,nlevi)
     ! pressure on midpoint grid [Pa]
     real(rtype), intent(in) :: pres(shcol,nlev)
     ! vertical heat flux (total: ED+MF) on interface grid [K m/s]
     real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
     ! vertical moisture flux (total: ED+MF) on interface grid [kg/kg m/s]
     real(rtype), intent(in) :: wqw_sec(shcol,nlevi)
     ! liquid water flux from PDF on midpoint grid [kg/kg m/s]
     real(rtype), intent(in) :: wqls_sec(shcol,nlev)
     ! MF liquid water flux on midpoint grid [kg/kg m/s]
     real(rtype), intent(in) :: mf_qlflx_zt(shcol,nlev)

     ! buoyancy flux on midpoint grid [K m/s]
     real(rtype), intent(out) :: wthv_out(shcol,nlev)

     ! Local variables
     integer :: i, k
     real(rtype) :: wthl_sec_zt(shcol,nlev)
     real(rtype) :: wqw_sec_zt(shcol,nlev)
     real(rtype) :: epsterm
     real(rtype) :: exner_inv  ! inverse Exner function: (p0/p)^(R/cp)
     real(rtype) :: wql_total  ! total liquid water flux (PDF + MF)
     real(rtype) :: liq_coeff  ! coefficient for liquid water flux term

     epsterm = rgas/rv

     ! Interpolate fluxes from interface grid to midpoint grid
     call linear_interp(zi_grid,zt_grid,wthl_sec,wthl_sec_zt,nlevi,nlev,shcol,largeneg)
     call linear_interp(zi_grid,zt_grid,wqw_sec,wqw_sec_zt,nlevi,nlev,shcol,largeneg)

     do k=1,nlev
        do i=1,shcol
           ! Total liquid water flux: PDF environment + MF updrafts
           wql_total = wqls_sec(i,k) + mf_qlflx_zt(i,k)

           ! Inverse Exner function
           exner_inv = bfb_pow(basepres/pres(i,k),(rgas/cp))

           ! Coefficient for liquid water flux term (from eq 7.15)
           liq_coeff = (lcond/cp)*exner_inv - (1._rtype/epsterm)*basetemp

           ! Buoyancy flux: same linearization as the PDF (constant reference
           ! state T0 = basetemp), applied to the total (ED+MF) fluxes
           wthv_out(i,k) = wthl_sec_zt(i,k) &
              + ((1._rtype-epsterm)/epsterm)*basetemp*wqw_sec_zt(i,k) &
              + liq_coeff * wql_total
        enddo
     enddo

   end subroutine buoyancy_total_fluxes

!==============================================================
! Calculates the PBL height: TKE and wthl method (criterium: TKE and wthl vanishes)

subroutine edmf_pblh(&
       shcol,nlev,nlevi, &     ! Input
       zt_grid,zi_grid,  &     ! Input
       wthl_sec,tke,     &     ! Input
       pblh_wthl, pblh_tke)    ! Output

  implicit none

! INPUT VARIABLES
  ! number of columns [-]
  integer, intent(in) :: shcol
  ! number of mid-point levels [-]
  integer, intent(in) :: nlev
  ! number of interface levels [-]
  integer, intent(in) :: nlevi
  ! mid-point grid heights [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! mid-point grid heights [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! vertical heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
  ! TKE [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)

  
  ! Output = PBL height using wthl
  real(rtype), intent(out) :: pblh_wthl(shcol)
  ! Output = PBL height using TKE
  real(rtype), intent(out) :: pblh_tke(shcol)
  
  ! local variables
  integer :: i, k, imax, kstart
  real(rtype) :: tke_threshold
  logical :: check

  ! Scan ceiling for PBL detection: start below this height to skip the stratosphere
  ! and upper troposphere. 20 km covers the tallest tropical convective towers.
  real(rtype), parameter :: pblh_scan_top = 20000._rtype  ! [m]

  ! PBL height: wthl method
  do i=1,shcol
    pblh_wthl(i) = 0._rtype
    check  = .true.
    ! Find first level at or below pblh_scan_top (k=1 is model top, heights decrease with k)
    ! Fallback kstart=1 means scan whole column if threshold somehow not found
    kstart = 1
    do k = 1, nlev
      if (zt_grid(i,k) <= pblh_scan_top) then
        kstart = k
        exit
      end if
    end do
    do k = kstart,nlev  ! from scan ceiling to surface
        if (check .and. abs(abs(wthl_sec(i,k)) - abs(wthl_sec(i,k+1))) > 0.001_rtype) then
         pblh_wthl(i) = zi_grid(i,k)
         check = .false.
       endif
    enddo
  enddo
         
  ! PBL height: TKE method
  do i=1,shcol
    if (tke(i,nlev) .NE. 0._rtype ) then
      tke_threshold = mintke*1.1_rtype
      check  = .true.
      do k = nlev,3,-1    
          if (zt_grid(i,k) < 100._rtype) cycle
          if (check .and. tke(i,k) < tke_threshold .and. tke(i,k-1) < tke_threshold) then
            pblh_tke(i) = zt_grid(i,k)
            check = .false.
          endif
      enddo  
    else
      pblh_tke(i) = 0._rtype
    endif
  enddo
  
end subroutine edmf_pblh

!==============================================================
                         
subroutine update_thermody_mf_source(&
  shcol,nlev,mf_auSthl_zt,mf_auSqt_zt, & ! Input
  thetal,qw)
 
   implicit none
 
 ! INPUT VARIABLES
   ! number of columns [-]
   integer, intent(in) :: shcol
   ! number of mid-point levels [-]
   integer, intent(in) :: nlev
   ! mid-point MF thl microphysics source term [K/s]
   real(rtype), intent(in) :: mf_auSthl_zt(shcol,nlev)
   ! mid-point MF qt microphysics source term [1/s]
   real(rtype), intent(in) :: mf_auSqt_zt(shcol,nlev)

   real(rtype), intent(inout) :: thetal(shcol,nlev)
   real(rtype), intent(inout) :: qw(shcol,nlev)

   ! local variables
   integer :: i
 
   do i=1,shcol
    thetal(i,:) = thetal(i,:) + mf_auSthl_zt(i,:)
    qw(i,:) = qw(i,:) + mf_auSqt_zt(i,:)
   enddo
      
 end subroutine update_thermody_mf_source


!==============================================================
! MJC: This function needs to be tested. For now, this is just a place holder
subroutine integrate_mf(shcol, nz, nzi, dt,                        & ! input
                 rho_zt_in, rho_zi_in,                             & ! input
                 zt_in, zi_in, dz_zt_in, p_in, thv_zi_in,          & ! input - MKW 20200804 removed iex and dz_zi_in
                 u_in,   v_in,   thl_in,   thv_in, qt_in,          & ! input
                 ust,    wthl,   wqt,   qc_in,                     & ! input
                 pblh, pblh_tke, pblh_wthl, tke_in,                & ! input
                 dry_a_out,   moist_a_out,                         & ! output: updraft properties for diagnostics
                 dry_w_out,   moist_w_out,                         & ! output: updraft properties for diagnostics
                 dry_qt_out,  moist_qt_out,                        & ! output: updraft properties for diagnostics
                 dry_thl_out, moist_thl_out,                       & ! output: updraft properties for diagnostics
                 dry_u_out,   moist_u_out,                         & ! output: updraft properties for diagnostics
                 dry_v_out,   moist_v_out,                         & ! output: updraft properties for diagnostics
                              moist_qc_out,                        & ! output: updraft properties for diagnostics
                 ae_out, aw_out,                                   & ! output: variables needed for  diffusion solver
                 awthv_out,                                        & ! output: variable needed for total wthv
                 awthl_out, awqt_out,                              & ! output: variables needed for  diffusion solver
                 awql_out, awqi_out,                               & ! output: variables needed for  diffusion solver
                 awu_out, awv_out,                                 & ! output: variables needed for  diffusion solver
                 mf_w_cp, &
                 auSthl_out, auSqt_out, auRRun_out,                    & ! output: sources terms microphysics
                 freq_dry, freq_moist, plumeheight,                & ! output
                 plume_dry_height, cfl,                            & ! output: frequency of plume activation (2D)
                 ent_ensemble_mean, ztop, dynamic_L0, &
                 wstar,     qstar,   thstar ) 
  ! ================================================================================= !
  ! Original author: Marcin Kurowski, JPL
  ! Modified heavily by Mikael Witte and Maria Chinita, UCLA/JPL for implementation in E3SM
  !
  !
  ! Variables needed for solver:
  ! ae = sum_i (1-a_i)
  ! aw = sum (a_i w_i)
  ! awthl = sum(a_i w_i*thl_i)
  ! awqt  = sum(a_i w_i*qt_i)
  ! awql,awqi,awu,awv similar to above except for different variables - not currently coupled to SHOC diffusion solver
  !
  !
  ! - mass flux variables are computed on edges (i.e. momentum grid):
  !  upa,upw,upqt,... 1:nzi
  !  dry_a,moist_a,dry_w,moist_w, ... 1:nzi
  ! ================================================================================= !
  
     ! ============================================================================== ! 
     ! INPUTS   
     ! physics controls
     integer, intent(in) :: shcol,nz,nzi
     real(rtype), dimension(shcol,nz),  intent(in) :: zt_in,   dz_zt_in, rho_zt_in

     ! MKW TODO: remove zi_in as an argument, was only needed for linear_interp calls that were removed on 2020/09/01
     real(rtype), dimension(shcol,nzi), intent(in) :: zi_in, p_in, thv_zi_in, rho_zi_in
     real(rtype), dimension(shcol,nz),  intent(in) :: u_in,v_in,thl_in,qt_in,qc_in,thv_in, tke_in  ! all on thermodynamic/midpoint levels

     real(rtype), dimension(shcol), intent(in) :: ust,   wthl,   wqt
     real(rtype), dimension(shcol), intent(in) :: pblh, pblh_tke, pblh_wthl
     !time step [s]   
     real(rtype), intent(in) :: dt
     ! ============================================================================== !
     ! OUTPUTS
     ! updraft properties
     real(rtype),dimension(shcol,nzi), intent(out) :: dry_a_out,  &  !moist_a_out,     &
                                                      dry_w_out,    moist_w_out,     &
                                                      dry_qt_out,   moist_qt_out,    &
                                                      dry_thl_out,  moist_thl_out,   &
                                                      dry_u_out,    moist_u_out,     &
                                                      dry_v_out,    moist_v_out,     &
                                                      moist_qc_out
    
     real(rtype),dimension(shcol,nzi), intent(inout) :: moist_a_out  
     
     real(rtype),dimension(shcol), intent(inout) :: mf_w_cp


     ! variables needed for diffusion solver
     real(rtype),dimension(shcol,nzi), intent(out) :: ae_out,       aw_out,          &
                                                      awthv_out,    awthl_out,       &
                                                      awqt_out,     awql_out,        &
                                                      awqi_out,     awu_out,         &
                                                      awv_out,      cfl,             &
                                                      auSthl_out,  auSqt_out,        &
                                                      auRRun_out

     real(rtype),dimension(shcol,nz), intent(out) :: ent_ensemble_mean
      
     ! plume activation frequency
     real(rtype),dimension(shcol),     intent(out) :: freq_dry,     freq_moist
     
     ! plume height from one plume test
     real(rtype), dimension(shcol), intent(out) :: plumeheight, plume_dry_height
     ! ztop and L0
     real(rtype), dimension(shcol), intent(out) :: ztop, dynamic_L0

     real(rtype), dimension(shcol), intent(out) :: wstar,     qstar,   thstar
                                                        
     ! Entrainment variables
     !real(rtype),dimension(nz,mf_nup), intent(out) :: ent
     !integer,    dimension(nz,mf_nup), intent(out) :: enti             
                                      
     !! Flux diagnostics - currently diagnosed elsewhere
     !real(rtype),dimension(shcol,nzi), intent(out) :: thlflx_out, qtflx_out
     
     ! ============================================================================== !
     ! INTERNAL VARIABLES

     ! flipped variables (i.e. here index 1 is at surface)
     real(rtype), dimension(shcol,nz)  :: zt, dz_zt, rho_zt
     real(rtype), dimension(shcol,nzi) :: zi, p, thv_zi, rho_zi
     real(rtype), dimension(shcol,nz)  :: u, v, thl, qt, qc, thv, tke
     ! Auxiliar variable for cloud depth calculation needed for autoconversion time scale
     real(rtype), dimension(shcol,nzi) :: moist_a_aux
     
     ! flipped updraft properties (i.e. index 1 is at surface)
     real(rtype), dimension(shcol,nzi) :: dry_a,     moist_a,      &
                                          dry_w,     moist_w,      &
                                          dry_qt,    moist_qt,     &
                                          dry_thl,   moist_thl,    &
                                          dry_u,     moist_u,      &
                                          dry_v,     moist_v,      &
                                                     moist_qc,     &
                                          dry_th,    moist_th,     &           
                                          ae,        aw,           &
                                          awu,       awv,          &
                                          awthv,     awthl,        &
                                          awqt,      awql,         &
                                          awqv,      awth,         & 
                                          awqi,      awqc,         &
                                          auSthl,    auSqt, auRRun                    
                                                     
                                          
     !real(rtype), dimension(shcol,nzi) :: thlflx, qtflx

     ! sums over all plumes
     !real(rtype), dimension(shcol,nzi) :: moist_th, dry_th, awqv, awth
     

     ! updraft properties
     real(rtype), dimension(nzi,mf_nup) :: upw,      upa,      &
                                           upthl,    upthv,    &
                                                     upth,     &
                                           upqt,     upqc,     &
                                           upql,     upqv,     &
                                           upqi,     ups,      &
                                           upu,      upv,      &
                                           upqs, upLtot, RRun, &
                                           Sac, Sev, Smel,     & ! Sac actually needs to be an ouput so this is wrong, fix it in a bit
                                           Sthl, Sqt

     ! entrainment variables
     real(rtype), dimension(nz,mf_nup) :: entf, ent 
     integer,     dimension(nz,mf_nup) :: enti
     real(rtype), dimension(nz) :: ent_oneplume
     !real(rtype), dimension(mf_nup) :: tau_prec ! autoconversion time-scale following eq 16 of Suselj et al 2019 
     
     
     ! other variables
     integer     :: k,j,i,index_top
     real        :: cfl_zt
     
     
     real(rtype) :: wthv, & !      wstar,     qstar,   thstar,    &
                    sigmaw,   sigmaqt,   sigmath,       z0,    &
                    wmin,        wmax,       wlv,      wtv,    &
                    wp, wp_integral, wstar_aux, thstar_help, mf_a_wcp_local
                    
                    
     real(rtype) :: pbj,            B,       qtn,     thln,    &
                    thvn,         thn,       qcn,      qln,    &
                    qin,           un,        vn,      wn2,    &
                    entexp,   entexpu,      entw,      iexh,   &
                    eturb,  enturb, qsn 
                    
     real(rtype) :: Ltot, entexpqt, cldepth, inv_tau_prec, mf_w_cp_old
     integer :: first_ind_cld, last_ind_cld

     
     real(rtype) :: plumeheight_aux !, plume_dry_height_aux
     real(rtype), dimension(shcol) :: plume_top_height
     ! internal surface cont
     real(rtype) :: dzt(nz) !, dzi(nzi)
  
     ! w parameters
     ! virtual mass coefficients for w-eqn after Suselj etal 2019
     real(rtype),parameter :: wa = 1._rtype,     &
                              wb = 1.5_rtype


     ! parameters defining initial conditions for updrafts
     real(rtype),parameter :: pwmin = 1.5_rtype, &
                              pwmax = 3._rtype

     ! min values to avoid singularities
     real(rtype),parameter :: wstarmin = 1.e-3_rtype,  &
                              pblhmin  = 100._rtype
                              
    ! fixed entrainment rate 
     real(rtype),parameter :: fixent = 1.e-3_rtype
     
     ! threshold value for precipitation formation due to accretion (q0 in eq 15 in Suselj et al 2019)  
     real(rtype),parameter :: acc_thres = 1.25e-3_rtype, &
                              tau_ref   = 15._rtype,     &
                              cldep_min = 1500._rtype,   &
                              cldep_max = 5000._rtype,   &
                              kev       = 2.5e-4_rtype,  & 
                              frain_det = 0._rtype

     logical ::check
     
     integer  :: nstep_mff , nstep_mff_help                           ! current timestep number
    
     integer, parameter :: nmax = 3000
     real(rtype) :: mf_w_cp_help(nmax)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!! BEGIN CODE !!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     nstep_mff = get_nstep()
     nstep_mff_help = nstep_mff
     
     if (nstep_mff_help .eq. 0) then
       mf_w_cp_help = 0._rtype
      !print*,'nstep_mff_help = ',nstep_mff_help
     endif

     ! Flip vertical coordinates and all input variables
     do k=1,nz
       ! thermodynamic grid variables
       zt(:,k)    =  zt_in(:,nz-k+1)
       dz_zt(:,k) =  dz_zt_in(:,nz-k+1)

       u(:,k)     =  u_in(:,nz-k+1)
       v(:,k)     =  v_in(:,nz-k+1)
       
       thl(:,k)   =  thl_in(:,nz-k+1)
       thv(:,k)   =  thv_in(:,nz-k+1)
       
       qt(:,k)    =  qt_in(:,nz-k+1)
       qc(:,k)    =  qc_in(:,nz-k+1)

       tke(:,k)  =  tke_in(:,nz-k+1)

       rho_zt(:,k) =  rho_zt_in(:,nz-k+1)

       ! Interface grid nzi variables
       zi(:,k)    = zi_in(:,nzi-k+1)
       p(:,k)     =  p_in(:,nzi-k+1)
       rho_zi(:,k) =  rho_zi_in(:,nzi-k+1)
       thv_zi(:,k)   =  thv_zi_in(:,nzi-k+1)

       ! auxiliar var for cloud depth calculation
       moist_a_aux(:,k) = moist_a_out(:,nzi-k+1)
     enddo
     zi(:,nzi) = zi_in(:,1)
     p(:,nzi)  =  p_in(:,1)
     rho_zi(:,nzi) = rho_zi_in(:,1)
     thv_zi(:,nzi)   = thv_zi_in(:,1)     

     moist_a_aux(:,nzi) = moist_a_out(:,1)

     ! INITIALIZE OUTPUT VARIABLES
     ! set updraft properties to zero
     dry_a     = 0._rtype
     moist_a   = 0._rtype
     dry_w     = 0._rtype
     moist_w   = 0._rtype
     dry_qt    = 0._rtype
     moist_qt  = 0._rtype
     dry_thl   = 0._rtype
     moist_thl = 0._rtype
     dry_u     = 0._rtype
     moist_u   = 0._rtype
     dry_v     = 0._rtype
     moist_v   = 0._rtype
     moist_qc  = 0._rtype
     moist_th  = 0._rtype
     dry_th    = 0._rtype
     ! outputs - variables needed for solver
     aw        = 0._rtype
     awthv     = 0._rtype
     awthl     = 0._rtype
     awqt      = 0._rtype
     awqc      = 0._rtype
     awqv      = 0._rtype
     awql      = 0._rtype
     awqi      = 0._rtype
     awu       = 0._rtype
     awv       = 0._rtype
     awth      = 0._rtype
     auSthl   = 0._rtype 
     auSqt    = 0._rtype
     auRRun   = 0._rtype

     ! autoconversion timescale
     inv_tau_prec  = 0._rtype
     entexpqt  = 0._rtype
     ! outputs - diagnostics
     !thlflx    = 0._rtype
     !qtflx     = 0._rtype
     
     ! this is the environmental area - by default 1.
     ae = 1._rtype
     ! CFL number
     cfl = 0._rtype

     ! START MAIN COMPUTATION
     ! NOTE: SHOC does not invert the vertical coordinate, which by default is ordered from lowest to highest pressure
     ! (i.e. top of atmosphere to bottom) so surface-based do loops are performed in reverse (i.e. from nz to 1)
     do j=1,shcol
       ! zero out plume properties
       upw   = 0._rtype
       upthl = 0._rtype
       upthv = 0._rtype
       upqt  = 0._rtype
       upa   = 0._rtype
       upu   = 0._rtype
       upv   = 0._rtype
       upqc  = 0._rtype
       ent   = 0._rtype
       upth  = 0._rtype
       upql  = 0._rtype
       upqi  = 0._rtype
       upqs  = 0._rtype
       upLtot = 0._rtype
       upqv  = 0._rtype
       Sac   = 0._rtype
       Sev   = 0._rtype
       Smel  = 0._rtype
       RRun  = 0._rtype 
       Sthl  = 0._rtype
       Sqt   = 0._rtype
       
       ! MJC: 
       mf_w_cp_old = mf_w_cp(j)

       !wthv = wthl(j)+eps*thl(j,1)*wqt(j)
       ! From Kay SCM code (GetSurfaceFluxes.m)
       wthv = wthl(j)*(1._rtype+eps*qt(j,1)) + wqt(j)*eps*thl(j,1)
       
       ! If surface buoyancy is positive then do mass-flux, otherwise not
       if (wthv>0.0) then
         dzt = dz_zt(j,:)
                                 
         if (plume_dry_height(j) > 0._rtype .and. plume_dry_height(j) > pblhmin) then
            pbj = plume_dry_height(j)
         else
            ! For the first time step we don't have a plume_dry_height value
            pbj = pblhmin
         endif
         
         wstar(j)  = max( wstarmin, (ggr/thv(j,1)*wthv*pbj)**(1._rtype/3._rtype) )
         qstar(j)  = wqt(j) / wstar(j)
         thstar(j) = wthv/ wstar(j)
        
         sigmaw  = 0.572_rtype * wstar(j)    
         sigmaqt = 2.89_rtype * abs(qstar(j)) 
         sigmath = 2.89_rtype * abs(thstar(j))
         
         wmin = sigmaw * pwmin
         wmax = sigmaw * pwmax
              
         ! One plume test (height reached by a single plume with a fixed entrainment rate = 1.e-3)  
         ! At the moment, we are not using this as the ztop for dynamic L0    
         ent_oneplume = fixent    
         plumeheight_aux = 0._rtype        
         call oneplume( nz, nzi, zi(j,:), dzt, ent_oneplume, p(j,:), qt(j,:), thl(j,:), thv(j,:),  &
                       thv_zi(j,:), wmax, wmin, sigmaw, sigmaqt, sigmath, wa, wb, &
                       do_condensation, plumeheight_aux)
         
         plumeheight(j) = plumeheight_aux      
         
         ! Compute entrainment coefficient
         ! get dz/L0
         ztop(j) = max(pblh_wthl(j),pbj,pblhmin)
         
         ! MJC: if do_precip = true, change entrainment lenght scale with w_d
         !if (do_precip) then
         ! !dynamic_L0(j) = mf_a*(ztop(j)**mf_b)*max(mf_w_cp_old/wstar(j),1._rtype)
         ! if ( mf_w_cp_help(nstep_mff_help)*mf_a_wcp/wstar(j) > 1.1_rtype ) then 
         !   dynamic_L0(j) = mf_a*(ztop(j)**mf_b)*1.1_rtype !max(10.0_rtype/wstar(j),1._rtype)
         ! else
         !   dynamic_L0(j) = mf_a*(ztop(j)**mf_b)*max(mf_w_cp_help(nstep_mff_help)*mf_a_wcp/wstar(j),1._rtype) 
         ! endif          
         !else
         ! dynamic_L0(j) = mf_a*(ztop(j)**mf_b)
         !endif
         dynamic_L0(j) = mf_a*(ztop(j)**mf_b)

         ! Entf is the Poisson argument 
         do i=1,mf_nup
           do k=1,nz
             if (do_dynamic_L) then
               entf(k,i) = dzt(k) / dynamic_L0(j)
             else
               entf(k,i) = dzt(k) / mf_L0
             end if
           enddo
         enddo
        
         ! Enti is the Poisson random number drawn from the Poisson distribution that represents the
         !number of entrainment events for a given average event frequency equal to L0
         call poisson( nz, mf_nup, entf, enti, u(j,2:5))
          
         ! entrainment: Ent=Ent0/dz*P(dz/L0), Ent0 = mf_ent0 = namelist constant (ent_0 = 0.2)
         do i=1,mf_nup
           do k=1,nz
             ent(k,i) = real( enti(k,i))*mf_ent0/dzt(k)
           enddo
         enddo
         
         do k=1,nz
            ent_ensemble_mean(j,k) = sum(ent(k,:))/mf_nup
         enddo
         
        ! Calculate the autoconversion time scale (eq 16 from Suselj et al 2019)
        ! Here inv_tau_prec = 1/tau_u, i.e. we are calculating 1/tau instead of tau_u
        first_ind_cld = -1
        last_ind_cld = -1
        if (do_precip) then
          ! Calculate cloud depth from previous time step
          do k=2,nzi
            if (moist_a_aux(j,k) > 0._rtype) then
              if (first_ind_cld == -1) then
                first_ind_cld = k
              end if
              last_ind_cld = k
            end if
          enddo
          ! Check if clouds are detected
          if (first_ind_cld == -1) then
            cldepth = 0._rtype
          else
            cldepth = zt(j,last_ind_cld) - zt(j,first_ind_cld-1)
          endif
  
          ! Calculate inv_tau_prec
          if (cldepth < cldep_min ) then
            inv_tau_prec =  0._rtype
          else if (cldepth .ge. cldep_min .and. cldepth .le. cldep_max ) then
            inv_tau_prec = ( (cldepth - cldep_min )/(cldep_max - cldep_min) ) / tau_ref
            !print*,'inv_tau_prec = ',inv_tau_prec,', cldepth = ',cldepth
          else if (cldepth > cldep_max ) then
            inv_tau_prec = 1._rtype / tau_ref 
          end if  
 
          ! If cldepth depends on updraft:
          !do i=1,mf_nup
          !  ! If cloud depth (cldepth) < cldep_min (1500 m), tau_prec = 0
          !  if (cldepth(i) < cldep_min ) then
          !    tau_prec(i) = 0._rtype
          !  else if (cldepth(i) .ge. cldep_min .and. cldepth(i) .le. cldep_max ) then
          !    tau_prec(i) = ( (cldepth(i) - cldep_min )/(cldep_max - cldep_min) ) / tau_ref
          !    print*,'tau_prec = ',tau_prec,', cldepth = ',cldepth(i)
          !  else if (cldepth(i) > cldep_max ) then
          !    tau_prec(i) = 1._rtype / tau_ref 
          !  end if
          !enddo
        end if !do_precip

        ! Loop through number of updrafts mf_nup
         do i=1,mf_nup
           ! wlv = w_min_i
           wlv = wmin + (wmax-wmin) / (real(mf_nup)) * (real(i)-1._rtype)
           ! wtv = w_max_i
           wtv = wmin + (wmax-wmin) / (real(mf_nup)) * real(i)

           ! Surface vertical velocity of updraft i: w_i
           upw(1,i) = 0.5_rtype * (wlv+wtv)
           ! Surface area of updraft i: a_i
           upa(1,i) = 0.5_rtype * erf( wtv/(sqrt(2._rtype)*sigmaw) ) &
                      - 0.5_rtype * erf( wlv/(sqrt(2._rtype)*sigmaw) )

           upu(1, i) = u(j,1)
           upv(1, i) = v(j,1)

           upqt(1,i)  = qt(j,1)  + 0.32_rtype * upw(1,i) * sigmaqt/sigmaw
           upthv(1,i) = thv(j,1) + 0.58_rtype * upw(1,i) * sigmath/sigmaw
           upthl(1,i) = upthv(1,i) / (1._rtype+eps*upqt(1,i))
           
           upqv(1,i)  = upqt(1,i)
           
           if (do_condensation) then
                iexh = (1.e5_rtype / p(j,1))**(rgas/cp)
                call condensation_mf(upqt(1,i), upthl(1,i), p(j,1), iexh, &
                                     thvn, qcn, thn, qln, qin, qsn, Ltot)
                upthv(1,i) = thvn
                upqc(1,i) = qcn
                upql(1,i) = qln
                upqi(1,i) = qin
                upqs(1,i) = qsn
                upth(1,i) = thn
                upLtot(1,i) = Ltot
           else
                upqc(1,i) = 0._rtype
                upql(1,i) = 0._rtype
                upqi(1,i) = 0._rtype
                upqs(1,i) = 0._rtype
                upth(1,i)  = upthl(1,i)
                ! even if there is no condensation Ltot was calculated as a function of T instead of being equal just to lcond 
                upLtot(1,i) = lcond  
           end if

         enddo
                    
         ! Integrate updrafts
         do i=1,mf_nup
           do k=2,nzi

            ! Inverse of Exner function
            ! We should be using the Exner function from SHOC (one calculated in shoc_intr.F). Fix this soon
            iexh = (1.e5_rtype / p(j,k))**(rgas/cp)

             if (tke(j,k-1) == 0.0004_rtype) then
                tke(j,k-1) = 0._rtype
             endif                                       
             ! Original
             if (do_entr_tke) then
                ! New entrainment approximation
                eturb = 1._rtype + mf_c*sqrt(tke(j,k-1))/upw(k-1,i)
                !print*,'real(mf_nup)*upa(k-1,i) = ',real(mf_nup)*upa(k-1,i)
                !eturb = sqrt(1._rtype + 2._rtype*tke(j,k-1)/(real(mf_nup)*upa(k-1,i)*(upw(k-1,i)**2._rtype)))
                entexp = exp(-ent(k-1,i)*eturb*dzt(k-1))
             else
                entexp = exp(-ent(k-1,i)*dzt(k-1))
             end if
             
             entexpu  = exp( -ent(k-1,i)*dzt(k-1)/3._rtype)
             un   = u(j,k-1)  *(1._rtype-entexpu) + upu  (k-1,i)*entexpu
             vn   = v(j,k-1)  *(1._rtype-entexpu) + upv  (k-1,i)*entexpu

             ! If do_precip == true, then include autoconversion in qt and thetal  
             ! following Suselj et al 2019. Equations 15, B8 and B9. 
             if (do_precip) then
               entexpqt = exp( -dzt(k-1)*inv_tau_prec/upw(k-1,i) ) 

               ! Check if H(qc-q0)>0 -> 1 or H(qc-q0) <= 0 -> 0
               if ( (upqc(k-1,i) - acc_thres) > 0._rtype ) then
                 ! Autoconversion occurs
                 qtn  = qt(j,k-1) *(1._rtype-entexp ) + upqt (k-1,i)*entexp &
                        - (upqc(k-1,i) - acc_thres)*(1._rtype - entexpqt) 
                 thln = thl(j,k-1)*(1._rtype-entexp ) + upthl(k-1,i)*entexp &
                        + (upLtot(k-1,i)*iexh/cp)*(upqc(k-1,i) - acc_thres)*(1._rtype - entexpqt)  
                 ! Autoconversion term (eq 15) which will be needed for Rain rate term (eq 18)
                 ! Calculating it here for efficiency (one less loop/if block). Note that Sac(nzi)=0
                 Sac(k-1,i) = - (upqc(k-1,i) - acc_thres)*inv_tau_prec          
               else
                 ! Only entrainment process
                 qtn  = qt(j,k-1) *(1._rtype-entexp ) + upqt (k-1,i)*entexp
                 thln = thl(j,k-1)*(1._rtype-entexp ) + upthl(k-1,i)*entexp
                 Sac(k-1,i) = 0._rtype
                end if
             else
               ! There is no autoconversion, thus only process is entrainment
               qtn  = qt(j,k-1) *(1._rtype-entexp ) + upqt (k-1,i)*entexp
               thln = thl(j,k-1)*(1._rtype-entexp ) + upthl(k-1,i)*entexp
               Sac(k-1,i) = 0._rtype
             endif


             !qtn  = qt(j,k-1) *(1._rtype-entexp ) + upqt (k-1,i)*entexp
             !thln = thl(j,k-1)*(1._rtype-entexp ) + upthl(k-1,i)*entexp
              
             ! Condensation within updrafts, input/output at full levels:
             if (do_condensation) then
                call condensation_mf(qtn, thln, p(j,k), iexh, &
                                    thvn, qcn, thn, qln, qin, qsn, Ltot)
             else
                thvn = thln*(1._rtype+eps*qtn)
                ! MJC (05/16/24) thn is theta, this seems incorrect, as thln is thetal and theta != thetal
                thn = thln                    ! THIS NEEEDS TO BE FIXED!! CALCULATED CORRECTLY
                qcn = 0._rtype
                qin = 0._rtype
                qln = 0._rtype
                qsn = 0._rtype
                Ltot = lcond
             end if

             ! To avoid singularities w equation has to be computed diferently if wp==0
             B=ggr*(0.5_rtype*(thvn+upthv(k-1,i))/thv(j,k-1)-1._rtype)

             if (do_entr_tke) then
             
                if (do_explicit) then
                   
                   enturb = 2._rtype*wb*ent(k-1,i)*dzt(k-1)
                   
                   ! Limiter
                   if (enturb > 1._rtype) then
                      ! we force enturb = 1
                      wn2 = -mf_c*upw(k-1,i)*sqrt(tke(j,k-1))+ 2._rtype*wa*B*dzt(k-1)
                   else
                      wn2 = (upw(k-1,i)**2._rtype)*(1._rtype - enturb)&
                         - enturb*mf_c*upw(k-1,i)*sqrt(tke(j,k-1))+ 2._rtype*wa*B*dzt(k-1)
                   end if
                   
                elseif (do_integral) then
                   wp_integral = wb*ent(k-1,i)*eturb !eturb = (1 + mf_c*sqrt(tke(j,k-1))/upw(k-1,i))
                   
                   if (wp_integral==0._rtype) then
                      wn2 = upw(k-1,i)**2._rtype+2._rtype*wa*B*dzt(k-1)
                   else
                      entw = exp(-2._rtype*wp_integral*dzt(k-1))
                      wn2 = entw*upw(k-1,i)**2._rtype+(wa*B/wp_integral)*(1._rtype-entw)
                   end if 
                
                elseif (do_implicit) then
                   call endscreamrun('integrate_mf: do_implicit=.true. but the implicit w^2 scheme is not implemented. Set edmf_do_implicit=.false. in the namelist.')
                end if
                   
             else
                ! Original code
                wp = wb*ent(k-1,i)
                if (wp==0._rtype) then
                   wn2 = upw(k-1,i)**2._rtype+2._rtype*wa*B*dzt(k-1)
                else
                   entw = exp(-2._rtype*wp*dzt(k-1))
                   wn2 = entw*upw(k-1,i)**2._rtype+wa*B/(wb*ent(k-1,i))*(1._rtype-entw)
                end if   
             end if
    
             if (wn2>0._rtype) then
               upw(k,i)   = sqrt(wn2)
               upthv(k,i) = thvn
               upthl(k,i) = thln
               upqt(k,i)  = qtn
               upqc(k,i)  = qcn
               upu(k,i)   = un
               upv(k,i)   = vn
               upa(k,i)   = upa(k-1,i)
               upth(k,i)  = thn
               upql(k,i)  = qln
               upqi(k,i)  = qin
               upqs(k,i)  = qsn
               upqv(k,i)  = qtn - qcn
               upLtot(k,i) = Ltot
             else
               exit
             end if

           enddo  ! i=1,mf_nup
         enddo    ! k=2,nzi

         ! Calculate Rain Rate term from top to bottom
         ! FIX QS and CHECK IF QS > 0 FIRST TO NOT DIVIDE BY 0
         if (do_precip) then
           do i=1,mf_nup
             do k=nzi-1,1,-1
               if (upqs(k+1,i) > 0._rtype) then
                 Sev(k+1,i) = max( kev*(1._rtype - upqv(k+1,i)/upqs(k+1,i))*sqrt(RRun(k+1,i)), 0._rtype ) 
               else
                 Sev(k+1,i) = 0._rtype
               endif
               ! note that dzt(k) and not dzt(k+1)!
               ! RRun(top) = 0 (initialized at the beginning of subroutine)
               RRun(k,i)  = RRun(k+1,i) + (zi(j,k)-zi(j,k+1))*rho_zi(j,k)*( Sac(k+1,i)*(1._rtype-frain_det) + Sev(k+1,i) )
               ! Limiter in case RRun < 0 due to strong evaporation
               if (RRun(k,i) < 0._rtype) then
                RRun(k,i) = 0._rtype
               endif
               Smel(k,i)  = RRun(k+1,i)*iexh*(upLtot(k+1,i)-upLtot(k,i))/(dzt(k)*cp*rho_zi(j,k))
             enddo
           enddo
         endif

         ! Calculate the total source terms (eq 20 and 21)
         if (do_precip) then
          do i=1,mf_nup
           do k=1,nzi
            Sthl(k,i) = - upLtot(k,i)*iexh*( Sac(k,i) + Sev(k,i) )/cp + Smel(k,i)
            Sqt(k,i)  = Sac(k,i) + Sev(k,i)
           enddo
          enddo
         endif
 

         ! writing updraft properties for output
         ! all variables, except areas (moist_a and dry_a) are now multipled by the area
         do k=1,nzi

           ! first sum over all i-updrafts
           do i=1,mf_nup
             if (upqc(k,i) > 0._rtype) then
               moist_a(j,k)   = moist_a(j,k)   + upa(k,i)
               moist_w(j,k)   = moist_w(j,k)   + upa(k,i)*upw(k,i)
               moist_qt(j,k)  = moist_qt(j,k)  + upa(k,i)*upqt(k,i)
               moist_thl(j,k) = moist_thl(j,k) + upa(k,i)*upthl(k,i)
               moist_th(j,k)  = moist_th(j,k)  + upa(k,i)*upth(k,i)
               moist_u(j,k)   = moist_u(j,k)   + upa(k,i)*upu(k,i)
               moist_v(j,k)   = moist_v(j,k)   + upa(k,i)*upv(k,i)
               moist_qc(j,k)  = moist_qc(j,k)  + upa(k,i)*upqc(k,i)
             else
               dry_a(j,k)     = dry_a(j,k)     + upa(k,i)
               dry_w(j,k)     = dry_w(j,k)     + upa(k,i)*upw(k,i)
               dry_qt(j,k)    = dry_qt(j,k)    + upa(k,i)*upqt(k,i)
               dry_thl(j,k)   = dry_thl(j,k)   + upa(k,i)*upthl(k,i)
               dry_th(j,k)    = dry_th(j,k)    + upa(k,i)*upth(k,i)
               dry_u(j,k)     = dry_u(j,k)     + upa(k,i)*upu(k,i)
               dry_v(j,k)     = dry_v(j,k)     + upa(k,i)*upv(k,i)
             endif
           enddo

           if ( dry_a(j,k) > 0._rtype ) then
             dry_w(j,k)   = dry_w(j,k)   / dry_a(j,k)
             dry_qt(j,k)  = dry_qt(j,k)  / dry_a(j,k)
             dry_thl(j,k) = dry_thl(j,k) / dry_a(j,k)
             dry_th(j,k)  = dry_th(j,k)  / dry_a(j,k)
             dry_u(j,k)   = dry_u(j,k)   / dry_a(j,k)
             dry_v(j,k)   = dry_v(j,k)   / dry_a(j,k)
           else
             dry_w(j,k)   = 0._rtype
             dry_qt(j,k)  = 0._rtype
             dry_thl(j,k) = 0._rtype
             dry_th(j,k)  = 0._rtype
             dry_u(j,k)   = 0._rtype
             dry_v(j,k)   = 0._rtype
           endif

           if ( moist_a(j,k) > 0._rtype ) then
             moist_w(j,k)   = moist_w(j,k)   / moist_a(j,k)
             moist_qt(j,k)  = moist_qt(j,k)  / moist_a(j,k)
             moist_thl(j,k) = moist_thl(j,k) / moist_a(j,k)
             moist_th(j,k)  = moist_th(j,k)  / moist_a(j,k)
             moist_u(j,k)   = moist_u(j,k)   / moist_a(j,k)
             moist_v(j,k)   = moist_v(j,k)   / moist_a(j,k)
             moist_qc(j,k)  = moist_qc(j,k)  / moist_a(j,k)
           else
             moist_w(j,k)   = 0._rtype
             moist_qt(j,k)  = 0._rtype
             moist_thl(j,k) = 0._rtype
             moist_th(j,k)  = 0._rtype
             moist_u(j,k)   = 0._rtype
             moist_v(j,k)   = 0._rtype
             moist_qc(j,k)  = 0._rtype
           endif
             
         enddo

         do k=1,nzi
           do i=1,mf_nup
             ae  (j,k) = ae  (j,k) - upa(k,i)
             aw  (j,k) = aw  (j,k) + upa(k,i)*upw(k,i)
             awu (j,k) = awu (j,k) + upa(k,i)*upw(k,i)*upu(k,i)
             awv (j,k) = awv (j,k) + upa(k,i)*upw(k,i)*upv(k,i)
             awthv(j,k)= awthv(j,k)+ upa(k,i)*upw(k,i)*upthv(k,i)
             awthl(j,k)= awthl(j,k)+ upa(k,i)*upw(k,i)*upthl(k,i) !*cpair/iexh
             awth(j,k) = awth(j,k) + upa(k,i)*upw(k,i)*upth(k,i) !*cpair/iexh
             awqt(j,k) = awqt(j,k) + upa(k,i)*upw(k,i)*upqt(k,i)
             awqc(j,k) = awqc(j,k) + upa(k,i)*upw(k,i)*upqc(k,i)             
             awqv(j,k) = awqv(j,k) + upa(k,i)*upw(k,i)*upqv(k,i)
             awql(j,k) = awql(j,k) + upa(k,i)*upw(k,i)*upql(k,i)
             awqi(j,k) = awqi(j,k) + upa(k,i)*upw(k,i)*upqi(k,i)
             auSthl(j,k) = auSthl(j,k) + upa(k,i)*Sthl(k,i)
             auSqt(j,k)  = auSqt(j,k)  + upa(k,i)*Sqt(k,i)
             auRRun(j,k) = auRRun(j,k)  + upa(k,i)*RRun(k,i)
           enddo
         enddo

         ! MJC: Calculate cold pool convective velocity scale
         ! it's missing a_cp
         if (do_precip) then
           !mf_w_cp(j) = ( mf_w_cp_old + auRRun(j,1)*dt*mf_a_wcp ) / (1 + dt/mf_tau_wcp)  ! tau_cp = 4*3600 = 14400._rtype
          mf_w_cp(j) = mf_w_cp_help(nstep_mff_help);
          ! (auRRun(j,1)*0.001) converts auRRun from mm/s to meters/s 
          mf_a_wcp_local =  (1.0_rtype + dt/(mf_tau_wcp*3600.0_rtype))/(dt*2.778e-7_rtype)
          !mf_w_cp_help(nstep_mff_help+1) = ( mf_w_cp_help(nstep_mff_help) + (auRRun(j,1)*0.001_rtype)*dt*mf_a_wcp_local ) / (1.0_rtype + dt/(mf_tau_wcp*3600.0_rtype))  ! 
          !mf_w_cp_help(nstep_mff_help+1) = ( mf_w_cp_help(nstep_mff_help) + (auRRun(j,1)*0.001_rtype)*dt/(mf_tau_wcp*3600.0_rtype) ) / (1 + dt/(mf_tau_wcp*3600.0_rtype))  ! tau_cp = 4*3600 = 14400._rtype
          
          !mf_w_cp_help(nstep_mff_help+1) = ( mf_w_cp_help(nstep_mff_help) + (auRRun(j,1)*0.001_rtype)*dt*mf_a_wcp ) / (1.0_rtype + dt/(mf_tau_wcp*3600.0_rtype))  ! tau_cp = 4*3600 = 14400._rtype 
          mf_w_cp_help(nstep_mff_help+1) = ( mf_w_cp_help(nstep_mff_help) + auRRun(j,1)*dt ) / (1.0_rtype + dt/(mf_tau_wcp*3600.0_rtype))  ! tau_cp = 4*3600 = 14400._rtype

         endif 

         ! Find highest vertical level where the plume ensemble is dry (i.e., moist_qc = 0)
         check  = .true. 
         plume_dry_height(j) = 0._rtype
         index_top = 1
         do k=1,nzi
            if (check .and. aw(j,k) > 0._rtype .and. moist_qc(j,k) .EQ. 0._rtype) then
               plume_dry_height(j) = zi(j,k)
               index_top = k
            else
               check = .false.  
            endif
         enddo
         
         ! Highest vertical level reached by the moist plumes 
         !(analyzed from the top of the dry CBL given by plume_dry_height)
         check  = .true. 
         plume_top_height(j) = 0._rtype
         do k=index_top+1,nzi
            if (check .and. aw(j,k) > 0._rtype .and. moist_qc(j,k) > 0._rtype) then
               plume_top_height(j) = zi(j,k)
            else
               check = .false.  
            endif
         enddo
        
         !print*,'plume_dry_height = ',plume_dry_height(j)
         !print*,'plume_top_height = ',plume_top_height(j)

         ! Check CLF condition on mass-flux (aw)
         do k=1,nz
           cfl_zt = (2._rtype/dt)*rho_zt(j,k)*dz_zt(j,k)
           if (zi(j,k) < ztop(j)*1.5_rtype) then
              if (aw(j,k) > (cfl_zt/rho_zi(j,k)) ) then
                 !print*,'WARNING: aw > CFL'
                 !print*,'aw(j,k) = ',aw(j,k)
                 !print*,'CFL = ',cfl_zt/rho_zi(j,k)
                 !print*,'k index = ',k
              endif
                            
           endif
           cfl(j,k) = cfl_zt/rho_zi(j,k)   
         enddo
         
       end if  ! ( wthv > 0.0 )

       if (ANY(dry_a  (j,:)>0._rtype)) freq_dry(j)   = 1._rtype
       if (ANY(moist_a(j,:)>0._rtype)) freq_moist(j) = 1._rtype
     end do ! j=1,shcol

     ! flip output variables so index 1 = model top (i.e. lowest pressure)
     do k=1,nzi
       dry_a_out(:,nzi-k+1) = dry_a(:,k)
       dry_w_out(:,nzi-k+1) = dry_w(:,k)
       dry_qt_out(:,nzi-k+1) = dry_qt(:,k)
       dry_thl_out(:,nzi-k+1) = dry_thl(:,k)
       dry_u_out(:,nzi-k+1) = dry_u(:,k)
       dry_v_out(:,nzi-k+1) = dry_v(:,k)

       moist_a_out(:,nzi-k+1) = moist_a(:,k)
       moist_w_out(:,nzi-k+1) = moist_w(:,k)
       moist_qt_out(:,nzi-k+1) = moist_qt(:,k)
       moist_thl_out(:,nzi-k+1) = moist_thl(:,k)
       moist_u_out(:,nzi-k+1) = moist_u(:,k)
       moist_v_out(:,nzi-k+1) = moist_v(:,k)
       moist_qc_out(:,nzi-k+1) = moist_qc(:,k)

       ae_out(:,nzi-k+1) = ae(:,k)
       aw_out(:,nzi-k+1) = aw(:,k)
       awthv_out(:,nzi-k+1) = awthv(:,k)
       awthl_out(:,nzi-k+1) = awthl(:,k)
       awqt_out(:,nzi-k+1) = awqt(:,k)
       awql_out(:,nzi-k+1) = awql(:,k)
       awqi_out(:,nzi-k+1) = awqi(:,k)
       awu_out(:,nzi-k+1) = awu(:,k)
       awv_out(:,nzi-k+1) = awv(:,k)
       auSthl_out(:,nzi-k+1) =  auSthl(:,k) 
       auSqt_out(:,nzi-k+1) =  auSqt(:,k) 
       auRRun_out(:,nzi-k+1) =  auRRun(:,k) 

       !thlflx_out(:,nzi-k+1) = thlflx(:,k)
       !qtflx_out(:,nzi-k+1) = qtflx(:,k)
     end do


  end subroutine integrate_mf



                       
  subroutine oneplume( nz, nzi, zi, dzt, ent, p, qt, thl, thv,   &
                       thv_zi, wmax, wmin, sigmaw, sigmaqt, sigmath, wa, wb, &
                       do_condensation, plumeheight )
  !**********************************************************************
  ! Calculate a single plume with zero entrainment
  ! to be used for a dynamic mixing length calculation
  ! By Rachel Storer
  !**********************************************************************

    integer,  intent(in)                    :: nz, nzi

    real(rtype), intent(in)                 :: wmax, wmin, sigmaw, sigmaqt, sigmath, wa, wb

    real(rtype), dimension(nz),  intent(in) ::  dzt, qt, thl, thv, ent
    real(rtype), dimension(nzi), intent(in) ::  zi, p, thv_zi

                                                      
    logical, intent(in)                       :: do_condensation

    real(rtype), intent(inout) :: plumeheight
     
    !local variables
    integer                        :: k
    real(rtype)                    :: thvn, qtn, thln, qcn, thn, qln, qin, qsn, wn2
    real(rtype)                    :: Ltot ! this is an output of condensation_mf() but not used in oneplume()
    real(rtype)                    :: iexh, entexp, entexpu, wp, entw
    real(rtype), dimension(nzi)    :: upw, upa, upqt, upthv, upthl, upth, &
                                      upqc, upql, upqi, b, thvflx
                                      
        

    thvflx  = 0._rtype
    b     = 0._rtype
    upw   = 0._rtype
    upthl = 0._rtype
    upthv = 0._rtype
    upqt  = 0._rtype
    upa   = 0._rtype
    upqc  = 0._rtype
    upth  = 0._rtype
    upql  = 0._rtype
    upqi  = 0._rtype
    Ltot  = 0._rtype

    upw(1) = 0.5_rtype * (wmax+wmin)
    upa(1) = 0.5_rtype * erf( wmax/(sqrt(2.5_rtype)*sigmaw) ) &
                      - 0.5_rtype * erf( wmin/(sqrt(2._rtype)*sigmaw) )
  
    upqt(1)  = qt(1)  + 0.32_rtype * upw(1) * sigmaqt/sigmaw
    upthv(1) = thv(1) + 0.58_rtype * upw(1) * sigmath/sigmaw
           
    upthl(1) = upthv(1) / (1._rtype+eps*upqt(1))
    upth(1)  = upthl(1)
  
    ! get cloud, lowest momentum level 
    if (do_condensation) then
      iexh = (1.e5_rtype / p(1))**(rgas/cp)
      call condensation_mf(upqt(1), upthl(1), p(1), iexh, &
                           thvn, qcn, thn, qln, qin, qsn, Ltot)
      upthv(1) = thvn
      upqc(1)  = qcn
      upql(1)  = qln
      upqi(1)  = qin
      upth(1)  = thn
    else
      ! assume no cldliq
      upthv(1) = upthl(1)*(1._rtype+eps*upqt(1))
      upth(1)  = upthl(1)

    end if
  
    do k=2,nzi
   
      entexp  = exp(-ent(k-1)*dzt(k-1))
      entexpu = exp(-ent(k-1)*dzt(k-1)/3._rtype)
              
      ! integrate updraft
      qtn  = qt(k-1) *(1._rtype-entexp ) + upqt (k-1)*entexp
      thln = thl(k-1)*(1._rtype-entexp ) + upthl(k-1)*entexp
                      
      ! get cloud, momentum levels
      if (do_condensation) then
        iexh = (1.e5_rtype / p(k))**(rgas/cp)
        call condensation_mf(qtn, thln, p(k), iexh, &
                             thvn, qcn, thn, qln, qin, qsn, Ltot)
      else
        thvn = thln*(1._rtype+eps*qtn)
        thn = thln                       ! THIS NEEEDS TO BE FIXED!! CALCULATED CORRECTLY
        qcn = 0._rtype
        qin = 0._rtype
        qln = 0._rtype
      end if
      ! get buoyancy
      b(k)=ggr*(0.5_rtype*(thvn+upthv(k-1))/thv(k-1)-1._rtype)
      wp = wb*ent(k-1)
      if (wp==0._rtype) then
         wn2 = upw(k-1)**2._rtype+2._rtype*wa*b(k)*dzt(k-1)
      else
         entw = exp(-2._rtype*wp*dzt(k-1))
         wn2 = entw*upw(k-1)**2._rtype+wa*b(k)/(wb*ent(k-1))*(1._rtype-entw)
      endif   
  
      if (wn2>0._rtype) then
        upw(k)   = sqrt(wn2)
        upthv(k) = thvn
        upthl(k) = thln
        upqt(k)  = qtn
        upqc(k)  = qcn
        upa(k)   = upa(k-1)
        upql(k)  = qln
        upqi(k)  = qin
        upth(k)  = thn
        plumeheight = zi(k)
      else
        exit
      end if
      
    enddo
     

  end subroutine oneplume
  


  subroutine condensation_mf( qt, thl, p, iex, thv, qc, th, ql, qi, qs, Ltot)
  !
  ! zero or one condensation for edmf: calculates thv and qc
  !
      ! use wv_saturation,      only : qsat

       real(rtype),intent(in) :: qt,thl,p,iex
       real(rtype),intent(out):: thv,qc,th,ql,qi,qs,Ltot

       !local variables
       integer :: niter,i
       real(rtype) :: diff,t,qcold,es
       ! ice_wt is the ICE weight returned by get_Ltot_rl:
       !   ice_wt = 0 at T >= 273.15 K (all liquid)
       !   ice_wt = 1 at T <= 253.15 K (all ice)
       !   linear in between. Liquid fraction is (1 - ice_wt).
       real(rtype) :: ice_wt

       ! max number of iterations
       niter=50
       ! minimum difference
       diff=2.e-5_rtype

       qc=0._rtype
       t=thl/iex

  !by definition:
  ! T   = Th*Exner, Exner=(p/p0)^(R/cp)   (1)
  ! Thl = Th - L/cp*ql/Exner              (2)
  !so:
  ! Th  = Thl + L/cp*ql/Exner             (3)
  ! T   = Th*Exner=(Thl+L/cp*ql/Exner)*Exner    (4)
  !     = Thl*Exner + L/cp*ql
       do i=1,niter
         call get_Ltot_rl(t,Ltot,ice_wt)
         t = thl/iex+Ltot/cp*qc   !as in (4)

         ! qsat, p is in pascal (check!)
         call qsat(t,p,es,qs)
         qcold = qc
         qc = max(0.5_rtype*qc+0.5_rtype*(qt-qs),0._rtype)
         if (abs(qc-qcold)<diff) exit
       enddo

       ! Update T and calculate the other variables
       call get_Ltot_rl(t,Ltot,ice_wt)
       t = thl/iex+Ltot/cp*qc
       thv = (thl+Ltot/cp*iex*qc)*(1.+eps*(qt-qc)-qc)
       th = t*iex
       ! ice_wt returned from get_Ltot_rl is the ice fraction (1 at cold T, 0 at warm T).
       ! Split total condensate qc into liquid and ice accordingly.
       qi = qc*ice_wt
       ql = qc*(1._rtype - ice_wt)

       ! Save qs needed for precip evaporation Sev term
       call qsat(t,p,es,qs)

  end subroutine condensation_mf
  ! MJ/Kay: The subroutine get_Ltot_rl outputs the total latent heat of vaporization 
  !i.e. taking into account liquid and ice phases if temperature is below freezing.
  ! This subroutine was adapted from the one in wv_saturation.F90 called calc_hltalt()    
  subroutine get_Ltot_rl(t, hltalt, weight)
    !------------------------------------------------------------------!
    ! Purpose:                                                         !
    !   Calculate latent heat of vaporization of water at a given      !
    !   temperature, taking into account the ice phase if temperature  !
    !   is below freezing.                                             !
    !   Optional argument also calculates a term used to calculate     !
    !   d(es)/dT within the water-ice transition range.                !
    !------------------------------------------------------------------!
  
    ! Inputs
    real(rtype), intent(in) :: t        ! Temperature
    ! Outputs
    real(rtype), intent(out) :: hltalt  ! Appropriately modified hlat
    real(rtype), intent(out) :: weight  ! Weight for es transition from water to ice

  
    ! Local variables
    real(rtype) :: tc      ! Temperature in degrees C
    real(rtype) :: ttrice = 20.00_rtype  ! transition range from es over H2O to es over ice in C
    real(rtype) :: tmelting = 273.15_rtype ! freezing T of fresh water          ~ K

    ! Loop iterator
    integer :: i
    
    weight = 0.0_rtype
    ! At the top of SHOCs module we added: use wv_saturation, only : qsat, no_ip_hltalt
    ! Subroutine no_ip_hltalt calculates latent heat of vaporization of pure liquid water at a given T
    !instead of just assuming the latent heat of evaporation at 100C
    call no_ip_hltalt(t,hltalt)

    if (t < tmelting) then
       ! Weighting of hlat accounts for transition from water to ice.
       tc = t - tmelting
  
       if (tc >= -ttrice) then
          weight = -tc/ttrice  
       else
          weight = 1.0_rtype
       end if
  
       hltalt = hltalt + weight*lice 
  
    end if
  
  end subroutine get_Ltot_rl

  subroutine calc_mf_vertflux(shcol,nlev,nlevi,aw,awvar,var,var_zi,varflx)

    implicit none

  ! INPUT VARIABLES
    ! number of SHOC columns
    integer, intent(in) :: shcol
    ! number of midpoint levels
    integer, intent(in) :: nlev
    ! number of interface levels
    integer, intent(in) :: nlevi
    ! Sum plume (a_i*w_i) [m/s]
    real(rtype), intent(in) :: aw(shcol,nlevi)
    ! Sum plume vertical flux of generic variable var (a_i*w_i*var_i) [units vary]
    real(rtype), intent(in) :: awvar(shcol,nlevi)
    ! Input variable on thermo/full grid [units vary]
    real(rtype), intent(in) :: var_zi(shcol,nlevi) ! NOTE: var is interpolated to zi, so has dim nzi
    real(rtype), intent(in) :: var(shcol,nlev)

  ! OUTPUT VARIABLE
    real(rtype), intent(out) :: varflx(shcol,nlevi)

  ! INTERNAL VARIABLES
    integer :: i,k

    ! MKW TODO: SHOC has separate subroutines for lower (k=nlevi) and
    !   upper (k=1) boundary conditions. Make these later if SCREAM
    !   folks want that. Should be very quick.

    ! diagnose MF fluxes
    varflx(:shcol,1) = 0._rtype
    do k=2,nlev
      do i=1,shcol
        varflx(i,k)= awvar(i,k) - aw(i,k)*0.5*(var(i,k-1)+var(i,k)) ! centered differences
        !varflx(i,k)= awvar(i,k) - aw(i,k)*var(i,k) ! upwind scheme (in reference to the surface)
        !varflx(i,k)= awvar(i,k) - aw(i,k)*var(i,k-1) ! downwind scheme (in reference to the surface)
        !varflx(i,k)= awvar(i,k) - aw(i,k)*var_zi(i,k)

      end do
    end do
    varflx(:shcol,nlevi) = 0._rtype
    
  end subroutine calc_mf_vertflux

  subroutine compute_tmpi3(nlevi, shcol, dtime, rho_zi, tmpi3)

    !intent-ins
    integer,     intent(in) :: nlevi, shcol
    !time step [s]
    real(rtype), intent(in) :: dtime
    !air density at interfaces [kg/m3]
    real(rtype), intent(in) :: rho_zi(shcol,nlevi)

    !intent-out
    real(rtype), intent(out) :: tmpi3(shcol,nlevi)

    !local vars
    integer :: i, k

    tmpi3(:,1) = 0._rtype
    ! eqn: tmpi3 = dt*g*rho
    do k = 2, nlevi
      do i = 1, shcol
         tmpi3(i,k) = dtime *  ggr*rho_zi(i,k)
      enddo
    enddo

  end subroutine compute_tmpi3


!==============================================================
! Define grid variables needed for the parameterization

subroutine shoc_grid( &
          shcol,nlev,nlevi,&           ! Input
          zt_grid,zi_grid,pdel,&       ! Input
          dz_zt,dz_zi,rho_zt)          ! Output

  ! Purpose of this subroutine is to define the thickness
  !  arrays of each column, to be used for finite differencing
  !  throughout the SHOC parameterization, also define air
  !  density in SHOC

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_grid_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of columns [-]
  integer, intent(in) :: shcol
  ! number of mid-point levels [-]
  integer, intent(in) :: nlev
  ! number of interface levels [-]
  integer, intent(in) :: nlevi
  ! mid-point grid heights [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! interface grid heights [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! pressure differences centered on mid-point grid [Pa]
  real(rtype), intent(in) :: pdel(shcol,nlev)

! OUTPUT VARIABLES
  ! thickness (dz) on the thermo grid [m]
  real(rtype), intent(out) :: dz_zt(shcol,nlev)
  ! thickness (dz) on the interface grid [m]
  real(rtype), intent(out) :: dz_zi(shcol,nlevi)
  ! air density on the thermo grid [kg/m3]
  real(rtype), intent(out) :: rho_zt(shcol,nlev)

  ! local variables
  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call shoc_grid_f(shcol,nlev,nlevi,&
                     zt_grid,zi_grid,pdel,&
                     dz_zt,dz_zi,rho_zt)
     return
  endif
#endif

  do k=1,nlev
    do i=1,shcol
      ! define thickness of the thermodynamic gridpoints
      dz_zt(i,k) = zi_grid(i,k) - zi_grid(i,k+1)

      ! define thickness of the interface grid points
      if (k .eq. 1) then
        dz_zi(i,k) = 0._rtype ! never used
      else
        dz_zi(i,k) = zt_grid(i,k-1) - zt_grid(i,k)
      endif

      ! Define the air density on the thermo grid
      rho_zt(i,k) = (1._rtype/ggr)*(pdel(i,k)/dz_zt(i,k))

    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)

  ! Set lower condition for dz_zi
  dz_zi(:shcol,nlevi) = zt_grid(:shcol,nlev)

  return

end subroutine shoc_grid

!==============================================================
! Compute vapor from SHOC prognostic/diagnostic variables

subroutine compute_shoc_vapor( &
          shcol,nlev,qw,ql,&           ! Input
          qv)                          ! Output

  ! Purpose of this subroutine is to compute water vapor
  !   based on SHOC's prognostic total water mixing ratio
  !   and diagnostic cloud water mixing ratio.

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_shoc_vapor_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of columns [-]
  integer, intent(in) :: shcol
  ! number of mid-point levels [-]
  integer, intent(in) :: nlev
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: qw(shcol,nlev)
  ! cloud water mixing ratio [kg/kg]
  real(rtype), intent(in) :: ql(shcol,nlev)

! OUTPUT VARIABLES
  ! water vapor mixing ratio [kg/kg]
  real(rtype), intent(out) :: qv(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call compute_shoc_vapor_f(shcol,nlev,qw,ql,qv)
     return
  endif
#endif

  do k = 1, nlev
    do i = 1, shcol
      qv(i,k) = qw(i,k) - ql(i,k)
    enddo
  enddo

  return

end subroutine compute_shoc_vapor

!==============================================================
! Compute temperature from SHOC prognostic/diagnostic variables

subroutine compute_shoc_temperature( &
          shcol,nlev,thetal,ql,inv_exner,& ! Input
          tabs)                            ! Output

  ! Purpose of this subroutine is to compute temperature
  !   based on SHOC's prognostic liquid water potential
  !   temperature.

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_shoc_temperature_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of columns [-]
  integer, intent(in) :: shcol
  ! number of mid-point levels [-]
  integer, intent(in) :: nlev
  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! cloud water mixing ratio [kg/kg]
  real(rtype), intent(in) :: ql(shcol,nlev)
  ! inverse exner function [-]
  real(rtype), intent(in) :: inv_exner(shcol,nlev)

! OUTPUT VARIABLES
  ! absolute temperature [K]
  real(rtype), intent(out) :: tabs(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
     call compute_shoc_temperature_f(shcol,nlev,thetal,ql,inv_exner,tabs)
     return
  endif
#endif

  do k = 1, nlev
    do i = 1, shcol
      tabs(i,k) = thetal(i,k)/inv_exner(i,k)+(lcond/cp)*ql(i,k)
    enddo
  enddo

  return

end subroutine compute_shoc_temperature

!==============================================================
! Update T, q, tracers, tke, u, and v based on implicit diffusion
! Here we use a backward Euler scheme.

subroutine update_prognostics_implicit( &
         shcol,nlev,nlevi,num_tracer,&    ! Input
         dtime,dz_zt,dz_zi,rho_zt,&       ! Input
         zt_grid,zi_grid,tk,tkh,&         ! Input
         uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,& ! Input
         wtracer_sfc,&                    ! Input
         do_mf,mf_ae,mf_aw,mf_awu,mf_awv,&! EDMF Input
         mf_awthl,mf_awqt,&               ! EDMF Input         
         thetal,qw,tracer,tke,&           ! Input/Output
         u_wind,v_wind)                   ! Input/Output

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: update_prognostics_implicit_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of vertical levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! number of tracers
  integer, intent(in) :: num_tracer

  ! SHOC timestep [s]
  real(rtype), intent(in) :: dtime
  ! Eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)
  ! Eddy coefficient for heat [m2/s]
  real(rtype), intent(in) :: tkh(shcol,nlev)
  ! Air density on thermo grid [kg/m3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)
  ! height thickness centered on thermo grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! height thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! vertical zonal momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! vertical meridional momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! vertical heat flux at surface [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! vertical moisture flux at surface [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! vertical tracer flux at surface [units vary m/s]
  real(rtype), intent(in) :: wtracer_sfc(shcol,num_tracer)
  ! heights of mid-point [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights at interfaces [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  
  ! MJC: EDMF inputs
  ! If .true., diagnose MF plumes and include in vertical diffusion solver
  logical, intent(in) ::  do_mf
  ! Fractional area of of nonconvective environment
  real(rtype), intent(in) :: mf_ae(shcol,nlevi)
  ! MF plume sum(a_i*w_i)
  real(rtype), intent(in) :: mf_aw(shcol,nlevi)
  ! Total MF plume u turbulent flux
  real(rtype), intent(in) :: mf_awu(shcol,nlevi)
  ! Total MF plume v turbulent flux
  real(rtype), intent(in) :: mf_awv(shcol,nlevi)
  ! Total MF plume theta_l turbulent flux
  real(rtype), intent(in) :: mf_awthl(shcol,nlevi)
  ! Total MF plume q_t turbulent flux
  real(rtype), intent(in) :: mf_awqt(shcol,nlevi)

! IN/OUT VARIABLES
  ! liquid water potential temperature [K]
  real(rtype), intent(inout) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(inout) :: qw(shcol,nlev)
  ! tracers [varies]
  real(rtype), intent(inout) :: tracer(shcol,nlev,num_tracer)
  ! zonal wind [m/s]
  real(rtype), intent(inout) :: u_wind(shcol,nlev)
  ! meridional wind [m/s]
  real(rtype), intent(inout) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)

! LOCAL VARIABLES
  integer     :: p
  real(rtype) :: rdp_zt(shcol,nlev)
  real(rtype) :: tmpi(shcol,nlevi)
  real(rtype) :: tkh_zi(shcol,nlevi)
  real(rtype) :: tk_zi(shcol,nlevi)
  real(rtype) :: rho_zi(shcol,nlevi)

  real(rtype) :: flux_dummy(shcol)
  real(rtype) :: ksrf(shcol), wtke_sfc(shcol)

!  real(rtype) :: du(shcol,nlev) ! Superdiagonal for solver
!  real(rtype) :: dl(shcol,nlev) ! Factorized subdiagonal for solver
!  real(rtype) :: d(shcol,nlev)  ! Factorized diagonal for solver

! MJC: Variables for previous diffusion solver method (SCREAMv0)
  real(rtype) :: ca(shcol,nlev) ! superdiagonal for solver
  real(rtype) :: cc(shcol,nlev) ! subdiagonal for solver
  real(rtype) :: denom(shcol,nlev) ! denominator in solver
  real(rtype) :: ze(shcol,nlev)
! MJC: For EDMF in diffusion solver
  real(rtype) :: tmpi3(shcol,nlevi) 

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call update_prognostics_implicit_f(&
           shcol,nlev,nlevi,num_tracer,&    ! Input
           dtime,dz_zt,dz_zi,rho_zt,&       ! Input
           zt_grid,zi_grid,tk,tkh,&         ! Input
           uw_sfc,vw_sfc,wthl_sfc,wqw_sfc,& ! Input
           wtracer_sfc,&                    ! Input
           thetal,qw,tracer,tke,&           ! Input/Output
           u_wind,v_wind)                   ! Input/Output
     return
  endif
#endif

  ! linearly interpolate tkh, tk, and air density onto the interface grids
  call linear_interp(zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,shcol,0._rtype)

  tmpi(:,1) = 0._rtype
  ! Define the tmpi variable, which is really dt*(g*rho)**2/dp
  !  at interfaces. Substitue dp = g*rho*dz in the above equation
  call compute_tmpi(nlevi, shcol, dtime, rho_zi, dz_zi, tmpi)

  ! For MF component
  tmpi3(:,1) = 0._rtype
  call compute_tmpi3(nlevi, shcol, dtime, rho_zi, tmpi3)
  
  ! compute 1/dp term, needed in diffusion solver
  call dp_inverse(nlev, shcol, rho_zt, dz_zt, rdp_zt)

  ! compute terms needed for the implicit surface stress (ksrf)
  ksrf(1:shcol)      = impli_srf_stress_term(shcol, rho_zi(:,nlevi), &
                         uw_sfc, vw_sfc, u_wind(:,nlev), v_wind(:,nlev))

  !compute term needed for tke flux calc (wtke_sfc)
  wtke_sfc(1:shcol) = tke_srf_flux_term(shcol, uw_sfc, vw_sfc)

  ! compute surface fluxes for liq. potential temp, water and tke
  call sfc_fluxes(shcol, num_tracer, dtime, rho_zi(:,nlevi), rdp_zt(:,nlev), &
                  wthl_sfc, wqw_sfc, wtke_sfc, wtracer_sfc, &
                  thetal(:,nlev), qw(:,nlev), tke(:,nlev), tracer(:,nlev,:))

! MJC: Commenting this code that uses the Thomas factorization method.
!  ! Call decomp for momentum variables
!  call vd_shoc_decomp(shcol,nlev,nlevi,tk_zi,tmpi,rdp_zt,dtime,&
!     ksrf,du,dl,d)
!
!  ! march u_wind one step forward using implicit solver
!  call vd_shoc_solve(shcol,nlev,du,dl,d,u_wind)
!
!  ! march v_wind one step forward using implicit solver
!  call vd_shoc_solve(shcol,nlev,du,dl,d,v_wind)
!
!  ! Call decomp for thermo variables
!  flux_dummy(:) = 0._rtype ! fluxes applied explicitly, so zero fluxes out
!                           ! for implicit solver decomposition
!  call vd_shoc_decomp(shcol,nlev,nlevi,tkh_zi,tmpi,rdp_zt,dtime,&
!     flux_dummy,du,dl,d)
!
!  ! march temperature one step forward using implicit solver
!  call vd_shoc_solve(shcol,nlev,du,dl,d,thetal)
!
!  ! march total water one step forward using implicit solver
!  call vd_shoc_solve(shcol,nlev,du,dl,d,qw)
!
!  ! march tke one step forward using implicit solver
!  call vd_shoc_solve(shcol,nlev,du,dl,d,tke)
!
!  ! march tracers one step forward using implicit solver
!  do p=1,num_tracer
!    call vd_shoc_solve(shcol,nlev,du,dl,d,tracer(:shcol,:nlev,p))
!  enddo


! Call decomp for momentum variables
  call vd_shoc_decomp(shcol,nlev,nlevi,tk_zi,tmpi,rdp_zt,dtime,&
     ksrf,.false.,mf_ae,mf_aw,tmpi3,ca,cc,denom,ze)

  ! march u_wind one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,.false.,mf_awu,tmpi3,rdp_zt,u_wind)

  ! march v_wind one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,.false.,mf_awv,tmpi3,rdp_zt,v_wind)

! Call decomp for thermo variables
  flux_dummy(:) = 0._rtype ! fluxes applied explicitly, so zero fluxes out
                           ! for implicit solver decomposition
  call vd_shoc_decomp(shcol,nlev,nlevi,tkh_zi,tmpi,rdp_zt,dtime,&
     flux_dummy,do_mf,mf_ae,mf_aw,tmpi3,ca,cc,denom,ze)
     
  ! march temperature one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,do_mf,mf_awthl,tmpi3,rdp_zt,thetal)

  ! march total water one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,do_mf,mf_awqt,tmpi3,rdp_zt,qw)

  ! MKW: Call decomp one more time for TKE and tracers so they don't "see" MF plumes
  call vd_shoc_decomp(shcol,nlev,nlevi,tkh_zi,tmpi,rdp_zt,dtime,&
          flux_dummy,.false.,mf_ae,mf_aw,tmpi3,ca,cc,denom,ze)

  ! march tke one step forward using implicit solver
  call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,.false.,mf_aw,tmpi3,rdp_zt,tke)

  ! march tracers one step forward using implicit solver
  do p=1,num_tracer
    call vd_shoc_solve(shcol,nlev,nlevi,ca,cc,denom,ze,.false.,mf_aw,tmpi3,rdp_zt,tracer(:shcol,:nlev,p))
  enddo
  
  return

end subroutine update_prognostics_implicit

subroutine compute_tmpi(nlevi, shcol, dtime, rho_zi, dz_zi, tmpi)

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_tmpi_f
#endif

  !intent-ins
  integer,     intent(in) :: nlevi, shcol
  !time step [s]
  real(rtype), intent(in) :: dtime
  !air density at interfaces [kg/m3]
  real(rtype), intent(in) :: rho_zi(shcol,nlevi)
  !height thickness at interfaces [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)

  !intent-out
  real(rtype), intent(out) :: tmpi(shcol,nlevi)

  !local vars
  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call compute_tmpi_f(nlevi, shcol, dtime, rho_zi, dz_zi, tmpi)
     return
  endif
#endif

  tmpi(:,1) = 0._rtype
  ! eqn: tmpi = dt*(g*rho)**2/dp, where dp = g*rho*dz, therefore tmpi = dt*g*rho/dz
  do k = 2, nlevi
    do i = 1, shcol
       tmpi(i,k) = dtime * (ggr*rho_zi(i,k)) / dz_zi(i,k)
    enddo
  enddo

end subroutine compute_tmpi

subroutine dp_inverse(nlev, shcol, rho_zt, dz_zt, rdp_zt)

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: dp_inverse_f
#endif

  !intent-ins
  integer,     intent(in) :: nlev, shcol
  ! Air density on thermo grid [kg/m3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)
  ! height thickness centered on thermo grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)

  !intent-out
  real(rtype), intent(out) :: rdp_zt(shcol,nlev)

  !local vars
  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call dp_inverse_f(nlev, shcol, rho_zt, dz_zt, rdp_zt)
     return
  endif
#endif

  do k = 1, nlev
    do i = 1, shcol
      rdp_zt(i,k) = 1._rtype/(ggr*rho_zt(i,k)*dz_zt(i,k))
    enddo
  enddo

end subroutine dp_inverse

function impli_srf_stress_term(shcol, rho_zi_sfc, uw_sfc, &
     vw_sfc, u_wind_sfc, v_wind_sfc) result (ksrf)

  !intent-ins
  integer,     intent(in) :: shcol

  !air density at interfaces [kg/m3]
  real(rtype), intent(in) :: rho_zi_sfc(shcol)
  !vertical zonal momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: uw_sfc(shcol)
  !vertical meridional momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: vw_sfc(shcol)
  !zonal wind [m/s]
  real(rtype), intent(in) :: u_wind_sfc(shcol)
  !meridional wind [m/s]
  real(rtype), intent(in) :: v_wind_sfc(shcol)

  !function return value
  real(rtype) :: ksrf(shcol)

  !local vars
  integer :: i

  real(rtype) :: taux, tauy !stresses (N/m2)
  real(rtype) :: ws         !wind speed (m/s)
  real(rtype) :: rho, tau, uw, vw

  real(rtype), parameter :: wsmin    = 1._rtype    ! Minimum wind speed for ksrfturb computation [ m/s ]
  real(rtype), parameter :: ksrfmin  = 1.e-4_rtype ! Minimum surface drag coefficient  [ kg/s/m^2 ]

  !store surface values of rho in a 1d array

  do i = 1, shcol
     rho          = rho_zi_sfc(i)
     uw           = uw_sfc(i)
     vw           = vw_sfc(i)

     taux         = rho*uw ! stress in N/m2
     tauy         = rho*vw ! stress in N/m2
     ! compute the wind speed
     ws           = max(bfb_sqrt(bfb_square(u_wind_sfc(i)) + bfb_square(v_wind_sfc(i))),wsmin)
     tau          = bfb_sqrt(bfb_square(taux) + bfb_square(tauy))
     ksrf(i)      = max(tau/ws, ksrfmin)
  enddo

  return
end function impli_srf_stress_term

function tke_srf_flux_term(shcol, uw_sfc, vw_sfc) result(wtke_sfc)

  !intent-ins
  integer,     intent(in) :: shcol

  !vertical zonal momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: uw_sfc(shcol)
  !vertical meridional momentum flux at surface [m3/s3]
  real(rtype), intent(in) :: vw_sfc(shcol)

  !function return value
  real(rtype) :: wtke_sfc(shcol)

  !local vars
  integer :: i

  real(rtype) :: ustar, uw, vw

  real(rtype), parameter :: ustarmin = 0.01_rtype  ! Minimum ustar

  do i = 1, shcol
     uw           = uw_sfc(i)
     vw           = vw_sfc(i)
     ustar        = max(bfb_sqrt(bfb_sqrt(bfb_square(uw) + bfb_square(vw))),ustarmin)
     wtke_sfc(i) = bfb_cube(ustar)
  enddo

  return
end function tke_srf_flux_term


subroutine sfc_fluxes(shcol, num_tracer, dtime, rho_zi_sfc, rdp_zt_sfc, wthl_sfc,  &
                      wqw_sfc, wtke_sfc, wtracer_sfc, thetal, qw, tke, wtracer)

  implicit none

  !intent-ins
  integer,     intent(in) :: shcol
  !number of tracers
  integer,     intent(in) :: num_tracer
  !time step [s]
  real(rtype), intent(in) :: dtime
  !air density at interfaces [kg/m3]
  real(rtype), intent(in) :: rho_zi_sfc(shcol)
  !inverse of dp
  real(rtype), intent(in) :: rdp_zt_sfc(shcol)
  !vertical heat flux at surface [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  !vertical moisture flux at surface [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  !vertical tke flux at surface [m3/s3]
  real(rtype), intent(in) :: wtke_sfc(shcol)
  !vertical tracer flux at surface [units vary m/s]
  real(rtype), intent(in) :: wtracer_sfc(shcol,num_tracer)

  !intent-inouts
  !liquid water potential temperature [K]
  real(rtype), intent(inout) :: thetal(shcol)
  !total water mixing ratio [kg/kg]
  real(rtype), intent(inout) :: qw(shcol)
  !turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol)
  !tracers [units vary]
  real(rtype), intent(inout) :: wtracer(shcol,num_tracer)

  !local variables
  integer :: i, p
  real(rtype) :: cmnfac

  ! Apply the surface fluxes explicitly for temperature and moisture
  do i = 1, shcol
     cmnfac       =  dtime * (ggr * rho_zi_sfc(i) * rdp_zt_sfc(i)) !a common factor for the following 3 equations

     thetal(i) = thetal(i) + cmnfac * wthl_sfc(i)
     qw(i)     = qw(i)     + cmnfac * wqw_sfc(i)
     tke(i)    = tke(i)    + cmnfac * wtke_sfc(i)

     ! surface fluxes for tracers
     do p = 1, num_tracer
        wtracer(i,p) = wtracer(i,p) + cmnfac * wtracer_sfc(i,p)
     enddo
  enddo

end subroutine sfc_fluxes


!=======================================================
! SHOC Diagnose the second order moments,
!  main routine

subroutine diag_second_shoc_moments(&
         shcol,nlev,nlevi, &                    ! Input
         thetal,qw,u_wind,v_wind,tke, &         ! Input
         isotropy,tkh,tk, &                     ! Input
         dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
         wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &   ! Input
         do_mf, ae, aw, awthl, awqt, &          ! EDMF Input         
         thl_sec,qw_sec,wthl_sec,wqw_sec, &     ! Output
         qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Output
         w_sec, &                               ! Output
         mf_thlflx, mf_qtflx, &                 ! EDMF Output
         wthl_sec_ed, wthl_sec_mf, &            ! EDMF Output
         wqw_sec_ed, wqw_sec_mf)                ! EDMF Output
         
#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: diag_second_shoc_moments_f
#endif

  ! This is the main routine to compute the second
  !   order moments in SHOC.

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi

  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: qw(shcol,nlev)
  ! zonal wind component [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! meridional wind component [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! return to isotropy timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(in) :: tkh(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)
  ! heights of mid-point grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)

  ! MJC: EDMF inputs
  ! Logical flag: Include MF in fluxes?
  logical,     intent(in) :: do_mf
  ! EDMF environment area [-]
  real(rtype), intent(in) :: ae(shcol,nlevi)
  ! EDMF area-weighted plume mean updraft speed [m/s]
  real(rtype), intent(in) :: aw(shcol,nlevi)
  ! EDMF area-weighted plume temperature transport [Km/s]
  real(rtype), intent(in) :: awthl(shcol,nlevi)
  ! EDMF area_weighted plume moisture transport [kgm/kgs]
  real(rtype), intent(in) :: awqt(shcol,nlevi)
  
! OUTPUT VARIABLES
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(out) :: thl_sec(shcol,nlevi)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(out) :: qw_sec(shcol,nlevi)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(out) :: qwthl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(out) :: wthl_sec(shcol,nlevi)
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(out) :: wqw_sec(shcol,nlevi)
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(out) :: uw_sec(shcol,nlevi)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(out) :: vw_sec(shcol,nlevi)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(out) :: wtke_sec(shcol,nlevi)
  ! second order vertical velocity [m2/s2]
  real(rtype), intent(out) :: w_sec(shcol,nlev)
  
  ! MJC: EDMF outputs
  ! MF temperature turbulent flux [Km/s]
  real(rtype), intent(out) :: mf_thlflx(shcol,nlevi)
  ! MF moisture turbulent flux [kgm/kgs]
  real(rtype), intent(out) :: mf_qtflx(shcol,nlevi)
  
  real(rtype), intent(out) :: wthl_sec_ed(shcol,nlevi)
  real(rtype), intent(out) :: wthl_sec_mf(shcol,nlevi)
  real(rtype), intent(out) :: wqw_sec_ed(shcol,nlevi)
  real(rtype), intent(out) :: wqw_sec_mf(shcol,nlevi)

! LOCAL VARIABLES
  real(rtype) :: wstar(shcol)
  real(rtype) :: ustar2(shcol)

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
     call  diag_second_shoc_moments_f(shcol,nlev,nlevi,  &
              thetal,qw,u_wind,v_wind,tke, isotropy,tkh,tk, dz_zi,zt_grid,zi_grid,shoc_mix, &
              wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,thl_sec,qw_sec,wthl_sec,wqw_sec,&
              qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec)
      return
   endif
#endif

  ! Calculate surface properties needed for lower
  !  boundary conditions
  call diag_second_moments_srf(&
     shcol,       &                         ! Input
     wthl_sfc, uw_sfc, vw_sfc, &            ! Input
     ustar2,wstar)                          ! Output

  ! Diagnose the second order moments flux,
  !  for the lower boundary
  call diag_second_moments_lbycond(&
     shcol,&                                             ! Input
     wthl_sfc, wqw_sfc, uw_sfc, vw_sfc,&                 ! Input
     ustar2,wstar,&                                      ! Input
     wthl_sec(:shcol,nlevi),wqw_sec(:shcol,nlevi),&      ! Output
     uw_sec(:shcol,nlevi), vw_sec(:shcol,nlevi),&        ! Output
     wtke_sec(:shcol,nlevi), thl_sec(:shcol,nlevi),&     ! Output
     qw_sec(:shcol,nlevi), qwthl_sec(:shcol,nlevi),&     ! Output
     wthl_sec_ed(:shcol,nlevi),wqw_sec_ed(:shcol,nlevi)) ! EDMF Output

  ! Diagnose the second order moments,
  !  for points away from boundaries.  this is
  !  the main computation for the second moments
  call diag_second_moments(&
     shcol, nlev, nlevi, &                  ! Input
     thetal, qw, u_wind, v_wind, tke, &     ! Input
     isotropy, tkh, tk,&                    ! Input
     dz_zi, zt_grid, zi_grid, shoc_mix, &   ! Input
     do_mf, ae, aw, awthl, awqt, &          ! EDMF Input
     thl_sec, qw_sec,wthl_sec,wqw_sec,&     ! Input/Output
     qwthl_sec, uw_sec, vw_sec, wtke_sec, & ! Input/Output
     w_sec, &                               ! Output
     mf_thlflx, mf_qtflx, &                 ! EDMF Output
     wthl_sec_ed, wthl_sec_mf, & 			      ! EDMF Output
     wqw_sec_ed, wqw_sec_mf) 				        ! EDMF Output

  ! Diagnose the second order moments,
  !  calculate the upper boundary conditions
  call diag_second_moments_ubycond(&
     shcol,                              &        ! Input
     thl_sec(:shcol,1), qw_sec(:shcol,1),&        ! Output
     wthl_sec(:shcol,1),wqw_sec(:shcol,1),& 	    ! Output
     qwthl_sec(:shcol,1), uw_sec(:shcol,1),&	    ! Output
     vw_sec(:shcol,1), wtke_sec(:shcol,1),&   	  ! Output
     wthl_sec_ed(:shcol,1), wqw_sec_ed(:shcol,1)) ! EDMF Output

  return
end subroutine diag_second_shoc_moments

!==============================================================
! SHOC Diagnose the second order moments,
!  lower boundary conditions

subroutine diag_second_moments_srf(&
         shcol,       &                         ! Input
         wthl_sfc, uw_sfc, vw_sfc, &            ! Input
         ustar2,wstar)                          ! Output

  ! Purpose of this subroutine is to diagnose surface
  !  properties needed for the the lower
  !  boundary condition for the second order moments needed
  !  for the SHOC parameterization.
#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_diag_second_moments_srf_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol

  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)

! OUTPUT VARIABLES
  ! Surface friction velocity [m4/s4]
  real(rtype), intent(out) :: ustar2(shcol)
  ! Surface convective velocity [m/s]
  real(rtype), intent(out) :: wstar(shcol)

! LOCAL VARIABLES
  integer :: i

  ! Constants to parameterize surface variances
  real(rtype), parameter :: z_const = 1.0_rtype

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_diag_second_moments_srf_f(shcol,wthl_sfc, uw_sfc, vw_sfc, &            ! Input
            ustar2,wstar)                          ! Output
      return
   endif
#endif

  ! apply the surface conditions to diagnose turbulent
  !  moments at the surface
  do i=1,shcol

    ! Parameterize thermodyanmics variances via Andre et al. 1978
    ustar2(i) = bfb_sqrt(uw_sfc(i) * uw_sfc(i) + vw_sfc(i) * vw_sfc(i))
    if (wthl_sfc(i) > 0._rtype) then
      wstar(i) = bfb_pow((1._rtype/basetemp * ggr * wthl_sfc(i) * z_const), (1._rtype/3._rtype))
    else
      wstar(i) = 0._rtype
    endif

  enddo ! end i loop (column loop)
  return
end subroutine diag_second_moments_srf

!==============================================================
! SHOC Diagnose the second order moments flux,
!  lower boundary conditions

subroutine diag_second_moments_lbycond(&
         shcol, &                                     ! Input
         wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &         ! Input
         ustar2,wstar, wthl_sec, wqw_sec,&            ! Output
         uw_sec, vw_sec, wtke_sec,&                   ! Output
         thl_sec, qw_sec, qwthl_sec,&                 ! Output
         wthl_sec_ed, wqw_sec_ed)                     ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: diag_second_moments_lbycond_f
#endif

  ! Purpose of this subroutine is to diagnose the lower
  !  boundary condition for the second order moments needed
  !  for the SHOC parameterization.
  ! The thermodymnamic, tracer, and momentum fluxes are set
  !  to the surface fluxes for the host model, while the
  !  thermodynamic variances and covariances are computed
  !  according to that of Andre et al. 1978.

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol

  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! Surface friction velocity squared [m4/s4]
  real(rtype), intent(in) :: ustar2(shcol)
  ! Surface convective velocity scale [m/s]
  real(rtype), intent(in) :: wstar(shcol)

! OUTPUT VARIABLES
  ! vertical flux of heat [K m/s]
  real(rtype), intent(out) :: wthl_sec(shcol)
  real(rtype), intent(out) :: wthl_sec_ed(shcol)
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(out) :: wqw_sec(shcol)
  real(rtype), intent(out) :: wqw_sec_ed(shcol)  
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(out) :: uw_sec(shcol)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(out) :: vw_sec(shcol)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(out) :: wtke_sec(shcol)
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(out) :: thl_sec(shcol)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(out) :: qw_sec(shcol)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(out) :: qwthl_sec(shcol)

! LOCAL VARIABLES
  integer :: i
  real(rtype) :: uf

  ! Constants to parameterize surface variances
  real(rtype), parameter :: a_const = 1.8_rtype
  real(rtype), parameter :: ufmin = 0.01_rtype

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call diag_second_moments_lbycond_f(     &
                shcol,           &                           ! Input
                wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, &         ! Input
                ustar2,wstar,            &                   ! Input
                wthl_sec,wqw_sec,&                           ! Output
                uw_sec, vw_sec, wtke_sec,&                   ! Output
                thl_sec,qw_sec,qwthl_sec)                    ! Output
      return
   endif
#endif

  ! apply the surface conditions to diagnose turbulent
  !  moments at the surface
  do i=1,shcol

    uf = bfb_sqrt(ustar2(i) + 0.3_rtype * wstar(i) * wstar(i))
    uf = max(ufmin,uf)

    ! Diagnose thermodynamics variances and covariances
    thl_sec(i) = 0.4_rtype * a_const * bfb_square(wthl_sfc(i)/uf)
    qw_sec(i) = 0.4_rtype * a_const * bfb_square(wqw_sfc(i)/uf)
    qwthl_sec(i) = 0.2_rtype * a_const * (wthl_sfc(i)/uf) * &
                         (wqw_sfc(i)/uf)

    ! Vertical fluxes of heat and moisture, simply
    !  use the surface fluxes given by host model
    wthl_sec(i) = wthl_sfc(i)
    wqw_sec(i) = wqw_sfc(i)
    uw_sec(i) = uw_sfc(i)
    vw_sec(i) = vw_sfc(i)
    wtke_sec(i) = bfb_cube(max(bfb_sqrt(ustar2(i)),0.01_rtype))

	wthl_sec_ed(i) = wthl_sfc(i)
    wqw_sec_ed(i) = wqw_sfc(i)
    
  enddo ! end i loop (column loop)
  return
end subroutine diag_second_moments_lbycond

subroutine diag_second_moments(&
         shcol,nlev,nlevi, &                    ! Input
         thetal,qw,u_wind,v_wind,tke, &         ! Input
         isotropy,tkh,tk, &                     ! Input
         dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
         do_mf, ae, aw, awthl, awqt, &          ! Input - EDMF         
         thl_sec,qw_sec,wthl_sec,wqw_sec, &     ! Input/Output
         qwthl_sec,uw_sec,vw_sec,wtke_sec, &    ! Input/Output
         w_sec, &                               ! Output
         mf_thlflx, mf_qtflx, &  				        ! Output - EDMF
         wthl_sec_ed,wthl_sec_mf, &   	        ! Output - EDMF
         wqw_sec_ed, wqw_sec_mf)                ! Output - EDMF
         
         
#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: diag_second_moments_f
#endif

  ! Purpose of this subroutine is to diagnose the second
  !  order moments needed for the SHOC parameterization.
  !  Namely these are variances of thetal, qw, and vertical
  !  velocity.  In addition the vertical fluxes of thetal, qw,
  !  u, v, TKE, and tracers are computed here as well as the
  !  correlation of qw and thetal.

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi

  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: qw(shcol,nlev)
  ! zonal wind component [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! meridional wind component [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! return to isotropy timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(in) :: tkh(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)
  ! heights of mid-point grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  
  ! MJC: EDMF inputs
  logical    , intent(in) :: do_mf
  real(rtype), intent(in) :: ae(shcol,nlevi)
  real(rtype), intent(in) :: aw(shcol,nlevi)
  real(rtype), intent(in) :: awthl(shcol,nlevi)
  real(rtype), intent(in) :: awqt(shcol,nlevi)
  
! INPUT/OUTPUT VARIABLES
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(inout) :: thl_sec(shcol,nlevi)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(inout) :: qw_sec(shcol,nlevi)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(inout) :: qwthl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(inout) :: wthl_sec(shcol,nlevi)
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(inout) :: wqw_sec(shcol,nlevi)
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(inout) :: uw_sec(shcol,nlevi)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(inout) :: vw_sec(shcol,nlevi)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(inout) :: wtke_sec(shcol,nlevi)

! OUTPUT VARIABLES
  ! second order vertical velocity [m2/s2]
  real(rtype), intent(out) :: w_sec(shcol,nlev)

  ! MJC: EDMF Output
  real(rtype), intent(out) :: wthl_sec_ed(shcol,nlevi)
  real(rtype), intent(out) :: wthl_sec_mf(shcol,nlevi)
  real(rtype), intent(out) :: wqw_sec_ed(shcol,nlevi)
  real(rtype), intent(out) :: wqw_sec_mf(shcol,nlevi)  
  ! vertical flux of heat from mass flux plumes [K m/s]
  real(rtype), intent(out) :: mf_thlflx(shcol,nlevi)
  ! vertical flux of moisture from mass flux plumes [kg/kg m/s]
  real(rtype), intent(out) :: mf_qtflx(shcol,nlevi)
  
  ! LOCAL VARIABLES
  real(rtype) :: isotropy_zi(shcol,nlevi)
  real(rtype) :: tkh_zi(shcol,nlevi)
  real(rtype) :: tk_zi(shcol,nlevi)
  
  ! MJC: Extra local variables
  real(rtype) :: thl_zi(shcol,nlevi)
  real(rtype) :: qw_zi(shcol,nlevi)
  real(rtype) :: wtke_sec_ed(shcol,nlevi)
  real(rtype) :: wtke_sec_mf(shcol,nlevi)
  real(rtype) :: uw_sec_ed(shcol,nlevi)
  real(rtype) :: uw_sec_mf(shcol,nlevi)
  real(rtype) :: vw_sec_ed(shcol,nlevi)
  real(rtype) :: vw_sec_mf(shcol,nlevi)
  
  ! Determines if total fluxes (ED+MF) are computed or not
  logical :: do_total_fluxes
  !Dummy variable
  real(rtype) :: aw_dummy(shcol,nlevi)

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
     call diag_second_moments_f(shcol,nlev,nlevi,thetal,qw,u_wind,v_wind,tke, &         ! Input
                                isotropy,tkh,tk,dz_zi,zt_grid,zi_grid,shoc_mix, &      ! Input
                                thl_sec,qw_sec,wthl_sec,wqw_sec,qwthl_sec,uw_sec,vw_sec,wtke_sec, &    ! Input/Output
                                w_sec)
      return
   endif
#endif

  ! Interpolate some variables from the midpoint grid to the interface grid
  call linear_interp(zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,shcol,0._rtype)
  ! MJC: 
  call linear_interp(zt_grid,zi_grid,thetal,thl_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,qw,qw_zi,nlev,nlevi,shcol,0._rtype)
  
  ! Vertical velocity variance is assumed to be propotional
  !  to the TKE
  w_sec = w2tune*(2._rtype/3._rtype)*tke

  ! Calculate the temperature variance
  call calc_shoc_varorcovar(&
         shcol,nlev,nlevi,thl2tune,&              ! Input
         isotropy_zi,tkh_zi,dz_zi,thetal,thetal,& ! Input
         thl_sec)                                 ! Input/Output

  ! Calculate the moisture variance
  call calc_shoc_varorcovar(&
         shcol,nlev,nlevi,qw2tune,&               ! Input
         isotropy_zi,tkh_zi,dz_zi,qw,qw,&         ! Input
         qw_sec)                                  ! Input/Output

  ! Calculate the temperature and moisture covariance
  call calc_shoc_varorcovar(&
         shcol,nlev,nlevi,qwthl2tune,&            ! Input
         isotropy_zi,tkh_zi,dz_zi,thetal,qw,&     ! Input
         qwthl_sec)                               ! Input/Output

  ! Calculate vertical flux for heat
  ! MJC: Flag for MF fluxes
  do_total_fluxes = do_mf
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,thetal,&   ! Input
         zt_grid,zi_grid,ae,aw,awthl,&            ! Input
         do_total_fluxes,&                        ! Input - Logical variable         
         wthl_sec,wthl_sec_ed,wthl_sec_mf)        ! Input/Output

  ! MJC: MF flux component
  call calc_mf_vertflux(shcol,nlev,nlevi,aw,awthl,thetal,thl_zi,mf_thlflx)

  ! Calculate vertical flux for moisture
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,qw,&       ! Input
         zt_grid,zi_grid,ae,aw,awqt,&             ! Input
         do_total_fluxes,&                        ! Input - Logical variable         
         wqw_sec,wqw_sec_ed,wqw_sec_mf)           ! Input/Output
  
  ! MJC: MF flux component       
  call calc_mf_vertflux(shcol,nlev,nlevi,aw,awqt,qw,qw_zi,mf_qtflx)
      
  ! Calculate vertical flux for TKE
  do_total_fluxes = .false.
  aw_dummy = 0._rtype  
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,tke,&     ! Input
        zt_grid,zi_grid,ae,aw,aw_dummy,&         ! Input
         do_total_fluxes,&                       ! Input - Logical variable
         wtke_sec,wtke_sec_ed,wtke_sec_mf)       ! Input/Output

  ! Calculate vertical flux for momentum (zonal wind)
  do_total_fluxes = .false.
  aw_dummy = 0._rtype  
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tk_zi,dz_zi,u_wind,&    ! Input
         zt_grid,zi_grid,ae,aw,aw_dummy,&         ! Input
         do_total_fluxes,&                        ! Input - Logical variable
         uw_sec,uw_sec_ed,uw_sec_mf)              ! Input/Output
         
  ! Calculate vertical flux for momentum (meridional wind)
  do_total_fluxes = .false.
  aw_dummy = 0._rtype  
  call calc_shoc_vertflux(&
         shcol,nlev,nlevi,tk_zi,dz_zi,v_wind,&    ! Input
         zt_grid,zi_grid,ae,aw,aw_dummy,&         ! Input
         do_total_fluxes,&                        ! Input - Logical variable
         vw_sec,vw_sec_ed,vw_sec_mf)              ! Input/Output
         
  return
end subroutine diag_second_moments

subroutine calc_shoc_varorcovar(&
         shcol,nlev,nlevi,tunefac,&                ! Input
         isotropy_zi,tkh_zi,dz_zi,invar1,invar2,&  ! Input
         varorcovar)                               ! Input/Output

  ! Compute either the variance or covariance
  !  (depending on if invar1 is the same as invar2)
  !  for a given set of inputs

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: calc_shoc_varorcovar_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! tuning factor for (co)variance []
  real(rtype), intent(in) :: tunefac
  ! Return to isotopic timescale [s]
  real(rtype), intent(in) :: isotropy_zi(shcol,nlevi)
  ! Eddy diffusivity for heat [ms-2]
  real(rtype), intent(in) :: tkh_zi(shcol,nlevi)
  ! delta z centerend on zi grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Input variable 1 [units vary]
  real(rtype), intent(in) :: invar1(shcol,nlev)
  ! Input variable 2 [units vary]
  real(rtype), intent(in) :: invar2(shcol,nlev)

! INPUT/OUTPUT VARIABLES
  ! variance or covariance [units vary]
  real(rtype), intent(inout) :: varorcovar(shcol,nlevi)

! LOCAL VARIABLES
  integer :: i, k, kt
  real(rtype) :: sm, grid_dz2

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call calc_shoc_varorcovar_f(shcol,nlev,nlevi,tunefac,isotropy_zi,tkh_zi,dz_zi,invar1,invar2,&  ! Input
           varorcovar)                              ! Input/Output)
      return
   endif
#endif

  do k=2,nlev

    kt=k-1 ! define upper grid point indicee
    do i=1,shcol

      grid_dz2=bfb_square(1._rtype/dz_zi(i,k)) ! vertical grid diff squared
      sm=isotropy_zi(i,k)*tkh_zi(i,k) ! coefficient for variances

      ! Compute the variance or covariance
      varorcovar(i,k)=tunefac*sm*grid_dz2*(invar1(i,kt)-invar1(i,k))*&
        (invar2(i,kt)-invar2(i,k))

    enddo
  enddo

  return
end subroutine calc_shoc_varorcovar

subroutine calc_shoc_vertflux(&
         shcol,nlev,nlevi,tkh_zi,dz_zi,invar,&     ! Input
         zt_grid,zi_grid,mf_ae,mf_aw,mf_aw_invar,& ! EDMF Input
         do_total_fluxes,&                         ! EDMF Input - Logical flag         
         vertflux,vertflux_ed,vertflux_mf)         ! Input/Output

  ! Compute either the vertical flux via
  !  downgradient diffusion for a given set of
  !  input variables

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: calc_shoc_vertflux_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! Eddy diffusivity for heat [ms-2]
  real(rtype), intent(in) :: tkh_zi(shcol,nlevi)
  ! delta z centerend on zi grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! Heights of mid-point grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! Heights of interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! Input variable [units vary]
  real(rtype), intent(in) :: invar(shcol,nlev)
  
  ! MJC: Input from the MF subroutine
  real(rtype), intent(in) :: mf_ae(shcol,nlevi)       ! EDMF: Environmental area
  real(rtype), intent(in) :: mf_aw(shcol,nlevi)       ! EDMF: sum(a_i * w_i) [m/s]
  real(rtype), intent(in) :: mf_aw_invar(shcol,nlevi) ! EDMF: sum(a_i * w_i * var_in) [? m/s]

  ! MJC: EDMF logical variable
  logical, intent(in) :: do_total_fluxes
  
! INPUT/OUTPUT VARIABLES
  real(rtype), intent(out) :: vertflux(shcol,nlevi)
  real(rtype), intent(out) :: vertflux_ed(shcol,nlevi)
  real(rtype), intent(out) :: vertflux_mf(shcol,nlevi)
  
! LOCAL VARIABLES
  integer :: i, k, kt
  real(rtype) :: grid_dz

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call calc_shoc_vertflux_f(shcol,nlev,nlevi,tkh_zi,dz_zi,invar,&  ! Input
           vertflux)                                                   ! Input/Output
      return
   endif
#endif

  if (do_total_fluxes) then
    do k=2,nlev

      kt=k-1 ! define upper grid point indice
      do i=1,shcol

        grid_dz = 1._rtype/dz_zi(i,k) ! vertical grid diff squared
            
        !! Centered
        vertflux(i,k)=(-1._rtype*tkh_zi(i,k)*grid_dz*(invar(i,kt)-invar(i,k))) + &
          mf_aw_invar(i,k) - mf_aw(i,k)*0.5_rtype*(invar(i,kt)+invar(i,k)) 
          
        vertflux_ed(i,k)= -1._rtype*tkh_zi(i,k)*grid_dz*(invar(i,kt)-invar(i,k))
        
        vertflux_mf(i,k)= mf_aw_invar(i,k) - mf_aw(i,k)*0.5_rtype*(invar(i,kt)+invar(i,k)) 
        
        !! MF upwind
        !vertflux(i,k)=(-1._rtype*mf_ae(i,k)*tkh_zi(i,k)*grid_dz*(invar(i,kt)-invar(i,k))) + &
        !  mf_aw_invar(i,k) - mf_aw(i,k)*invar(i,k)   
        !! MF downwind
        !vertflux(i,k)=(-1._rtype*mf_ae(i,k)*tkh_zi(i,k)*grid_dz*(invar(i,kt)-invar(i,k))) + &
        !  mf_aw_invar(i,k) - mf_aw(i,k)*invar(i,kt)
      enddo
    enddo
  else

    do k=2,nlev
 
      kt=k-1 ! define upper grid point indice
      do i=1,shcol

        grid_dz=1._rtype/dz_zi(i,k) ! vertical grid diff squared

        ! Compute the vertical flux via downgradient diffusion
        vertflux(i,k)=-1._rtype*tkh_zi(i,k)*grid_dz*&
          (invar(i,kt)-invar(i,k))

        vertflux_ed(i,k)= -1._rtype*tkh_zi(i,k)*grid_dz*(invar(i,kt)-invar(i,k))
        vertflux_mf(i,k)= 0._rtype  
        
      enddo
    enddo
    
  endif

  return
end subroutine calc_shoc_vertflux

subroutine diag_second_moments_ubycond(&
         shcol, &                               ! Input
         thl_sec, qw_sec,&                      ! Output
         wthl_sec,wqw_sec,&                     ! Output
         qwthl_sec, uw_sec, vw_sec, wtke_sec,&  ! Output
         wthl_sec_ed, wqw_sec_ed)               ! Output

  ! Purpose of this subroutine is to diagnose the upper
  !  boundary condition for the second order moments
  !  needed for the SHOC parameterization.  Currently
  !  set all to zero.

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_diag_second_moments_ubycond_f
#endif
  implicit none

  ! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol

  ! OUTPUT VARIABLES
  ! second order liquid wat. potential temp. [K^2]
  real(rtype), intent(out) :: thl_sec(shcol)
  ! second order total water mixing rat. [kg^2/kg^2]
  real(rtype), intent(out) :: qw_sec(shcol)
  ! covariance of temp and moisture [K kg/kg]
  real(rtype), intent(out) :: qwthl_sec(shcol)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(out) :: wthl_sec(shcol)
  real(rtype), intent(out) :: wthl_sec_ed(shcol)  
  ! vertical flux of total water [kg/kg m/s]
  real(rtype), intent(out) :: wqw_sec(shcol)
  real(rtype), intent(out) :: wqw_sec_ed(shcol)  
  ! vertical flux of zonal wind [m2/s2]
  real(rtype), intent(out) :: uw_sec(shcol)
  ! vertical flux of meridional wind [m2/s2]
  real(rtype), intent(out) :: vw_sec(shcol)
  ! vertical flux of tke [m3/s3]
  real(rtype), intent(out) :: wtke_sec(shcol)

  ! LOCAL VARIABLES
  integer :: i

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
       call shoc_diag_second_moments_ubycond_f(&
                             shcol, &                    ! Input
                             thl_sec, qw_sec,&                      ! Output
                             wthl_sec,wqw_sec,&                     ! Output
                             qwthl_sec, uw_sec, vw_sec, wtke_sec)   ! Output
      return
   endif
#endif

  ! apply the upper boundary condition
  do i=1,shcol
    wthl_sec(i) = 0._rtype
    wqw_sec(i) = 0._rtype
    uw_sec(i) = 0._rtype
    vw_sec(i) = 0._rtype
    wtke_sec(i) = 0._rtype

    thl_sec(i) = 0._rtype
    qw_sec(i) = 0._rtype
    qwthl_sec(i) = 0._rtype
    
    wthl_sec_ed(i) = 0._rtype
    wqw_sec_ed(i) = 0._rtype    
  enddo ! end i loop (column loop)
  return
end subroutine diag_second_moments_ubycond

!==============================================================
! SHOC Diagnose the third order moment of vertical velocity

subroutine diag_third_shoc_moments(&
          shcol,nlev,nlevi, &                 ! Input
          w_sec, thl_sec, &                   ! Input
          wthl_sec, isotropy, brunt,&         ! Input
          thetal,tke,&                        ! Input
          dz_zt, dz_zi, zt_grid, zi_grid,&    ! Input
          w3)                                 ! Output

  ! Purpose of this subroutine is to diagnose the third
  !  order moment of the vertical velocity, needed
  !  for the skewness calculation in the PDF.
  !  This calculation follows that of Canuto et al. (2001)

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: diag_third_shoc_moments_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi

  ! second order vertical velocity [m2/s2]
  real(rtype), intent(in) :: w_sec(shcol,nlev)
  ! second order liquid wat. potential temperature [K^2]
  real(rtype), intent(in) :: thl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
  ! return to isotropy timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! brunt vaisallia frequency [s]
  real(rtype), intent(in) :: brunt(shcol,nlev)
  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! thickness centered on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! heights of thermodynamics points [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights of interface points [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)

! OUTPUT VARIABLES
  ! third moment of vertical velocity
  real(rtype), intent(out) :: w3(shcol,nlevi)

! LOCAL VARIABLES
  real(rtype) :: w_sec_zi(shcol,nlevi)    ! second order vertical velocity
  real(rtype) :: isotropy_zi(shcol,nlevi)
  real(rtype) :: brunt_zi(shcol,nlevi)
  real(rtype) :: thetal_zi(shcol,nlevi)

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
      call diag_third_shoc_moments_f(&
        shcol,nlev,nlevi, &                 ! Input
        w_sec, thl_sec, &                   ! Input
        wthl_sec, isotropy, brunt,&         ! Input
        thetal,tke,&                        ! Input
        dz_zt, dz_zi, zt_grid, zi_grid,&    ! Input
        w3)                                 ! Output
     return
  endif
#endif

  ! Interpolate variables onto the interface levels
  call linear_interp(zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,shcol,0._rtype)
  call linear_interp(zt_grid,zi_grid,brunt,brunt_zi,nlev,nlevi,shcol,largeneg)
  call linear_interp(zt_grid,zi_grid,w_sec,w_sec_zi,nlev,nlevi,shcol,(2._rtype/3._rtype)*mintke)
  call linear_interp(zt_grid,zi_grid,thetal,thetal_zi,nlev,nlevi,shcol,0._rtype)

  !Diagnose the third moment of the vertical-velocity
  call compute_diag_third_shoc_moment(&
          shcol,nlev,nlevi, w_sec,thl_sec, & ! Input
          wthl_sec, tke, dz_zt, dz_zi,&      ! Input
          isotropy_zi, brunt_zi,w_sec_zi,&   ! Input
          thetal_zi,&                        ! Input
          w3)                                ! Output

  ! perform clipping to prevent unrealistically large values from occuring
  call clipping_diag_third_shoc_moments(&
          nlevi,shcol,w_sec_zi,&    !Input
          w3)                       !Input/Output

  return

end subroutine diag_third_shoc_moments

subroutine compute_diag_third_shoc_moment(&
          shcol,nlev,nlevi, w_sec, thl_sec,& ! Input
          wthl_sec, tke, dz_zt, dz_zi,&      ! Input
          isotropy_zi, brunt_zi,w_sec_zi,&   ! Input
          thetal_zi,&                        ! Input
          w3)                                ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_diag_third_shoc_moment_f
#endif

  implicit none
! INPUT VARIABLES
  ! number of SHOC columns
  integer, intent(in) :: shcol
  ! number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi
  ! second order vertical velocity [m2/s2]
  real(rtype), intent(in) :: w_sec(shcol,nlev)
  ! second order liquid wat. potential temperature [K^2]
  real(rtype), intent(in) :: thl_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! thickness centered on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! thickness centered on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)

  !Interpolated varaibles
  real(rtype), intent(in) :: isotropy_zi(shcol,nlevi)
  real(rtype), intent(in) :: brunt_zi(shcol,nlevi)
  real(rtype), intent(in) :: w_sec_zi(shcol,nlevi)
  real(rtype), intent(in) :: thetal_zi(shcol,nlevi)
  ! third moment of vertical velocity
  real(rtype), intent(out) :: w3(shcol,nlevi)
! LOCAL VARIABLES
  integer i, k, kb, kc
  real(rtype) :: omega0, omega1, omega2
  real(rtype) :: X0, Y0, X1, Y1, AA0, AA1
  real(rtype) :: iso, thedz, thedz2
  real(rtype) :: isosqrd
  real(rtype) :: buoy_sgs2, bet2
  real(rtype) :: f0, f1, f2, f3, f4, f5

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call compute_diag_third_shoc_moment_f(shcol,nlev,nlevi,w_sec,thl_sec, &   ! Input
                                          wthl_sec, tke, dz_zt, dz_zi, &      ! Input
                                          isotropy_zi,brunt_zi,w_sec_zi, &    ! Input
                                          thetal_zi, &                        ! Input
                                          w3)                                 ! Output
    return
  endif
#endif

  ! set lower condition
  w3(:,nlevi) = 0._rtype

  do k=2,nlev

     kb=k+1
     kc=k-1
     do i=1,shcol

        !Compute inputs for computing f0 to f5 terms
        call fterms_input_for_diag_third_shoc_moment(&
             dz_zi(i,k), dz_zt(i,k), dz_zt(i,kc), &               ! Input
             isotropy_zi(i,k), brunt_zi(i,k), thetal_zi(i,k), &   ! Input
             thedz, thedz2, iso, isosqrd, buoy_sgs2, bet2)        ! Output

        !Compute f0 to f5 terms
        call f0_to_f5_diag_third_shoc_moment(&
             thedz, thedz2, bet2, iso, isosqrd, &                 ! Input
             wthl_sec (i,k), wthl_sec(i,kc), wthl_sec(i,kb), &    ! Input
             thl_sec(i,kc), thl_sec(i,kb), &                      ! Input
             w_sec(i,k), w_sec(i,kc), w_sec_zi(i,k), &            ! Input
             tke(i,k), tke(i,kc), &                               ! Input
             f0, f1, f2, f3, f4, f5)                              ! Output

        !Compute the omega terms
        call omega_terms_diag_third_shoc_moment(&
             buoy_sgs2, f3, f4, &       ! Input
             omega0, omega1, omega2)    ! Output

        !Compute the X0, Y0, X1, Y1 terms
        call x_y_terms_diag_third_shoc_moment(&
             buoy_sgs2, f0, f1, f2, &   ! Input
             x0, y0, x1, y1)            ! Output

        !Compute the AA0, AA1 terms
        call aa_terms_diag_third_shoc_moment(&
             omega0, omega1, omega2, &  ! Input
             x0, x1, y0, y1, &          ! Input
             aa0, aa1)                  ! Output

        !Finally, we have the third moment of w
        w3(i,k) = w3_diag_third_shoc_moment(aa0, aa1, x0, x1, f5)

     enddo  ! end i loop (column loop)
  enddo  ! end k loop (vertical loop)

  ! set upper condition
  w3(:,1) = 0._rtype

end subroutine compute_diag_third_shoc_moment

subroutine fterms_input_for_diag_third_shoc_moment(&
     dz_zi, dz_zt, dz_zt_kc, &                      ! Input
     isotropy_zi, brunt_zi, thetal_zi, &            ! Input
     thedz, thedz2, iso, isosqrd, buoy_sgs2, bet2)  ! Output

  !Compute inputs for computing f0 to f5 terms

  implicit none

  !intent-ins
  real(rtype), intent(in) :: dz_zi, dz_zt, dz_zt_kc
  real(rtype), intent(in) :: isotropy_zi, brunt_zi, thetal_zi

  !intent-outs
  real(rtype), intent(out) :: thedz, thedz2, iso, isosqrd
  real(rtype), intent(out) :: buoy_sgs2, bet2

  thedz  = 1._rtype/dz_zi
  thedz2 = 1._rtype/(dz_zt+dz_zt_kc)

  iso       = isotropy_zi
  isosqrd   = bfb_square(iso)
  buoy_sgs2 = isosqrd*brunt_zi
  bet2      = ggr/thetal_zi

  return
end subroutine fterms_input_for_diag_third_shoc_moment

subroutine f0_to_f5_diag_third_shoc_moment(&
     thedz, thedz2, bet2, iso, isosqrd, &    ! Input
     wthl_sec, wthl_sec_kc, wthl_sec_kb, &   ! Input
     thl_sec_kc, thl_sec_kb, &               ! Input
     w_sec, w_sec_kc,w_sec_zi, &             ! Input
     tke, tke_kc, &                          ! Input
     f0, f1, f2, f3, f4, f5)                 ! Output

  !Compute f0 to f5 terms

  implicit none

  !intent-ins
  real(rtype), intent(in) :: thedz, thedz2, bet2, iso, isosqrd
  real(rtype), intent(in) :: wthl_sec, wthl_sec_kc, wthl_sec_kb
  real(rtype), intent(in) :: thl_sec_kc, thl_sec_kb
  real(rtype), intent(in) :: w_sec, w_sec_kc, w_sec_zi, tke, tke_kc

  !intent-out
  real(rtype), intent(out) :: f0, f1, f2, f3, f4, f5

  !local variables
  real(rtype) :: thl_sec_diff, wthl_sec_diff, wsec_diff, tke_diff

  !Some common factors
  thl_sec_diff  = thl_sec_kc  - thl_sec_kb
  wthl_sec_diff = wthl_sec_kc - wthl_sec_kb
  wsec_diff     = w_sec_kc    - w_sec
  tke_diff      = tke_kc      - tke

  f0 = thedz2 * bfb_cube(bet2) * bfb_quad(iso) * wthl_sec * &
       thl_sec_diff

  f1 = thedz2 * bfb_square(bet2) * bfb_cube(iso) * (wthl_sec * &
       wthl_sec_diff + 0.5_rtype * &
       w_sec_zi*thl_sec_diff)

  f2 = thedz * bet2 * isosqrd * wthl_sec * &
       wsec_diff+ 2._rtype * thedz2 * bet2 * &
       isosqrd * w_sec_zi * wthl_sec_diff

  f3 = thedz2 * bet2 * isosqrd * w_sec_zi * &
       wthl_sec_diff + thedz * &
       bet2 * isosqrd * (wthl_sec * tke_diff)

  f4 = thedz * iso * w_sec_zi * (wsec_diff + &
       tke_diff)

  f5 = thedz * iso * w_sec_zi * wsec_diff

  return
end subroutine f0_to_f5_diag_third_shoc_moment

subroutine omega_terms_diag_third_shoc_moment(&
           buoy_sgs2, f3, f4, &    ! Input
           omega0, omega1, omega2) ! Output

  implicit none

  !Compute the omega terms

  !initent-ins
  real(rtype), intent(in) :: buoy_sgs2, f3, f4

  !intent-out
  real(rtype), intent(out) :: omega0, omega1, omega2

  real(rtype) :: a4, a5

  a4=2.4_rtype/(3._rtype*c_diag_3rd_mom+5._rtype)
  a5=0.6_rtype/(c_diag_3rd_mom*(3._rtype+5._rtype*c_diag_3rd_mom))

  omega0 = a4 / (1._rtype - a5 * buoy_sgs2)
  omega1 = omega0/(2._rtype * c_diag_3rd_mom)
  omega2 = omega1 * f3 + (5._rtype/4._rtype) * omega0 * f4

  return
end subroutine omega_terms_diag_third_shoc_moment

subroutine x_y_terms_diag_third_shoc_moment(&
           buoy_sgs2, f0, f1, f2,&  ! Input
           x0, y0, x1, y1)          ! Output

  implicit none

  !Compute the X0, Y0, X1, Y1 terms

  !intent-ins
  real(rtype), intent(in) :: buoy_sgs2, f0, f1, f2

  !intent-outs
  real(rtype), intent(out) :: x0, y0, x1, y1

! local variables
  real(rtype) :: a0, a1, a2, a3

  a0=(0.52_rtype*(1._rtype/bfb_square(c_diag_3rd_mom)))/(c_diag_3rd_mom-2._rtype)
  a1=0.87_rtype/bfb_square(c_diag_3rd_mom)
  a2=0.5_rtype/c_diag_3rd_mom
  a3=0.6_rtype/(c_diag_3rd_mom*(c_diag_3rd_mom-2._rtype))

  x0 = (a2 * buoy_sgs2 * (1._rtype - a3 * buoy_sgs2)) / &
       (1._rtype - (a1 + a3) * buoy_sgs2)
  y0 = (2._rtype * a2 * buoy_sgs2 * x0) / (1._rtype - a3 * buoy_sgs2)
  x1 = (a0 * f0 + a1 * f1 + a2 * (1._rtype - a3 * buoy_sgs2) * f2) / &
       (1._rtype - (a1 + a3) * buoy_sgs2)
  y1 = (2._rtype * a2 * (buoy_sgs2 * x1 + (a0/a1) * f0 + f1)) / &
       (1._rtype - a3* buoy_sgs2)

  return
end subroutine x_y_terms_diag_third_shoc_moment

subroutine aa_terms_diag_third_shoc_moment(&
           omega0, omega1, omega2, & ! Input
           x0, x1, y0, y1, &         ! Input
           aa0, aa1)                 ! Output

  implicit none

  !Compute the AA0, AA1 terms

  !intent-ins
  real(rtype), intent(in) :: omega0, omega1, omega2, x0, x1, y0, y1

  !intent-outs
  real(rtype), intent(out) :: aa0, aa1

  aa0 = omega0 * x0 + omega1 * y0
  aa1 = omega0 * x1 + omega1 * y1 + omega2

  return
end subroutine aa_terms_diag_third_shoc_moment

pure function w3_diag_third_shoc_moment(aa0, aa1, x0, x1, f5) result(w3)

  implicit none

  !Compute third moment of w

  !intent-ins
  real(rtype), intent(in) :: aa0, aa1, x0, x1, f5

  !return type
  real(rtype) :: w3

  w3 = (aa1-1.2_rtype*x1-1.5_rtype*f5)/(c_diag_3rd_mom-1.2_rtype*x0+aa0)

  return
end function w3_diag_third_shoc_moment

subroutine clipping_diag_third_shoc_moments(&
           nlevi,shcol,w_sec_zi,& ! Input
           w3)                    ! Input/Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: clipping_diag_third_shoc_moments_f
#endif

  ! perform clipping to prevent unrealistically large values from occuring

  implicit none

  integer, intent(in) :: nlevi
  integer, intent(in) :: shcol

  real(rtype), intent(in) :: w_sec_zi(shcol,nlevi)
  real(rtype), intent(inout) :: w3(shcol,nlevi)

  real(rtype) :: tsign
  real(rtype) :: cond
  real(rtype) :: theterm

  real(rtype), parameter :: w3clipdef = 0.02_rtype

  integer k, i

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call clipping_diag_third_shoc_moments_f(nlevi,shcol,w_sec_zi, & ! Input
                                              w3)                     ! Input/Output
      return
   endif
#endif

  do k=1, nlevi
    do i=1, shcol

      tsign = 1._rtype
      theterm = w_sec_zi(i,k)
      cond = w3clip * bfb_sqrt(2._rtype * bfb_cube(theterm))
      if (w3(i,k) .lt. 0) tsign = -1._rtype
      if (tsign * w3(i,k) .gt. cond) w3(i,k) = w3clipdef

    enddo !end i loop (column loop)
  enddo ! end k loop (vertical loop)

end subroutine clipping_diag_third_shoc_moments

!==============================================================
! Assumed PDF closure for the SHOC scheme

subroutine shoc_assumed_pdf(&
         shcol,nlev,nlevi, &                ! Input
         thetal,qw,w_field,thl_sec,qw_sec,& ! Input
         wthl_sec,w_sec, &                  ! Input
         wqw_sec,qwthl_sec,w3,pres, &       ! Input
         mf_qlflx_zt, &
         zt_grid,zi_grid,&                  ! Input
         shoc_cldfrac,shoc_ql,&             ! Output
         wqls,wthv_sec,shoc_ql2,&           ! Output
         a1_out,C1_out,C2_out,ql1_out, ql2_out)              ! Output

  ! Purpose of this subroutine is calculate the
  !  double Gaussian PDF of SHOC, which is the centerpiece
  !  of the scheme.  The main outputs are the SGS cloud
  !  fraction and liquid water amount, in addition to the
  !  SGS buoyancy flux which is needed to close the SGS
  !  TKE equation.  This code follows the appendix of
  !  Larson et al. (2002) for Analytic Double Gaussian 1

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: shoc_assumed_pdf_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of midpoint layers
  integer, intent(in) :: nlev
  ! number of interface layers
  integer, intent(in) :: nlevi

  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thetal(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: qw(shcol,nlev)
  ! thetal variance [K^2]
  real(rtype), intent(in) :: thl_sec(shcol,nlevi)
  ! qw variance [kg/kg^2]
  real(rtype), intent(in) :: qw_sec(shcol,nlevi)
  ! vertical flux of heat [K m/s]
  real(rtype), intent(in) :: wthl_sec(shcol,nlevi)
  ! vertical velocity variance [m2/s2]
  real(rtype), intent(in) :: w_sec(shcol,nlev)
  ! vertical flux of moisture [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sec(shcol,nlevi)
  ! qw and thetal correlation [K kg/kg]
  real(rtype), intent(in) :: qwthl_sec(shcol,nlevi)
  ! third moment vertical velocity [m^3/s^3]
  real(rtype), intent(in) :: w3(shcol,nlevi)
  ! large scale vertical velocity [m/s]
  real(rtype), intent(in) :: w_field(shcol,nlev)
  ! pressure [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)

  real(rtype), intent(in) :: mf_qlflx_zt(shcol,nlev)


! OUTPUT VARIABLES
  ! SGS cloud fraction [-]
  real(rtype), intent(out) :: shoc_cldfrac(shcol,nlev)
  ! SGS liquid water mixing ratio [kg/kg]
  real(rtype), intent(out) :: shoc_ql(shcol,nlev)
  ! SGS buoyancy flux [K m/s]
  real(rtype), intent(out) :: wthv_sec(shcol,nlev)
  ! SGS liquid water flux [kg/kg m/s]
  real(rtype), intent(out) :: wqls(shcol,nlev)
  ! SGS liquid water mixing ratio variance [kg/kg]
  real(rtype), intent(out) :: shoc_ql2(shcol,nlev)
  ! MJC [11/03/24]: Adding variables to the output file for diagnostic purposes
  real(rtype), intent(out) :: a1_out(shcol,nlev)
  real(rtype), intent(out) :: C1_out(shcol,nlev)
  real(rtype), intent(out) :: C2_out(shcol,nlev)
  real(rtype), intent(out) :: ql1_out(shcol,nlev)
  real(rtype), intent(out) :: ql2_out(shcol,nlev)

! LOCAL VARIABLES
  integer i,k
  real(rtype) Skew_w,a
  real(rtype) w1_1,w1_2,w2_1,w2_2,w3var
  real(rtype) thl1_1,thl1_2,thl2_1,thl2_2
  real(rtype) qw1_1,qw1_2,qw2_1,qw2_2
  real(rtype) r_qwthl_1
  real(rtype) s1,s2,std_s1,std_s2,C1,C2
  real(rtype) ql1,ql2
  real(rtype) thl_first,qw_first,w_first
  real(rtype) Tl1_1,Tl1_2,pval
  real(rtype) thlsec,qwsec,qwthlsec,wqwsec,wthlsec
  real(rtype) qn1,qn2
  real(rtype) beta1, beta2, qs1, qs2
  real(rtype) sqrtw2, sqrtthl, sqrtqt
  real(rtype) epsterm
  real(rtype) sqrtqw2_1, sqrtqw2_2, sqrtthl2_1, sqrtthl2_2
  real(rtype) thl_tol, rt_tol, w_tol_sqd, w_thresh

  ! variables on thermo grid
  real(rtype) :: wthl_sec_zt(shcol,nlev)
  real(rtype) :: wqw_sec_zt(shcol,nlev)
  real(rtype) :: w3_zt(shcol,nlev)
  real(rtype) :: thl_sec_zt(shcol,nlev)
  real(rtype) :: qwthl_sec_zt(shcol,nlev)
  real(rtype) :: qw_sec_zt(shcol,nlev)

  real(rtype), parameter :: Tl_min = 100._rtype

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
     call shoc_assumed_pdf_f(&
       shcol,nlev,nlevi, &                ! Input
       thetal,qw,w_field,thl_sec,qw_sec,& ! Input
       wthl_sec,w_sec, &                  ! Input
       wqw_sec,qwthl_sec,w3,pres, &       ! Input
       zt_grid,zi_grid,&                  ! Input
       shoc_cldfrac,shoc_ql,&             ! Output
       wqls,wthv_sec,shoc_ql2)            ! Output
    return
  endif
#endif

  epsterm=rgas/rv

  thl_tol=1.e-2_rtype
  rt_tol=1.e-4_rtype
  w_tol_sqd=bfb_square(2.e-2_rtype)
  w_thresh=0.0_rtype

  ! Initialize cloud variables to zero
  shoc_cldfrac(:,:)=0._rtype
  shoc_ql(:,1)=0._rtype
  shoc_ql2(:,:) = 0._rtype

  ! Interpolate many variables from interface grid to themo grid
  call linear_interp(zi_grid,zt_grid,w3,w3_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,thl_sec,thl_sec_zt,nlevi,nlev,shcol,0._rtype)
  call linear_interp(zi_grid,zt_grid,wthl_sec,wthl_sec_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,qwthl_sec,qwthl_sec_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,wqw_sec,wqw_sec_zt,nlevi,nlev,shcol,largeneg)
  call linear_interp(zi_grid,zt_grid,qw_sec,qw_sec_zt,nlevi,nlev,shcol,0._rtype)

  do k=1,nlev
    do i=1,shcol

      pval = pres(i,k)

      ! Get all needed input moments for the PDF
      !  at this particular point
      thl_first = thetal(i,k)
      w_first = w_field(i,k)
      qw_first = qw(i,k)

      w3var = w3_zt(i,k)
      thlsec = thl_sec_zt(i,k)
      qwsec = qw_sec_zt(i,k)
      qwthlsec = qwthl_sec_zt(i,k)
      wqwsec = wqw_sec_zt(i,k)
      wthlsec = wthl_sec_zt(i,k)

      ! Compute square roots of some variables so we don't
      !  have to compute these again
      sqrtw2 = bfb_sqrt(w_sec(i,k))
      sqrtthl = max(thl_tol,bfb_sqrt(thlsec))
      sqrtqt = max(rt_tol,bfb_sqrt(qwsec))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR VERTICAL VELOCITY

      call shoc_assumed_pdf_vv_parameters(&
         w_first,w_sec(i,k),w3var,&    ! Input
         Skew_w,w1_1,w1_2,w2_1,w2_2,a) ! Output
        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR THETAL

      call shoc_assumed_pdf_thl_parameters(&
         wthlsec,sqrtw2,sqrtthl,thlsec,thl_first,& ! Input
         w1_1,w1_2,Skew_w,a,dothetal_skew,&        ! Input
         thl1_1,thl1_2,thl2_1,thl2_2,sqrtthl2_1,&  ! Output
         sqrtthl2_2)                               ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

      call shoc_assumed_pdf_qw_parameters(&
         wqwsec,sqrtw2,Skew_w,sqrtqt,& ! Input
         qwsec,w1_2,w1_1,qw_first,a,&  ! Input
         qw1_1,qw1_2,qw2_1,&           ! Output
         qw2_2,sqrtqw2_1,sqrtqw2_2)    ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  CONVERT FROM TILDE VARIABLES TO "REAL" VARIABLES

      call shoc_assumed_pdf_tilde_to_real(&
         w_first,sqrtw2,& ! Input
         w1_1)            ! Output
      call shoc_assumed_pdf_tilde_to_real(&
         w_first,sqrtw2,& ! Input
         w1_2)            ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  FIND WITHIN-PLUME CORRELATIONS

      call shoc_assumed_pdf_inplume_correlations(&
        sqrtqw2_1,sqrtthl2_1,a,sqrtqw2_2,sqrtthl2_2,&           ! Input
        qwthlsec,qw1_1,qw_first,thl1_1,thl_first,qw1_2,thl1_2,& ! Input
        r_qwthl_1)                                              ! Output

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS

      call shoc_assumed_pdf_compute_temperature(&
        thl1_1,basepres,pval,& ! Input
        Tl1_1)                 ! Output
      call shoc_assumed_pdf_compute_temperature(&
        thl1_2,basepres,pval,& ! Input
        Tl1_2)                 ! Output

      ! Check to ensure Tl1_1 and Tl1_2 are not excessively small.
      !  Temporary fix to set to minimum value if so.
      if (Tl1_1 .le. Tl_min) then
         Tl1_1 = Tl_min
      endif

      if (Tl1_2 .le. Tl_min) then
         Tl1_2 = Tl_min
      endif

      ! Now compute qs
      call shoc_assumed_pdf_compute_qs(&
        Tl1_1,Tl1_2,pval,&   ! Input
        qs1,beta1,qs2,beta2) ! Output

      !!!!!  Now compute cloud stuff
      !!!!!!  compute s term
      call shoc_assumed_pdf_compute_s(&
        qw1_1,qs1,beta1,pval,thl2_1,&          ! Input
        qw2_1,sqrtthl2_1,sqrtqw2_1,r_qwthl_1,& ! Input
        s1,std_s1,qn1,C1)                      ! Output

      !!!!! now compute non-precipitating cloud condensate

      ! If two plumes exactly equal, then just set many of these
      ! variables to themselves to save on computation.
      if (qw1_1 .eq. qw1_2 .and. thl2_1 .eq. thl2_2 .and. qs1 .eq. qs2) then
        s2=s1
        std_s2=std_s1
        C2=C1
        qn2=qn1
      else
        call shoc_assumed_pdf_compute_s(&
        qw1_2,qs2,beta2,pval,thl2_2,&          ! Input
        qw2_2,sqrtthl2_2,sqrtqw2_2,r_qwthl_1,& ! Input
        s2,std_s2,qn2,C2)                      ! Output
      endif

      ql1=min(qn1,qw1_1)
      ql2=min(qn2,qw1_2)

      ! Finally, compute SGS cloud fraction
      shoc_cldfrac(i,k) = min(1._rtype,a*C1+(1._rtype-a)*C2)
      C1_out(i,k) = C1
      C2_out(i,k) = C2
      a1_out(i,k) = a  

      ! Compute SGS liquid water mixing ratio
      call shoc_assumed_pdf_compute_sgs_liquid(&
        a,ql1,ql2,&   ! Input
        shoc_ql(i,k)) ! Output

      ql1_out(i,k) = ql1
      ql2_out(i,k) = ql2

      ! Compute cloud liquid variance (CLUBB formulation, adjusted to SHOC parameters based)
      call shoc_assumed_pdf_compute_cloud_liquid_variance(&
        a,s1,ql1,C1,std_s1,s2,ql2,C2,std_s2,shoc_ql(i,k),& ! Input
        shoc_ql2(i,k))                                     ! Output

      ! Compute liquid water flux
      call shoc_assumed_pdf_compute_liquid_water_flux(&
        a,w1_1,w_first,ql1,w1_2,ql2,& ! Input
        wqls(i,k))                    ! Output

      ! Compute the SGS buoyancy flux
      call shoc_assumed_pdf_compute_buoyancy_flux(&
        wthlsec, epsterm, wqwsec, pval, wqls(i,k),& ! Input
        mf_qlflx_zt(i,k),&
        wthv_sec(i,k))                              ! Output

    enddo  ! end i loop here
  enddo ! end k loop here

  return

end subroutine shoc_assumed_pdf

subroutine shoc_assumed_pdf_vv_parameters(&
   w_first,w_sec,w3var,&          ! Input
   Skew_w,w1_1,w1_2,w2_1,w2_2,a)  ! Output
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND PARAMETERS FOR VERTICAL VELOCITY

  implicit none

  ! intent-ins
  real(rtype), intent(in) :: w_first
  real(rtype), intent(in) :: w_sec
  real(rtype), intent(in) :: w3var

  ! intent-out
  real(rtype), intent(out) :: Skew_w
  real(rtype), intent(out) :: w1_1
  real(rtype), intent(out) :: w1_2
  real(rtype), intent(out) :: w2_1
  real(rtype), intent(out) :: w2_2
  real(rtype), intent(out) :: a

  ! local vars
  real(rtype) :: sqrtw2t

  ! parameters
  real(rtype) :: w_tol_sqd
  w_tol_sqd = bfb_square(2.e-2_rtype)

  Skew_w=w3var/bfb_sqrt(bfb_cube(w_sec))

  if (w_sec .le. w_tol_sqd) then
    Skew_w=0._rtype
    w1_1=w_first
    w1_2=w_first
    w2_1=0._rtype
    w2_2=0._rtype
    a=0.5_rtype
  else

    w2_1=0.4_rtype
    w2_2=0.4_rtype

    a=max(0.01_rtype,min(0.5_rtype*(1._rtype-Skew_w*bfb_sqrt(1._rtype/(4._rtype*bfb_cube(1._rtype-w2_1)+bfb_square(Skew_w)))),0.99_rtype))

    sqrtw2t=bfb_sqrt(1._rtype-w2_1)

    w1_1=bfb_sqrt((1._rtype-a)/a)*sqrtw2t
    w1_2=-1._rtype*bfb_sqrt(a/(1._rtype-a))*sqrtw2t

    w2_1=w2_1*w_sec
    w2_2=w2_2*w_sec

  endif


end subroutine shoc_assumed_pdf_vv_parameters

subroutine shoc_assumed_pdf_thl_parameters(&
  wthlsec,sqrtw2,sqrtthl,thlsec,thl_first,& ! Input
  w1_1,w1_2,Skew_w,a,dothetal_skew,&        ! Input
  thl1_1,thl1_2,thl2_1,thl2_2,sqrtthl2_1,&  ! Output
  sqrtthl2_2)                               ! Output

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO
  implicit none

  ! intent-ins
  real(rtype), intent(in) :: wthlsec
  real(rtype), intent(in) :: sqrtw2
  real(rtype), intent(in) :: sqrtthl
  real(rtype), intent(in) :: thlsec
  real(rtype), intent(in) :: thl_first
  real(rtype), intent(in) :: w1_1
  real(rtype), intent(in) :: w1_2
  real(rtype), intent(in) :: Skew_w
  real(rtype), intent(in) ::  a
  logical(btype), intent(in)  :: dothetal_skew

  ! intent-outs
  real(rtype), intent(out) :: thl1_1
  real(rtype), intent(out) :: thl1_2
  real(rtype), intent(out) :: thl2_1
  real(rtype), intent(out) :: thl2_2
  real(rtype), intent(out) :: sqrtthl2_1
  real(rtype), intent(out) :: sqrtthl2_2

  ! local vars
  real(rtype) :: corrtest1, tsign, Skew_thl

  ! parameters
  real(rtype), parameter :: thl_tol = 1.e-2_rtype
  real(rtype), parameter :: w_thresh = 0.0_rtype

  corrtest1=max(-1._rtype,min(1._rtype,wthlsec/(sqrtw2*sqrtthl)))

  if (thlsec .le. bfb_square(thl_tol) .or. abs(w1_2-w1_1) .le. w_thresh) then
    thl1_1=thl_first
    thl1_2=thl_first
    thl2_1=0._rtype
    thl2_2=0._rtype
    sqrtthl2_1=0._rtype
    sqrtthl2_2=0._rtype
  else

    thl1_1=(-1._rtype*corrtest1)/w1_2
    thl1_2=(-1._rtype*corrtest1)/w1_1

    if (dothetal_skew) then
      tsign=abs(thl1_2-thl1_1)

      if (tsign .gt. 0.4_rtype) then
        Skew_thl=1.2_rtype*Skew_w
      else if (tsign .le. 0.2_rtype) then
        Skew_thl=0.0_rtype
      else
        Skew_thl=((1.2_rtype*Skew_w)/0.2_rtype)*(tsign-0.2_rtype)
      endif
    else
      Skew_thl = 0.0_rtype
    endif

    thl2_1=min(100._rtype,max(0._rtype,(3._rtype*thl1_2*(1._rtype-a*bfb_square(thl1_1)-(1._rtype-a)*bfb_square(thl1_2)) &
            -(Skew_thl-a*bfb_cube(thl1_1)-(1._rtype-a)*bfb_cube(thl1_2)))/ &
            (3._rtype*a*(thl1_2-thl1_1))))*thlsec

    thl2_2=min(100._rtype,max(0._rtype,(-3._rtype*thl1_1*(1._rtype-a*bfb_square(thl1_1)-(1._rtype-a)*bfb_square(thl1_2)) &
      +(Skew_thl-a*bfb_cube(thl1_1)-(1._rtype-a)*bfb_cube(thl1_2)))/ &
      (3._rtype*(1._rtype-a)*(thl1_2-thl1_1))))*thlsec


    thl1_1=thl1_1*sqrtthl+thl_first
    thl1_2=thl1_2*sqrtthl+thl_first

    sqrtthl2_1=bfb_sqrt(thl2_1)
    sqrtthl2_2=bfb_sqrt(thl2_2)

  endif

end subroutine shoc_assumed_pdf_thl_parameters

subroutine shoc_assumed_pdf_qw_parameters(&
  wqwsec, sqrtw2, Skew_w, sqrtqt, qwsec, w1_2, w1_1, qw_first, a, &
  qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

  implicit none

  !intent-in
  real(rtype), intent(in) :: wqwsec
  real(rtype), intent(in) :: sqrtw2
  real(rtype), intent(in) :: Skew_w
  real(rtype), intent(in) :: sqrtqt
  real(rtype), intent(in) :: qwsec
  real(rtype), intent(in) :: w1_2
  real(rtype), intent(in) :: w1_1
  real(rtype), intent(in) :: qw_first
  real(rtype), intent(in) :: a

  ! intent-out
  real(rtype), intent(out) :: qw1_1
  real(rtype), intent(out) :: qw1_2
  real(rtype), intent(out) :: qw2_1
  real(rtype), intent(out) :: qw2_2
  real(rtype), intent(out) :: sqrtqw2_1
  real(rtype), intent(out) :: sqrtqw2_2

  ! local vars
  real(rtype) :: corrtest2, tsign, Skew_qw

  ! Parameters
  real(rtype), parameter :: rt_tol=1.e-4_rtype
  real(rtype), parameter :: w_thresh=0.0_rtype

  corrtest2=max(-1.0_rtype,min(1.0_rtype,wqwsec/(sqrtw2*sqrtqt)))


  if (qwsec .le. bfb_square(rt_tol) .or. abs(w1_2-w1_1) .le. w_thresh) then
    qw1_1=qw_first
    qw1_2=qw_first
    qw2_1=0._rtype
    qw2_2=0._rtype
    sqrtqw2_1=0._rtype
    sqrtqw2_2=0._rtype
  else

    qw1_1=(-1._rtype*corrtest2)/w1_2
    qw1_2=(-1._rtype*corrtest2)/w1_1

    tsign=abs(qw1_2-qw1_1)

    if (tsign .gt. 0.4_rtype) then
      Skew_qw=1.2_rtype*Skew_w
    else if (tsign .le. 0.2_rtype) then
      Skew_qw=0._rtype
    else
      Skew_qw=((1.2_rtype*Skew_w)/0.2_rtype)*(tsign-0.2_rtype)
    endif
    qw2_1=min(100._rtype,max(0._rtype,(3._rtype*qw1_2*(1._rtype-a*bfb_square(qw1_1)-(1._rtype-a)*bfb_square(qw1_2)) &
      -(Skew_qw-a*bfb_cube(qw1_1)-(1._rtype-a)*bfb_cube(qw1_2)))/ &
      (3._rtype*a*(qw1_2-qw1_1))))*qwsec

    qw2_2=min(100._rtype,max(0._rtype,(-3._rtype*qw1_1*(1._rtype-a*bfb_square(qw1_1)-(1._rtype-a)*bfb_square(qw1_2)) &
      +(Skew_qw-a*bfb_cube(qw1_1)-(1._rtype-a)*bfb_cube(qw1_2)))/ &
      (3._rtype*(1._rtype-a)*(qw1_2-qw1_1))))*qwsec

    qw1_1=qw1_1*sqrtqt+qw_first
    qw1_2=qw1_2*sqrtqt+qw_first

    sqrtqw2_1=bfb_sqrt(qw2_1)
    sqrtqw2_2=bfb_sqrt(qw2_2)

  endif

end subroutine shoc_assumed_pdf_qw_parameters

subroutine shoc_assumed_pdf_tilde_to_real(&
  w_first,sqrtw2,& ! intent-in
  w1)              ! intent-out

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  CONVERT FROM TILDE VARIABLES TO "REAL" VARIABLES
  implicit none

  ! intent-ins
  real(rtype), intent(in) :: w_first
  real(rtype), intent(in) :: sqrtw2

  !intent-inouts
  real(rtype), intent(inout) :: w1

  w1 = w1 * sqrtw2 + w_first

end subroutine shoc_assumed_pdf_tilde_to_real

subroutine shoc_assumed_pdf_inplume_correlations(&
  sqrtqw2_1, sqrtthl2_1, a, sqrtqw2_2, sqrtthl2_2,&              ! Input
  qwthlsec, qw1_1, qw_first, thl1_1, thl_first, qw1_2, thl1_2,&  ! Input
  r_qwthl_1)                                                     ! Output

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  FIND WITHIN-PLUME CORRELATIONS
  implicit none

  ! intent out
  real(rtype), intent(in)  :: sqrtqw2_1
  real(rtype), intent(in)  :: sqrtthl2_1
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: sqrtqw2_2
  real(rtype), intent(in)  :: sqrtthl2_2
  real(rtype), intent(in)  :: qwthlsec
  real(rtype), intent(in)  :: qw1_1
  real(rtype), intent(in)  :: qw_first
  real(rtype), intent(in)  :: thl1_1
  real(rtype), intent(in)  :: thl_first
  real(rtype), intent(in)  :: qw1_2
  real(rtype), intent(in)  :: thl1_2

  ! intent out
  real(rtype), intent(out) :: r_qwthl_1

  real(rtype) :: testvar
  testvar=(a*sqrtqw2_1*sqrtthl2_1+(1._rtype-a)*sqrtqw2_2*sqrtthl2_2)

  if (testvar .eq. 0._rtype) then
    r_qwthl_1=0._rtype
  else
    r_qwthl_1=max(-1.0_rtype,min(1.0_rtype,(qwthlsec-a*(qw1_1-qw_first) &
      *(thl1_1-thl_first)-(1._rtype-a)*(qw1_2-qw_first) &
      *(thl1_2-thl_first))/testvar))
  endif

end subroutine shoc_assumed_pdf_inplume_correlations


subroutine shoc_assumed_pdf_compute_temperature(&
  thl1,basepres,pval,& ! Input
  Tl1)                 ! Output

  implicit none
  ! intent-in
  real(rtype), intent(in)  :: thl1
  real(rtype), intent(in)  ::basepres
  real(rtype), intent(in)  ::pval

  ! intent-out
  real(rtype), intent(out) :: Tl1

  TL1 = thl1/(bfb_pow(basepres/pval,(rgas/cp)))

end subroutine shoc_assumed_pdf_compute_temperature

subroutine shoc_assumed_pdf_compute_qs(&
  Tl1_1,Tl1_2,pval,&   ! Input
  qs1,beta1,qs2,beta2) ! Ouput

  use wv_sat_scream, only: MurphyKoop_svp
  implicit none

  ! intent-in
  real(rtype), intent(in) :: Tl1_1
  real(rtype), intent(in) :: Tl1_2
  real(rtype), intent(in) :: pval

  ! intent-out
  real(rtype), intent(out) :: qs1
  real(rtype), intent(out) ::   beta1
  real(rtype), intent(out) ::   qs2
  real(rtype), intent(out) ::   beta2

  ! local vars
  integer, parameter :: liquid = 0   ! liquid flag for MurphyKoop
  real(rtype) :: esval1_1
  real(rtype) :: esval1_2
  real(rtype) :: lstarn1
  real(rtype) :: lstarn2

  esval1_1=0._rtype
  esval1_2=0._rtype

  esval1_1=MurphyKoop_svp(Tl1_1, liquid)
  lstarn1=lcond

  qs1=0.622_rtype*esval1_1/max(esval1_1,pval-esval1_1)
  beta1=(rgas/rv)*(lstarn1/(rgas*Tl1_1))*(lstarn1/(cp*Tl1_1))

  ! Are the two plumes equal?  If so then set qs and beta
  ! in each column to each other to save computation
  lstarn2=lcond
  if (Tl1_1 .eq. Tl1_2) then
    qs2=qs1
    beta2=beta1
  else
    esval1_2=MurphyKoop_svp(Tl1_2, liquid)
    qs2=0.622_rtype*esval1_2/max(esval1_2,pval-esval1_2)
    beta2=(rgas/rv)*(lstarn2/(rgas*Tl1_2))*(lstarn2/(cp*Tl1_2))
  endif

end subroutine shoc_assumed_pdf_compute_qs

subroutine shoc_assumed_pdf_compute_s(&
  qw1,qs1,beta,pval,thl2,&        ! Input
  qw2,sqrtthl2,sqrtqw2,r_qwthl,&  ! Input
  s,std_s,qn,C)                   ! Ouput

  !!!!!!  compute s term
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: qw1
  real(rtype), intent(in)  :: qs1
  real(rtype), intent(in)  :: beta
  real(rtype), intent(in)  :: pval
  real(rtype), intent(in)  :: thl2
  real(rtype), intent(in)  :: qw2
  real(rtype), intent(in)  :: sqrtthl2
  real(rtype), intent(in)  :: sqrtqw2
  real(rtype), intent(in)  :: r_qwthl

  ! intent-out
  real(rtype), intent(out) :: s
  real(rtype), intent(out) :: std_s
  real(rtype), intent(out) :: qn
  real(rtype), intent(out) :: C

  ! local variables
  real(rtype) :: cthl, cqt, tmp_val

  ! Parameters
  real(rtype) :: sqrt2, sqrt2pi
  sqrt2 = bfb_sqrt(2._rtype)
  sqrt2pi = bfb_sqrt(2._rtype*pi)

  s=qw1-qs1*((1._rtype+beta*qw1)/(1._rtype+beta*qs1))
  cthl=((1._rtype+beta*qw1)/bfb_square(1._rtype+beta*qs1))*(cp/lcond) &
    *beta*qs1*bfb_pow(pval/basepres, (rgas/cp))

  cqt=1._rtype/(1._rtype+beta*qs1)
  tmp_val=max(0._rtype,bfb_square(cthl)*thl2+bfb_square(cqt)*qw2-2._rtype*cthl &
    *sqrtthl2*cqt*sqrtqw2*r_qwthl)
  std_s=bfb_sqrt(tmp_val)

  qn=0._rtype
  C=0._rtype

  if (std_s .gt. bfb_sqrt(tiny(1._rtype)) * 100) then
    C=0.5_rtype*(1._rtype+ bfb_erf(s/(sqrt2*std_s)))
    if (C .ne. 0._rtype) qn=s*C+(std_s/sqrt2pi)*bfb_exp(-0.5_rtype*bfb_square(s/std_s))
  else
    if (s .gt. 0._rtype) then
      C=1.0_rtype
      qn=s
    endif
  endif

  ! Prevent possibility of empty clouds or rare occurence of
  !  cloud liquid less than zero
  if (qn .le. 0._rtype) then
    C=0._rtype
    qn=0._rtype
  endif

end subroutine shoc_assumed_pdf_compute_s

subroutine shoc_assumed_pdf_compute_sgs_liquid(&
  a, ql1, ql2,& ! Input
  shoc_ql)      ! Output

  ! Compute SGS liquid water mixing ratio
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: ql1
  real(rtype), intent(in)  :: ql2

  ! intent-out
  real(rtype), intent(out) :: shoc_ql

  shoc_ql = max(0._rtype,a*ql1+(1._rtype-a)*ql2)

end subroutine shoc_assumed_pdf_compute_sgs_liquid

subroutine shoc_assumed_pdf_compute_cloud_liquid_variance(&
  a,s1,ql1,C1,std_s1,s2,ql2,C2,std_s2,shoc_ql,& ! Input
  shoc_ql2)                                     ! Output

  ! Compute cloud liquid variance (CLUBB formulation, adjusted to SHOC parameters based)
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: s1
  real(rtype), intent(in)  :: ql1
  real(rtype), intent(in)  :: C1
  real(rtype), intent(in)  :: std_s1
  real(rtype), intent(in)  :: s2
  real(rtype), intent(in)  :: ql2
  real(rtype), intent(in)  :: C2
  real(rtype), intent(in)  :: std_s2
  real(rtype), intent(in)  :: shoc_ql

  ! intent-out
  real(rtype), intent(out) :: shoc_ql2

  shoc_ql2 = a * ( s1*ql1 + C1*bfb_square(std_s1) )                  &
    + ( 1._rtype-a ) * ( s2*ql2 + C2*bfb_square(std_s2) ) - bfb_square(shoc_ql)
  shoc_ql2 = max( 0._rtype, shoc_ql2 )

end subroutine shoc_assumed_pdf_compute_cloud_liquid_variance

subroutine shoc_assumed_pdf_compute_liquid_water_flux(&
  a,w1_1,w_first,ql1,w1_2,ql2,& ! Input
  wqls)                         ! Output
  ! Compute liquid water flux
  implicit none

  ! intent-in
  real(rtype), intent(in)  :: a
  real(rtype), intent(in)  :: w1_1
  real(rtype), intent(in)  :: w_first
  real(rtype), intent(in)  :: ql1
  real(rtype), intent(in)  :: w1_2
  real(rtype), intent(in)  :: ql2

  ! intent-out
  real(rtype), intent(out) :: wqls

  wqls =a*((w1_1-w_first)*ql1)+(1._rtype-a)*((w1_2-w_first)*ql2)

end subroutine shoc_assumed_pdf_compute_liquid_water_flux
        
subroutine shoc_assumed_pdf_compute_buoyancy_flux(&
  wthlsec,epsterm,wqwsec,pval,wqls,& ! Input
  mf_qlflx_zt,&
  wthv_sec)                          ! Output
  ! Compute the SGS buoyancy flux
  implicit none

  ! intent-in
  real(rtype), intent(in) :: wthlsec
  real(rtype), intent(in) :: epsterm
  real(rtype), intent(in) :: wqwsec
  real(rtype), intent(in) :: pval
  real(rtype), intent(in) :: wqls
  real(rtype), intent(in) :: mf_qlflx_zt

  ! intent-out
  real(rtype), intent(out) :: wthv_sec

  wthv_sec=wthlsec+((1._rtype-epsterm)/epsterm)*basetemp*wqwsec &
  +((lcond/cp)*bfb_pow(basepres/pval,(rgas/cp))-(1._rtype/epsterm)*basetemp)*wqls
  
  !wthv_sec=wthlsec+((1._rtype-epsterm)/epsterm)*basetemp*wqwsec &
  ! +((lcond/cp)*bfb_pow(basepres/pval,(rgas/cp))-(1._rtype/epsterm)*basetemp)*(wqls+mf_qlflx_zt)

end subroutine shoc_assumed_pdf_compute_buoyancy_flux

!==============================================================
! Advance turbulent kinetic energy equation

subroutine shoc_tke(&
         shcol,nlev,nlevi,dtime,&       ! Input
         wthv_sec,shoc_mix,&            ! Input
         dz_zi,dz_zt,pres,tabs,&        ! Input
         u_wind,v_wind,brunt,&          ! Input
         zt_grid,zi_grid,pblh,&         ! Input
         tke,tk,tkh, &                  ! Input/Output
         isotropy, &                    ! Output
         a_diss, a_prod_bu, a_prod_sh)  ! Output

  ! Purpose of this subroutine is to advance the SGS
  !  TKE equation due to shear production, buoyant
  !  production, and dissipation processes.

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels on midpoint grid
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi

  ! timestep [s]
  real(rtype), intent(in) :: dtime
  ! SGS buoyancy flux [K m/s]
  real(rtype), intent(in) :: wthv_sec(shcol,nlev)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! Meridional wind [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! Zonal wind [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! thickness on interface grid [m]
  real(rtype), intent(in) :: dz_zi(shcol,nlevi)
  ! thickness on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! pressure [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! absolute temperature [K]
  real(rtype), intent(in) :: tabs(shcol,nlev)
  ! Brunt Vaisalla frequncy [/s]
  real(rtype), intent(in) :: brunt(shcol,nlev)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! PBLH height
  real(rtype), intent(in) :: pblh(shcol)

! INPUT/OUTPUT VARIABLES
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(inout) :: tk(shcol,nlev)
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(inout) :: tkh(shcol,nlev)

! OUTPUT VARIABLES
  ! Return to isotropic timescale [s]
  real(rtype), intent(out) :: isotropy(shcol,nlev)
  real(rtype), intent(out) :: a_prod_bu(shcol,nlev)
  real(rtype), intent(out) :: a_prod_sh(shcol,nlev)
  real(rtype), intent(out) :: a_diss(shcol,nlev)

! LOCAL VARIABLES
  real(rtype) :: sterm(shcol,nlevi), sterm_zt(shcol,nlev)
  ! Dissipation term
  !real(rtype) :: a_diss(shcol,nlev)
  !column integrated stability
  real(rtype) :: brunt_int(shcol)

  ! Compute integrated column stability in lower troposphere
  call integ_column_stability(nlev, shcol, dz_zt, pres, brunt, brunt_int)

  ! Compute shear production term, which is on interface levels
  ! This follows the methods of Bretheron and Park (2010)

  call compute_shr_prod(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm)

  ! Interpolate shear term from interface to thermo grid
  call linear_interp(zi_grid,zt_grid,sterm,sterm_zt,nlevi,nlev,shcol,0._rtype)

  !advance sgs TKE
  call adv_sgs_tke(nlev, shcol, dtime, shoc_mix, wthv_sec, &
       sterm_zt, tk, tke, a_diss, a_prod_bu, a_prod_sh)

  !Compute isotropic time scale [s]
  call isotropic_ts(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy)

  !Compute eddy diffusivity for heat and momentum
  call eddy_diffusivities(nlev, shcol, pblh, zt_grid, tabs, &
       shoc_mix, sterm_zt, isotropy, tke, tkh, tk)

  return

end subroutine shoc_tke

!==============================================================
! Compute column integrated stability in lower troposphere

subroutine integ_column_stability(nlev, shcol, dz_zt, pres, brunt, brunt_int)

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: integ_column_stability_f
#endif

  implicit none
  !intent-ins
  integer,     intent(in) :: nlev, shcol
  ! thickness on thermodynamic grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! pressure [Pa]
  real(rtype), intent(in) :: pres(shcol,nlev)
  ! Brunt Vaisalla frequncy [/s]
  real(rtype), intent(in) :: brunt(shcol,nlev)

  !intent-out
  !column integrated stability
  real(rtype), intent(out) :: brunt_int(shcol)

  !local variables
  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
     call integ_column_stability_f(nlev, shcol, dz_zt, pres, brunt, brunt_int)
     return
  endif
#endif

  brunt_int(1:shcol) = 0._rtype
  do k = 1, nlev
     do i = 1, shcol
        if (pres(i,k) .gt. troppres) then
           brunt_int(i) = brunt_int(i) + dz_zt(i,k)*brunt(i,k)
        endif
     enddo
  enddo

  return

end subroutine integ_column_stability

!==============================================================
! Compute shear production term, which is on interface levels
! This follows the methods of Bretheron and Park (2010)

subroutine compute_shr_prod(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm)

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_shr_prod_f
#endif

  implicit none

  integer,     intent(in)  :: nlevi, nlev, shcol
  ! thickness on interface grid [m]
  real(rtype), intent(in)  :: dz_zi(shcol,nlevi)
  ! Meridional wind [m/s]
  real(rtype), intent(in)  :: u_wind(shcol,nlev)
  ! Zonal wind [m/s]
  real(rtype), intent(in)  :: v_wind(shcol,nlev)

  !intent-outs
  real(rtype), intent(out) :: sterm(shcol,nlevi)

  !local variables
  integer :: i, k, km1
  real(rtype) :: grid_dz, u_grad, v_grad

  ! Turbulent coefficient
  real(rtype), parameter :: Ck_sh = 0.1_rtype

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
     call compute_shr_prod_f(nlevi, nlev, shcol, dz_zi, u_wind, v_wind, sterm)
     return
  endif
#endif

  !compute shear production term
  do k = 2, nlev
     km1 = k - 1
     do i = 1, shcol
        grid_dz = 1._rtype/dz_zi(i,k)

        ! calculate vertical gradient of u&v wind
        u_grad = grid_dz*(u_wind(i,km1)-u_wind(i,k))
        v_grad = grid_dz*(v_wind(i,km1)-v_wind(i,k))
        sterm(i,k) = Ck_sh*(bfb_square(u_grad)+bfb_square(v_grad))
     enddo
  enddo

  ! Set lower and upper boundary for shear production
  ! Note that the lower bound for shear production has already
  ! been taken into account for the TKE boundary condition,
  ! thus zero out here
  sterm(1:shcol,1)     = 0._rtype
  sterm(1:shcol,nlevi) = 0._rtype

  return

end subroutine compute_shr_prod

!==============================================================
! Advance SGS TKE

subroutine adv_sgs_tke(nlev, shcol, dtime, shoc_mix, wthv_sec, &
     sterm_zt, tk, tke, a_diss, a_prod_bu, a_prod_sh)

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: adv_sgs_tke_f
#endif

  implicit none

  !intent -ins
  integer, intent(in) :: nlev, shcol

  ! timestep [s]
  real(rtype), intent(in) :: dtime
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! SGS buoyancy flux [K m/s]
  real(rtype), intent(in) :: wthv_sec(shcol,nlev)
  ! Interpolate shear production to thermo grid
  real(rtype), intent(in) :: sterm_zt(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(in) :: tk(shcol,nlev)

  ! intent-inout
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(inout) :: tke(shcol,nlev)

  !intent-out
  ! Dissipation term
  real(rtype), intent(out) :: a_diss(shcol,nlev)
  real(rtype), intent(out) :: a_prod_bu(shcol,nlev)
  real(rtype), intent(out) :: a_prod_sh(shcol,nlev)

  !local variables
  integer :: i, k
  !real(rtype) :: a_prod_bu, a_prod_sh

  real(rtype) :: Ck, Cs, Ce, Ce1, Ce2, Cee

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
     call adv_sgs_tke_f(nlev, shcol, dtime, shoc_mix, wthv_sec, &
          sterm_zt, tk, tke, a_diss)
     return
  endif
#endif

  Cs=0.15_rtype
  Ck=0.1_rtype
  Ce=bfb_cube(Ck)/bfb_quad(Cs)  ! MJC: 1.975

  Ce1=Ce/0.7_rtype*0.19_rtype   ! MJC: 0.536
  Ce2=Ce/0.7_rtype*0.51_rtype   ! MJC: 1.439
  Cee=Ce1+Ce2                   ! MJC: 1.975
  !print*,'Cee = ', Cee

  do k = 1, nlev
     do i = 1, shcol

        ! Compute buoyant production term
        a_prod_bu(i,k)=(ggr/basetemp)*wthv_sec(i,k)

        tke(i,k)=max(0._rtype,tke(i,k))

        ! Shear production term, use diffusivity from
        !  previous timestep
        a_prod_sh(i,k)=tk(i,k)*sterm_zt(i,k)

        ! Dissipation term
        !a_diss(i,k)=Cee/shoc_mix(i,k)*bfb_pow(tke(i,k),1.5_rtype)
        ! MJC: Cee as tunable constant
        a_diss(i,k)=Cee_const/shoc_mix(i,k)*bfb_pow(tke(i,k),1.5_rtype)
        
        ! March equation forward one timestep
        tke(i,k)=max(mintke,tke(i,k)+dtime* &
             (max(0._rtype,a_prod_sh(i,k)+a_prod_bu(i,k))-a_diss(i,k)))

        tke(i,k)=min(tke(i,k),maxtke)
     enddo
  enddo

  return

end subroutine adv_sgs_tke

subroutine isotropic_ts(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy)
  !------------------------------------------------------------
  ! Compute the return to isotropic timescale as per
  ! Canuto et al. 2004.  This is used to define the
  ! eddy coefficients as well as to diagnose higher
  ! moments in SHOC
  !------------------------------------------------------------

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: isotropic_ts_f
#endif

  implicit none

  !intent-ins
  integer, intent(in) :: nlev, shcol

  !column integrated stability
  real(rtype), intent(in) :: brunt_int(shcol)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! Dissipation term
  real(rtype), intent(in) :: a_diss(shcol,nlev)
  ! Brunt Vaisalla frequncy [/s]
  real(rtype), intent(in) :: brunt(shcol,nlev)

  ! intent-out
  ! Return to isotropic timescale [s]
  real(rtype), intent(out) :: isotropy(shcol,nlev)

  !local vars
  integer     :: i, k
  real(rtype) :: tscale, lambda, buoy_sgs_save

  !Parameters
  real(rtype), parameter :: maxiso       = 20000.0_rtype ! Return to isotropic timescale [s]

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
     call isotropic_ts_f(nlev, shcol, brunt_int, tke, a_diss, brunt, isotropy)
     return
  endif
#endif


  do k = 1, nlev
     do i = 1, shcol

        ! define the time scale
        tscale=(2.0_rtype*tke(i,k))/a_diss(i,k)

        ! define a damping term "lambda" based on column stability
        lambda=lambda_low+((brunt_int(i)/ggr)-lambda_thresh)*lambda_slope
        lambda=max(lambda_low,min(lambda_high,lambda))

        buoy_sgs_save=brunt(i,k)
        if (buoy_sgs_save .le. 0._rtype) lambda=0._rtype

        ! Compute the return to isotropic timescale
        isotropy(i,k)=min(maxiso,tscale/(1._rtype+lambda*buoy_sgs_save*bfb_square(tscale)))
     enddo
  enddo

  return

end subroutine isotropic_ts

subroutine eddy_diffusivities(nlev, shcol, pblh, zt_grid, tabs, &
     shoc_mix, sterm_zt, isotropy, tke, tkh, tk)

  !------------------------------------------------------------
  ! Compute eddy diffusivity for heat and momentum
  !------------------------------------------------------------

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: eddy_diffusivities_f
#endif

  implicit none

  !intent-ins
  integer, intent(in) :: nlev, shcol

  ! PBL height [m]
  real(rtype), intent(in) :: pblh(shcol)
  ! Heights on the mid-point grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! Absolute temperature [K]
  real(rtype), intent(in) :: tabs(shcol,nlev)
  ! Mixing length [m]
  real(rtype), intent(in) :: shoc_mix(shcol,nlev)
  ! Interpolate shear production to thermo grid
  real(rtype), intent(in) :: sterm_zt(shcol,nlev)
  ! Return to isotropic timescale [s]
  real(rtype), intent(in) :: isotropy(shcol,nlev)
  ! turbulent kinetic energy [m2/s2]
  real(rtype), intent(in) :: tke(shcol,nlev)

  !intent-outs
  ! eddy coefficient for heat [m2/s]
  real(rtype), intent(out) :: tkh(shcol,nlev)
  ! eddy coefficient for momentum [m2/s]
  real(rtype), intent(out) :: tk(shcol,nlev)

  !local vars
  integer     :: i, k

  !parameters
  ! Minimum absolute temperature threshold for which to apply extra mixing [K]
  real(rtype), parameter :: temp_crit = 182.0_rtype
  ! Transition depth [m] above PBL top to allow
  ! stability diffusivities
  real(rtype), parameter :: pbl_trans = 200.0_rtype

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call eddy_diffusivities_f(nlev, shcol, pblh, zt_grid, tabs, &
                                shoc_mix, sterm_zt, isotropy, tke, tkh, tk)
      return
   endif
#endif

  do k = 1, nlev
     do i = 1, shcol

        if (tabs(i,nlev) .lt. temp_crit .and. (zt_grid(i,k) .lt. pblh(i)+pbl_trans)) then
           ! If surface layer temperature is running away, apply extra mixing
	   !   based on traditional stable PBL diffusivities that are not damped
	   !   by stability functions.

           tkh(i,k) = Ckh_s*bfb_square(shoc_mix(i,k))*bfb_sqrt(sterm_zt(i,k))
           tk(i,k)  = Ckm_s*bfb_square(shoc_mix(i,k))*bfb_sqrt(sterm_zt(i,k))
        else
           ! Default definition of eddy diffusivity for heat and momentum
           tkh(i,k) = Ckh*isotropy(i,k)*tke(i,k)
           tk(i,k)  = Ckm*isotropy(i,k)*tke(i,k)
        endif

     enddo
  enddo

  return

end subroutine eddy_diffusivities

!==============================================================
! Check the TKE

subroutine check_tke(&
             shcol,nlev,& ! Input
             tke)         ! Input/Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: check_tke_f
#endif

  implicit none
  ! Make sure TKE falls within reasonable bounds
  ! If not, then clip

! INPUT VARIABLES
  integer, intent(in) :: shcol
  integer, intent(in) :: nlev

! IN/OUT VARIABLES
  real(rtype), intent(inout) :: tke(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call check_tke_f(shcol,nlev, & ! Input
           tke)                      ! Input/Output
      return
   endif
#endif

  do k=1,nlev
    do i=1,shcol
      tke(i,k)=max(mintke,tke(i,k))
    enddo
  enddo

end subroutine check_tke

!==============================================================
! Compute the turbulent length scale

subroutine shoc_length(&
         shcol,nlev,nlevi,&            ! Input
         host_dx,host_dy,&             ! Input
         zt_grid,zi_grid,dz_zt,&       ! Input
         tke,thv,&                     ! Input
         brunt,l_inf,shoc_mix)         ! Output

  ! Purpose of this subroutine is to compute the SHOC
  !  mixing length scale, which is used to compute the
  !  turbulent dissipation in the SGS TKE equation

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: shoc_length_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  !  number of midpoint levels
  integer, intent(in) :: nlev
  ! number of interface levels
  integer, intent(in) :: nlevi

  ! host model grid size [m]
  real(rtype), intent(in) :: host_dx(shcol)
  ! host model grid size [m]
  real(rtype), intent(in) :: host_dy(shcol)
  ! turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! dz on midpoint grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! virtual potential temperature [K]
  real(rtype), intent(in) :: thv(shcol,nlev)

  ! OUTPUT VARIABLES
  ! Brunt-Vaisala frequency [/s]
  real(rtype), intent(out) :: brunt(shcol,nlev)
  ! SHOC mixing length [m]
  real(rtype), intent(out) :: shoc_mix(shcol,nlev)
  real(rtype), intent(out) :: l_inf(shcol)

  ! LOCAL VARIABLES
  real(rtype) :: thv_zi(shcol,nlevi)
  !real(rtype) :: l_inf(shcol)

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_length_f(shcol,nlev,nlevi,host_dx,host_dy,&
                         zt_grid,zi_grid,dz_zt,tke,&
                         thv,brunt,shoc_mix)

      return
   endif
#endif

  ! Interpolate virtual potential temperature onto interface grid
  call linear_interp(zt_grid,zi_grid,thv,thv_zi,nlev,nlevi,shcol,0._rtype)

  ! Define the brunt vaisalia frequency
  call compute_brunt_shoc_length(nlev,nlevi,shcol,dz_zt,thv,thv_zi,brunt)

  ! Find L_inf
  call compute_l_inf_shoc_length(nlev,shcol,zt_grid,dz_zt,tke,l_inf)

  ! compute mixing-length
  call compute_shoc_mix_shoc_length(nlev,shcol,tke,brunt,zt_grid,l_inf,shoc_mix)

  ! Do checks on the length scale.  Make sure it is not
  !  larger than the grid mesh of the host model.
  call check_length_scale_shoc_length(nlev,shcol,host_dx,host_dy,shoc_mix)

  return

end subroutine shoc_length

!==============================================================
! Subroutine to determine superdiagonal and the factorized
! subdiagonal and diagonal coeffs for solving the
! tridiagonal diffusion matrix.

!subroutine vd_shoc_decomp( &
!         shcol,nlev,nlevi,&          ! Input
!         kv_term,tmpi,rdp_zt,dtime,& ! Input
!         flux, &                     ! Input
!         du,dl,d)                    ! Output
!
!  implicit none
!
!! INPUT VARIABLES
!  ! number of columns
!  integer, intent(in) :: shcol
!  ! number of mid-point levels
!  integer, intent(in) :: nlev
!  ! number of levels on the interface
!  integer, intent(in) :: nlevi
!
!  ! SHOC timestep [s]
!  real(rtype), intent(in) :: dtime
!  ! diffusion coefficent [m2/s]
!  real(rtype), intent(in) :: kv_term(shcol,nlevi)
!  ! dt*(g*rho)**2/dp at interfaces
!  real(rtype), intent(in) :: tmpi(shcol,nlevi)
!  ! 1/dp
!  real(rtype), intent(in) :: rdp_zt(shcol,nlev)
!  ! surface flux [varies]
!  real(rtype), intent(in) :: flux(shcol)
!
!! OUTPUT VARIABLES
!  ! superdiagonal
!  real(rtype), intent(out) :: du(shcol,nlev)
!  ! Factorized version of subdiagonal
!  real(rtype), intent(out) :: dl(shcol,nlev)
!  ! Factorized version of diagonal
!  real(rtype), intent(out) :: d(shcol,nlev)
!
!! LOCAL VARIABLES
!  integer :: i, k
!
!  ! Determine superdiagonal (du) and subdiagonal (dl) coeffs of the
!  ! tridiagonal diffusion matrix.
!  do k=1,nlev-1
!    do i=1,shcol
!      du(i,k)   = -1._rtype * kv_term(i,k+1) * tmpi(i,k+1) * rdp_zt(i,k)
!      dl(i,k+1) = -1._rtype * kv_term(i,k+1) * tmpi(i,k+1) * rdp_zt(i,k+1)
!    enddo
!  enddo
!
!  ! The bottom element of the superdiagonal (du) and the top element of
!  ! the subdiagonal (dl) is set to zero (not included in linear system).
!  du(:,nlev) = 0._rtype
!  dl(:,1)    = 0._rtype
!
!  ! Compute the diagonal and perform Thomas factorization. The diagonal
!  ! elements are a combination of du and dl (d=1-du-dl). Surface fluxes
!  ! are applied explicitly in the diagonal at the top level.
!  do i=1,shcol
!    d(i,1) = 1._rtype - du(i,1)
!  enddo
!  do k=2,nlev-1
!    do i=1,shcol
!      d(i,k) = 1._rtype - du(i,k) - dl(i,k)
!
!      dl(i,k) = dl(i,k)/d(i,k-1)
!      d (i,k) = d (i,k) - dl(i,k)*du(i,k-1)
!    enddo
!  enddo
!  do i=1,shcol
!    d(i,nlev) = 1._rtype - dl(i,nlev) + flux(i)*dtime*ggr*rdp_zt(i,nlev)
!
!    dl(i,nlev) = dl(i,nlev)/d(i,nlev-1)
!    d (i,nlev) = d (i,nlev) - dl(i,nlev)*du(i,nlev-1)
!  enddo
!
!  return
!
!end subroutine vd_shoc_decomp

! MJC: Previous version (SCREAMv0)
subroutine vd_shoc_decomp( &
         shcol,nlev,nlevi,&          ! Input
         kv_term,tmpi,rdp_zt,dtime,& ! Input
         flux, &                     ! Input
         do_mf,mf_ae,mf_aw,tmpi3,&   ! EDMF input
         ca,cc,denom,ze)             ! Output

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of mid-point levels
  integer, intent(in) :: nlev
  ! number of levels on the interface
  integer, intent(in) :: nlevi

  ! SHOC timestep [s]
  real(rtype), intent(in) :: dtime
  ! diffusion coefficent [m2/s]
  real(rtype), intent(in) :: kv_term(shcol,nlevi)
  ! dt*(g*rho)**2/dp at interfaces
  real(rtype), intent(in) :: tmpi(shcol,nlevi)
  ! 1/dp
  real(rtype), intent(in) :: rdp_zt(shcol,nlev)
  ! surface flux [varies]
  real(rtype), intent(in) :: flux(shcol)

  ! MJC: EDMF inputs
  ! Include mass flux contribution?
  logical,  intent(in)  :: do_mf
  ! Sum of environment area, i.e. 1-sum(a_i) [-]
  real(rtype), intent(in)  :: mf_ae(shcol,nlevi)
  ! Sum (a_i*w_i) [m/s]
  real(rtype), intent(in)  :: mf_aw(shcol,nlevi)
  ! dt*g*rho on interfaces
  real(rtype), intent(in)  :: tmpi3(shcol,nlevi)
! OUTPUT VARIABLES
  ! superdiagonal
  real(rtype), intent(out) :: ca(shcol,nlev)
  ! subdiagonal
  real(rtype), intent(out) :: cc(shcol,nlev)
  ! 1./(1.+ca(k)+cc(k)-cc(k)*ze(k-1))
  real(rtype), intent(out) :: denom(shcol,nlev)
  ! Term in tri-diag. matrix system
  real(rtype), intent(out) :: ze(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k

  ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the
  ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
  ! a combination of ca and cc; they are not required by the solver.

  do k=nlev-1,1,-1
    do i=1,shcol
      if ( do_mf ) then
        ca(i,k  ) = (mf_ae(i,k+1)*kv_term(i,k+1)*tmpi(i,k+1) &
             - 0.5_rtype*tmpi3(i,k+1)*mf_aw(i,k+1))*rdp_zt(i,k  )
        cc(i,k+1) = (mf_ae(i,k+1)*kv_term(i,k+1)*tmpi(i,k+1) &
             + 0.5_rtype*tmpi3(i,k+1)*mf_aw(i,k+1))*rdp_zt(i,k+1)
      else
        ca(i,k) = kv_term(i,k+1) * tmpi(i,k+1) * rdp_zt(i,k)
        cc(i,k+1) = kv_term(i,k+1) * tmpi(i,k+1) * rdp_zt(i,k+1)
      endif
    enddo
  enddo

  ! The bottom element of the upper diagonal (ca) is zero (not used).
  ! The subdiagonal (cc) is not needed in the solver.

  ca(:,nlev) = 0._rtype

  ! Calculate e(k). This term is required in the solution of the
  ! tridiagonal matrix as defined by the implicit diffusion equation.

  do i=1,shcol
    if ( do_mf ) then
      denom(i,nlev) = 1._rtype/ &
        (1._rtype + cc(i,nlev) - &
        tmpi3(i,nlev)*mf_aw(i,nlev)*rdp_zt(i,nlev) + &
        flux(i)*dtime*ggr*rdp_zt(i,nlev))
    else
      denom(i,nlev) = 1._rtype/ &
        (1._rtype + cc(i,nlev) + flux(i)*dtime*ggr*rdp_zt(i,nlev))
    endif
    ze(i,nlev) = cc(i,nlev) * denom(i,nlev)
  enddo

  do k=nlev-1,2,-1
    do i=1,shcol
      if ( do_mf ) then
        denom(i,k) = 1._rtype/ &
          (1._rtype + ca(i,k) + cc(i,k) - &
          ca(i,k)*ze(i,k+1) - &
          (tmpi3(i,k)*mf_aw(i,k) - tmpi3(i,k+1)*mf_aw(i,k+1))*rdp_zt(i,k))
      else
        denom(i,k) = 1._rtype/ &
          (1._rtype + ca(i,k) + cc(i,k) - &
          ca(i,k) * ze(i,k+1))
      endif
      ze(i,k) = cc(i,k) * denom(i,k)
    enddo
  enddo

  do i=1,shcol
    if ( do_mf ) then
      ! MJC: bug fix
      !denom(i,1) = 1._rtype/ &
      !  (1._rtype + ca(i,1) - ca(i,1)*ze(i,2) - &
      !  (tmpi3(i,1)*mf_aw(i,1) - tmpi3(i,2)*mf_aw(i,2))*rdp_zt(i,1))
      denom(i,1) = 1._rtype/ &
        (1._rtype + ca(i,1) - ca(i,1)*ze(i,2) + &
        tmpi3(i,2)*mf_aw(i,2)*rdp_zt(i,1))
    else
      denom(i,1) = 1._rtype/ &
        (1._rtype + ca(i,1) - ca(i,1) * ze(i,2))
    endif
  enddo

  return

end subroutine vd_shoc_decomp


!==============================================================
! Subroutine to solve the implicit vertical diffsion equation
! with zero flux boundary conditions. Actual surface fluxes
! should be applied explicitly. Procedure for solution of the
! implicit equation follows Richtmeyer and Morton (1967, pp 198-200).
!
! Here, du is the superdiagonal and dl and d are the factorized
! version of the subdiagonal and diagonal using the Thomas algorithm.
! This subroutine takes in an input profile var and solves using the
! remainder of the Thomas algorithm. The solution is stored in var.
!
! Note that the same routine is used for temperature, momentum and
! tracers.
! ---------------------------------------------------------------
!
!subroutine vd_shoc_solve(&
!         shcol,nlev,& ! Input
!         du,dl,d,&    ! Input
!         var)         ! Input/Output
!
!  implicit none
!
!! INPUT VARIABLES
!  ! number of columns
!  integer, intent(in) :: shcol
!  ! number of mid-point levels
!  integer, intent(in) :: nlev
!  ! superdiagonal
!  real(rtype), intent(in) :: du(shcol,nlev)
!  ! Factorized version of subdiagonal
!  real(rtype), intent(in) :: dl(shcol,nlev)
!  ! Factorized version of diagonal
!  real(rtype), intent(in) :: d(shcol,nlev)
!
!! IN/OUT VARIABLES
!  real(rtype), intent(inout) :: var(shcol,nlev)
!
!! LOCAL VARIABLES
!  integer :: i, k
!
!  ! Solve using Thomas algorithm
!  do k=2,nlev
!    do i=1,shcol
!      var(i,k) = var(i,k) - dl(i,k)*var(i,k-1)
!    enddo
!  enddo
!  do i=1,shcol
!    var(i,nlev) = var(i,nlev)/d(i,nlev)
!  enddo
!  do k=nlev,2,-1
!    do i=1,shcol
!      var(i,k-1) = (var(i,k-1) - du(i,k-1)*var(i,k))/d(i,k-1)
!    enddo
!  enddo
!
!  return
!
!end subroutine vd_shoc_solve

! MJC: Previous version (SCREAMv0)
subroutine vd_shoc_solve(&
         shcol,nlev,nlevi,&   ! Input
         ca,cc,denom,ze,&     ! Input
         do_mf,mf_awvar,&     ! EDMF Input
         tmpi3,rdp_zt,&       ! EDMF Input
         var)                 ! Input/Output

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of mid-point levels
  integer, intent(in) :: nlev
  ! number of levels on the interface
  integer, intent(in) :: nlevi
  ! superdiagonal
  real(rtype), intent(in) :: ca(shcol,nlev)
  ! subdiagonal
  real(rtype), intent(in) :: cc(shcol,nlev)
  ! 1./(1.+ca(k)+cc(k)-cc(k)*ze(k-1))
  real(rtype), intent(in) :: denom(shcol,nlev)
  ! Term in tri-diag. matrix system
  real(rtype), intent(in) :: ze(shcol,nlev)

  ! MJC: EDMF inputs
  ! Include mass flux contribution?
  logical,  intent(in)    :: do_mf
  ! Sum of plume (a_i*w_i*var_i)
  real(rtype), intent(in)    :: mf_awvar(shcol,nlevi)
  ! dt*g*rho on interfaces
  real(rtype), intent(in)    :: tmpi3(shcol,nlevi)
  ! 1/dp
  real(rtype), intent(in)    :: rdp_zt(shcol,nlev)

! IN/OUT VARIABLES
  real(rtype), intent(inout) :: var(shcol,nlev)

! LOCAL VARIABLES
  integer :: i, k
  ! Term in tri-diag solution
  real(rtype) :: zf(shcol,nlev)


  ! Calculate zf(k). Terms zf(k) and ze(k) are required in solution of
  ! tridiagonal matrix defined by implicit diffusion equation.
  ! Note that only levels ntop through nbot need be solved for.

  do i=1,shcol
    if (do_mf) then
      zf(i,nlev) = (var(i,nlev) - &
        tmpi3(i,nlev)*mf_awvar(i,nlev)*rdp_zt(i,nlev))*denom(i,nlev)
    else
      zf(i,nlev) = var(i,nlev) * denom(i,nlev)
    endif 
  enddo

  do k=nlev-1,1,-1
    do i=1,shcol
      if (do_mf) then
        zf(i,k) = (var(i,k) + &
          (tmpi3(i,k+1)*mf_awvar(i,k+1) - tmpi3(i,k)*mf_awvar(i,k))*rdp_zt(i,k) + &
          ca(i,k)*zf(i,k+1))*denom(i,k)
      else
        zf(i,k) = (var(i,k) + ca(i,k) * zf(i,k+1)) * denom(i,k)
      endif 
    enddo
  enddo

  ! Perform back substitution

  do i=1,shcol
    var(i,1) = zf(i,1)
  enddo

  do k=2,nlev
    do i=1,shcol
      var(i,k) = zf(i,k) + ze(i,k)*var(i,k-1)
    enddo
  enddo

  return

end subroutine vd_shoc_solve

!==============================================================
! Subroutine to compute integrals for SHOC conservation
!  with host model

subroutine shoc_energy_integrals(&
         shcol,nlev,host_dse,pdel,&     ! Input
         rtm,rcm,u_wind,v_wind,&        ! Input
         se_int,ke_int,wv_int,wl_int)   ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_energy_integrals_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! host model temperature [K]
  real(rtype), intent(in) :: host_dse(shcol,nlev)
  ! pressure layer thickness [Pa]
  real(rtype), intent(in) :: pdel(shcol,nlev)
  ! zonal wind [m/s]
  real(rtype), intent(in) :: u_wind(shcol,nlev)
  ! meridional wind [m/s]
  real(rtype), intent(in) :: v_wind(shcol,nlev)
  ! total water mixing ratio [kg/kg]
  real(rtype), intent(in) :: rtm(shcol,nlev)
  ! cloud liquid mixing ratio [kg/kg]
  real(rtype), intent(in) :: rcm(shcol,nlev)

! OUTPUT VARIABLES
  ! integrated static energy
  real(rtype), intent(out) :: se_int(shcol)
  ! integrated kinetic energy
  real(rtype), intent(out) :: ke_int(shcol)
  ! integrated water vapor
  real(rtype), intent(out) :: wv_int(shcol)
  ! integrated liquid water
  real(rtype), intent(out) :: wl_int(shcol)

! LOCAL VARIABLES
  integer :: i, k
  real(rtype) :: rvm

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_energy_integrals_f(shcol,nlev,host_dse,pdel,&   ! Input
                                   rtm,rcm,u_wind,v_wind,&      ! Input
                                   se_int,ke_int,wv_int,wl_int) ! Output
      return
   endif
#endif

  se_int(:) = 0._rtype
  ke_int(:) = 0._rtype
  wv_int(:) = 0._rtype
  wl_int(:) = 0._rtype
  do k=1,nlev
    do i=1,shcol
       rvm = rtm(i,k) - rcm(i,k) ! compute water vapor

!technically wrong, need to remove gz from geopotential
!but shoc does not change gz term
       se_int(i) = se_int(i) + host_dse(i,k)*pdel(i,k)/ggr
       ke_int(i) = ke_int(i) + 0.5_rtype*(bfb_square(u_wind(i,k))+bfb_square(v_wind(i,k)))*pdel(i,k)/ggr
       wv_int(i) = wv_int(i) + rvm*pdel(i,k)/ggr
       wl_int(i) = wl_int(i) + rcm(i,k)*pdel(i,k)/ggr
    enddo
  enddo

  return

end subroutine shoc_energy_integrals

!==============================================================
! Subroutine to update SHOC output to host model temperature

subroutine update_host_dse(&
         shcol,nlev,thlm,&                 ! Input
         shoc_ql,inv_exner,zt_grid,phis,&  ! Input
         host_dse)                         ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: update_host_dse_f
#endif

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! liquid water potential temperature [K]
  real(rtype), intent(in) :: thlm(shcol,nlev)
  ! cloud liquid water mixing ratio [kg/kg]
  real(rtype), intent(in) :: shoc_ql(shcol,nlev)
  ! inverse of exner function [-]
  real(rtype), intent(in) :: inv_exner(shcol,nlev)
  ! heights at mid point [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! surface geopotential height of host model [m]
  real(rtype), intent(in) :: phis(shcol)

  ! OUTPUT VARIABLES
  ! host model temperature [K]
  real(rtype), intent(out) :: host_dse(shcol,nlev)

  ! LOCAL VARIABLES
  ! Temperature [K]
  real(rtype) :: temp

  integer :: i, k

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call update_host_dse_f(shcol,nlev,thlm,shoc_ql,inv_exner,zt_grid,phis, &  ! Input
           host_dse)                              ! Input/Output)
      return
   endif
#endif

  do k=1,nlev
    do i=1,shcol
      temp = (thlm(i,k)/inv_exner(i,k))+(lcond/cp)*shoc_ql(i,k)
      host_dse(i,k) = cp*temp+ggr*zt_grid(i,k)+phis(i)
    enddo
  enddo

  return

end subroutine update_host_dse

!==============================================================
! Subroutine for SHOC energy fixer with host model temp

subroutine shoc_energy_fixer(&
         shcol,nlev,nlevi,dtime,nadv,&  ! Input
         zt_grid,zi_grid,&              ! Input
         se_b,ke_b,wv_b,wl_b,&          ! Input
         se_a,ke_a,wv_a,wl_a,&          ! Input
         wthl_sfc,wqw_sfc,&             ! Input
         rho_zt,tke,pint,&              ! Input
         host_dse)                      ! Input/Output

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: shoc_energy_fixer_f
#endif

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi
  ! SHOC timestep
  real(rtype), intent(in) :: dtime
  ! number of SHOC iterations
  integer, intent(in) :: nadv
  ! integrated static energy before
  real(rtype), intent(in) :: se_b(shcol)
  ! integrated kinetic energy before
  real(rtype), intent(in) :: ke_b(shcol)
  ! integrated water vapor before
  real(rtype), intent(in) :: wv_b(shcol)
  ! integrated liquid water before
  real(rtype), intent(in) :: wl_b(shcol)
  ! integrated static energy after
  real(rtype), intent(in) :: se_a(shcol)
  ! integrated kinetic energy after
  real(rtype), intent(in) :: ke_a(shcol)
  ! integrated water vapor after
  real(rtype), intent(in) :: wv_a(shcol)
  ! integrated liquid water after
  real(rtype), intent(in) :: wl_a(shcol)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! pressure on interface grid [Pa]
  real(rtype), intent(in) :: pint(shcol,nlevi)
  ! density on midpoint grid [kg/m^3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)
  !turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)

  ! INPUT VARIABLES
  !host temperature [K]
  real(rtype), intent(inout) :: host_dse(shcol,nlev)

  ! LOCAL VARIABLES
  real(rtype) :: se_dis(shcol), te_a(shcol), te_b(shcol)
  integer :: shoctop(shcol)

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_energy_fixer_f(shcol,nlev,nlevi,dtime,nadv,& ! Input
                              zt_grid,zi_grid,&              ! Input
                              se_b,ke_b,wv_b,wl_b,&          ! Input
                              se_a,ke_a,wv_a,wl_a,&          ! Input
                              wthl_sfc,wqw_sfc,&             ! Input
                              rho_zt,tke,pint,&              ! Input
                              host_dse)                      ! Input/Output
      return
   endif
#endif

  call shoc_energy_total_fixer(&
         shcol,nlev,nlevi,dtime,nadv,&  ! Input
         zt_grid,zi_grid,&              ! Input
         se_b,ke_b,wv_b,wl_b,&          ! Input
         se_a,ke_a,wv_a,wl_a,&          ! Input
         wthl_sfc,wqw_sfc,rho_zt,pint,& ! Input
         te_a, te_b)                    ! Output

  call shoc_energy_threshold_fixer(&
         shcol,nlev,nlevi,&             ! Input
         pint,tke,te_a,te_b,&           ! Input
         se_dis,shoctop)                ! Output

  call shoc_energy_dse_fixer(&
         shcol,nlev,&                  ! Input
         se_dis,shoctop,   &           ! Input
         host_dse)                     ! Input/Output

  return

end subroutine shoc_energy_fixer

!==============================================================
! Subroutine foe SHOC energy fixer with host model temp

subroutine shoc_energy_total_fixer(&
         shcol,nlev,nlevi,dtime,nadv,&  ! Input
         zt_grid,zi_grid,&              ! Input
         se_b,ke_b,wv_b,wl_b,&          ! Input
         se_a,ke_a,wv_a,wl_a,&          ! Input
         wthl_sfc,wqw_sfc,rho_zt,pint,& ! Input
         te_a, te_b)                    ! Output

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi
  ! SHOC timestep
  real(rtype), intent(in) :: dtime
  ! number of SHOC iterations
  integer, intent(in) :: nadv
  ! integrated static energy before
  real(rtype), intent(in) :: se_b(shcol)
  ! integrated kinetic energy before
  real(rtype), intent(in) :: ke_b(shcol)
  ! integrated water vapor before
  real(rtype), intent(in) :: wv_b(shcol)
  ! integrated liquid water before
  real(rtype), intent(in) :: wl_b(shcol)
  ! integrated static energy after
  real(rtype), intent(in) :: se_a(shcol)
  ! integrated kinetic energy after
  real(rtype), intent(in) :: ke_a(shcol)
  ! integrated water vapor after
  real(rtype), intent(in) :: wv_a(shcol)
  ! integrated liquid water after
  real(rtype), intent(in) :: wl_a(shcol)
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! heights on midpoint grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  ! heights on interface grid [m]
  real(rtype), intent(in) :: zi_grid(shcol,nlevi)
  ! density on midpoint grid [kg/m^3]
  real(rtype), intent(in) :: rho_zt(shcol,nlev)
  ! pressure on interface grid [Pa]
  real(rtype), intent(in) :: pint(shcol,nlevi)

  ! OUTPUT VARIABLES
  real(rtype), intent(out) :: te_a(shcol)
  real(rtype), intent(out) :: te_b(shcol)

  ! LOCAL VARIABLES
  ! density on interface grid [kg/m^3]
  real(rtype) :: rho_zi(shcol,nlevi)
  ! sensible and latent heat fluxes [W/m^2]
  real(rtype) :: shf, lhf, hdtime, exner_surf
  integer :: i

  ! compute the host timestep
  hdtime = dtime * float(nadv)

  call linear_interp(zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,shcol,0._rtype)

  ! Based on these integrals, compute the total energy before and after SHOC
  ! call
  do i=1,shcol
    ! convert shf and lhf
    exner_surf = bfb_pow(pint(i,nlevi)/p0, rgas/cp)
    shf=wthl_sfc(i)*cp*rho_zi(i,nlevi)*exner_surf

    lhf=wqw_sfc(i)*rho_zi(i,nlevi)
    te_a(i) = se_a(i) + ke_a(i) + (lcond+lice)*wv_a(i)+lice*wl_a(i)
    te_b(i) = se_b(i) + ke_b(i) + (lcond+lice)*wv_b(i)+lice*wl_b(i)
    te_b(i) = te_b(i)+(shf+(lhf)*(lcond+lice))*hdtime
  enddo

  return

end subroutine shoc_energy_total_fixer



!==============================================================
! Subroutine foe SHOC energy fixer with host model temp

subroutine shoc_energy_threshold_fixer(&
         shcol,nlev,nlevi,&             ! Input
         pint,tke,te_a,te_b,&           ! Input
         se_dis,shoctop)                ! Output

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev
  ! number of levels on interface grid
  integer, intent(in) :: nlevi
  ! pressure on interface grid [Pa]
  real(rtype), intent(in) :: pint(shcol,nlevi)
  !turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)

  real(rtype), intent(in) :: te_a(shcol)
  real(rtype), intent(in) :: te_b(shcol)


  ! OUTPUT VARIABLES
  real(rtype), intent(out) :: se_dis(shcol)
  integer, intent(out) :: shoctop(shcol)

  ! LOCAL VARIABLES
  ! sensible and latent heat fluxes [W/m^2]
  integer :: i

  ! Limit the energy fixer to find highest layer where SHOC is active
  ! Find first level where wp2 is higher than lowest threshold
  do i=1,shcol
    shoctop(i) = 1
    do while (tke(i,shoctop(i)) .eq. mintke .and. shoctop(i) .lt. nlev-1)
      shoctop(i) = shoctop(i) + 1
    enddo

    ! Compute the disbalance of total energy, over depth where SHOC is active
    se_dis(i) = (te_a(i) - te_b(i))/(pint(i,nlevi)-pint(i,shoctop(i)))
  enddo

  return

end subroutine shoc_energy_threshold_fixer

!==============================================================
! Subroutine foe SHOC energy fixer with host model temp

subroutine shoc_energy_dse_fixer(&
         shcol,nlev,&                  ! Input
         se_dis,shoctop,   &           ! Input
         host_dse)                     ! Input/Output

  implicit none

  ! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! number of levels
  integer, intent(in) :: nlev

  ! INPUT VARIABLES
  real(rtype), intent(in) :: se_dis(shcol)
  integer, intent(in) :: shoctop(shcol)

  !host temperature [K]
  real(rtype), intent(inout) :: host_dse(shcol,nlev)

  ! LOCAL VARIABLES
  integer :: i, k

  do i=1,shcol
    do k=shoctop(i),nlev
      host_dse(i,k) = host_dse(i,k) - se_dis(i)*ggr
    enddo
  enddo

  return

end subroutine shoc_energy_dse_fixer

!==============================================================
! Linear interpolation to get values on various grids

subroutine shoc_diag_obklen(&
         shcol,uw_sfc,vw_sfc,&      ! Input
         wthl_sfc,wqw_sfc,thl_sfc,& ! Input
         cldliq_sfc,qv_sfc,&        ! Input
         ustar,kbfs,obklen)         ! Ouput

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: shoc_diag_obklen_f
#endif

  implicit none

! INPUT VARIABLES
  ! number of columns
  integer, intent(in) :: shcol
  ! Surface sensible heat flux [K m/s]
  real(rtype), intent(in) :: wthl_sfc(shcol)
  ! Surface latent heat flux [kg/kg m/s]
  real(rtype), intent(in) :: wqw_sfc(shcol)
  ! Surface momentum flux (u-direction) [m2/s2]
  real(rtype), intent(in) :: uw_sfc(shcol)
  ! Surface momentum flux (v-direction) [m2/s2]
  real(rtype), intent(in) :: vw_sfc(shcol)
  ! Surface potential temperature [K]
  real(rtype), intent(in) :: thl_sfc(shcol)
  ! Surface cloud liquid water [kg /kg]
  real(rtype), intent(in) :: cldliq_sfc(shcol)
  ! Surface water vapor
  real(rtype), intent(in) :: qv_sfc(shcol)

! OUTPUT VARIABLES
  ! Obukhov length [m]
  real(rtype), intent(out) :: obklen(shcol)
  ! surface friction velocity [m/s]
  real(rtype), intent(out) :: ustar(shcol)
  ! surface kinematic buoyancy flux [m^s/s^3]
  real(rtype), intent(out) :: kbfs(shcol)

! LOCAL VARIABLES
  integer :: i
  real(rtype) :: th_sfc  ! potential temperature at surface
  real(rtype) :: thv_sfc ! virtual potential temperature at surface

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call shoc_diag_obklen_f(shcol,uw_sfc,vw_sfc,&      ! Input
                            wthl_sfc,wqw_sfc,thl_sfc,& ! Input
                            cldliq_sfc,qv_sfc,&        ! Input
                            ustar,kbfs,obklen)         ! Ouput
    return
  endif
#endif

  do i=1,shcol
    th_sfc = thl_sfc(i) + (lcond/cp)*cldliq_sfc(i)
    thv_sfc = th_sfc*(1._rtype+eps*qv_sfc(i)-cldliq_sfc(i))
    ustar(i) = max(bfb_sqrt(bfb_square(uw_sfc(i)) + bfb_square(vw_sfc(i))),ustar_min)
    kbfs(i) = wthl_sfc(i)+eps*th_sfc*wqw_sfc(i)
    obklen(i) = -thv_sfc*bfb_cube(ustar(i))/(ggr*vk*(kbfs(i)+sign(1.e-10_rtype,kbfs(i))))
  enddo

  return

end subroutine shoc_diag_obklen

  !
  !===============================================================================
subroutine pblintd(&
       shcol,nlev,nlevi,&             ! Input
       z,zi,thl,ql,&                  ! Input
       q,u,v,&                        ! Input
       ustar,obklen,kbfs,cldn,&       ! Input
       pblh)                          ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: pblintd_f
#endif

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Diagnose standard PBL variables
    !
    ! Method:
    ! Diagnosis of PBL depth.
    ! The PBL depth follows:
    !    Holtslag, A.A.M., and B.A. Boville, 1993:
    !    Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
    !    Model. J. Clim., vol. 6., p. 1825--1842.
    !
    ! Updated by Holtslag and Hack to exclude the surface layer from the
    ! definition of the boundary layer Richardson number. Ri is now defined
    ! across the outer layer of the pbl (between the top of the surface
    ! layer and the pbl top) instead of the full pbl (between the surface and
    ! the pbl top). For simiplicity, the surface layer is assumed to be the
    ! region below the first model level (otherwise the boundary layer depth
    ! determination would require iteration).
    !
    ! Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
    ! (Use ricr = 0.3 in this formulation)
    !
    ! Author: B. Stevens (extracted from pbldiff, August 2000)
    !
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: thl(shcol,nlev)         ! liquid water potential temp [K]
    real(rtype), intent(in)  :: ql(shcol,nlev)          ! cloud liquid mixing ratio [kg/kg]
    real(rtype), intent(in)  :: q(shcol,nlev)           ! water vapor [kg/kg]
    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: u(shcol,nlev)           ! windspeed x-direction [m/s]
    real(rtype), intent(in)  :: v(shcol,nlev)           ! windspeed y-direction [m/s]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    real(rtype), intent(in)  :: obklen(shcol)           ! Obukhov length
    real(rtype), intent(in)  :: kbfs(shcol)             ! sfc kinematic buoyancy flux [m^2/s^3]
    real(rtype), intent(in)  :: zi(shcol,nlevi)         ! height above surface [m]
    real(rtype), intent(in)  :: cldn(shcol,nlev)        ! new cloud fraction
    !
    ! Output arguments
    !
    real(rtype), intent(out) :: pblh(shcol)             ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !

    real(rtype) :: rino(shcol,nlev)        ! bulk Richardson no. from level to ref lev
    real(rtype) :: thv(shcol,nlev)         ! virtual potential temperature
    real(rtype) :: tlv(shcol)              ! ref. level pot tmp + tmp excess

    logical(btype)  :: check(shcol)            ! True=>chk if Richardson no.>critcal

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call pblintd_f(&
      shcol,nlev,nlevi,npbl,&        ! Input
      z,zi,thl,ql,&                  ! Input
      q,u,v,&                        ! Input
      ustar,obklen,kbfs,cldn,&       ! Input
      pblh)                          ! Output
    return
  endif
#endif

    !
    ! Compute Obukhov length virtual temperature flux and various arrays for use later:
    !

    ! Compute virtual potential temperature
    call pblintd_init_pot(&
               shcol,nlev,&             ! Input
               thl,ql,q,&               ! Input
               thv)                     ! Output

    call pblintd_init(&
           shcol,nlev,&             ! Input
           z,&                      ! Input
           check,rino,pblh)         ! Output

    !
    ! PBL height calculation
    !
    call pblintd_height(&
       shcol,nlev,&                   ! Input
       z,u,v,ustar,&                  ! Input
       thv,thv(:,nlev),&              ! Input
       pblh,rino,check)               ! Output
    !
    ! Estimate an effective surface temperature to account for surface
    ! fluctuations
    !
    call pblintd_surf_temp(&
       shcol,nlev,nlevi,&          ! Input
       z,ustar,obklen,kbfs,thv,&   ! Input
       tlv,&                       ! Output
       pblh,check,rino)            ! InOutput
    !
    ! Improve pblh estimate for unstable conditions using the convective
    ! temperature excess as reference temperature:
    !
    call pblintd_height(&
       shcol,nlev,&             ! Input
       z,u,v,ustar,&            ! Input
       thv,tlv,&                ! Input
       pblh,rino,check)         ! Output
    !
    ! Check PBL height
    !
    call pblintd_check_pblh(&
       shcol,nlev,nlevi,&             ! Input
       z,ustar,check,&                ! Input
       pblh)                          ! Output
    !
    ! PBL check over ocean
    !
    call pblintd_cldcheck(      &
                   shcol,nlev,nlevi, &                  ! Input
                   zi,cldn,          &                  ! Input
                   pblh)                                ! InOutput

    return
end subroutine pblintd

subroutine pblintd_init_pot(&
       shcol,nlev,&             ! Input
       thl,ql,q,&               ! Input
       thv)                     ! Output
#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_pblintd_init_pot_f
#endif
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers

    real(rtype), intent(in)  :: thl(shcol,nlev)         ! liquid water potential temp [K]
    real(rtype), intent(in)  :: ql(shcol,nlev)          ! cloud liquid mixing ratio [kg/kg]
    real(rtype), intent(in)  :: q(shcol,nlev)           ! water vapor [kg/kg]
    real(rtype), intent(out) :: thv(shcol,nlev)         ! virtual potential temperature
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index
    real(rtype) :: th

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_pblintd_init_pot_f(shcol,nlev,thl,ql,q,&               ! Input
                                   thv)                     ! Output
      return
   endif
#endif
    ! Compute virtual potential temperature
    do k=1,nlev
      do i=1,shcol
        th=thl(i,k)+(lcond/cp)*ql(i,k)
        thv(i,k)=th*(1._rtype+eps*q(i,k)-ql(i,k))
      enddo
    enddo

end subroutine pblintd_init_pot

subroutine pblintd_init(&
       shcol,nlev,&             ! Input
       z,&                      ! Input
       check,rino,pblh)         ! Output

    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    !
    ! Output arguments
    !
    real(rtype), intent(out) :: pblh(shcol)             ! boundary-layer height [m]
    real(rtype), intent(out) :: rino(shcol,nlev)        ! bulk Richardson no. from level to ref lev
    logical(btype), intent(out)     :: check(shcol)            ! True=>chk if Richardson no.>critcal

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index

    do i=1,shcol
       check(i)     = .true.
       rino(i,nlev) = 0.0_rtype
       pblh(i)      = z(i,nlev)
    end do
end subroutine pblintd_init

subroutine pblintd_height(&
       shcol,nlev,&              ! Input
       z,u,v,ustar,&             ! Input
       thv,thv_ref,&             ! Input
       pblh,rino,check)          ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: pblintd_height_f
#endif

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: u(shcol,nlev)           ! windspeed x-direction [m/s]
    real(rtype), intent(in)  :: v(shcol,nlev)           ! windspeed y-direction [m/s]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    real(rtype), intent(in)  :: thv(shcol,nlev)         ! virtual potential temperature
    real(rtype), intent(in)  :: thv_ref(shcol)          ! ref. level pot tmp

    !
    ! Output arguments
    !
    real(rtype), intent(inout)   :: pblh(shcol)             ! boundary-layer height [m]
    real(rtype), intent(inout) :: rino(shcol,nlev)                     ! bulk Richardson no. from level to ref lev
    logical(btype), intent(inout)     :: check(shcol)            ! True=>chk if Richardson no.>critcal

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    integer  :: k                       ! level index
    real(rtype) :: vvk                     ! velocity magnitude squared

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call pblintd_height_f(shcol,nlev,npbl,z,u,v,ustar,thv,thv_ref,& ! Input
                            pblh,rino,check)                          ! Output
      return
   endif
#endif

    !
    ! PBL height calculation:  Scan upward until the Richardson number between
    ! the first level and the current level exceeds the "critical" value.
    !
    do k=nlev-1,nlev-npbl+1,-1
       do i=1,shcol
          if (check(i)) then
             vvk = bfb_square((u(i,k) - u(i,nlev))) + bfb_square((v(i,k) - v(i,nlev))) + fac*bfb_square(ustar(i))
             vvk = max(vvk,tinyw)
             rino(i,k) = ggr*(thv(i,k) -thv_ref(i))*(z(i,k)-z(i,nlev))/(thv(i,nlev)*vvk)
             if (rino(i,k) >= ricr) then
                pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/(rino(i,k) -rino(i,k+1)) * &
                     (z(i,k) - z(i,k+1))
                check(i) = .false.
             endif
          endif
       end do
    end do
    return
end subroutine pblintd_height

subroutine pblintd_surf_temp(&
       shcol,nlev,nlevi,&          ! Input
       z,ustar,obklen,kbfs,thv,&   ! Input
       tlv,&                       ! Output
       pblh,check,rino)            ! InOutput

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: pblintd_surf_temp_f
#endif

    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    real(rtype), intent(in)  :: obklen(shcol)           ! Obukhov length
    real(rtype), intent(in)  :: kbfs(shcol)             ! sfc kinematic buoyancy flux [m^2/s^3]
    real(rtype), intent(in) :: thv(shcol,nlev)          ! virtual potential temperature

    real(rtype), intent(out) :: tlv(shcol)              ! ref. level pot tmp + tmp excess
    logical(btype), intent(inout)  :: check(shcol)             ! True=>chk if Richardson no.>critcal
    real(rtype), intent(inout) :: rino(shcol,nlev)      ! bulk Richardson no. from level to ref lev
    real(rtype), intent(inout) :: pblh(shcol)              ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !
    real(rtype) :: phiminv
    integer  :: i                       ! longitude index

    !===================
    ! const parameter for Diagnosis of PBL depth
    real(rtype), parameter :: fak   =  8.5_rtype      ! Constant in surface temperature excess
    real(rtype), parameter :: betam = 15.0_rtype      ! Constant in wind gradient expression
    real(rtype), parameter :: sffrac=  0.1_rtype      ! Surface layer fraction of boundary layer
    real(rtype), parameter :: binm  = betam*sffrac ! betam * sffrac

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call pblintd_surf_temp_f(shcol,nlev,nlevi,&          ! Input
                         z,ustar,obklen,kbfs,thv,&   ! Input
                         tlv,pblh,check,rino)            ! InOutput
      return
   endif
#endif

    !
    ! Estimate an effective surface temperature to account for surface
    ! fluctuations
    !
    do i=1,shcol
       if (check(i)) pblh(i) = z(i,nlevi-npbl)
       check(i)  = (kbfs(i) > 0._rtype)
       tlv(i)    = thv(i,nlev)
       if (check(i)) then
          phiminv      = bfb_cbrt(1._rtype - binm*pblh(i)/obklen(i))
          rino(i,nlev) = 0.0_rtype
          tlv(i)       = thv(i,nlev) + kbfs(i)*fak/( ustar(i)*phiminv )
       end if
    end do
    return
end subroutine pblintd_surf_temp

subroutine pblintd_check_pblh(&
       shcol,nlev,nlevi,&             ! Input
       z,ustar,check,&                ! Input
       pblh)                          ! Output

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: pblintd_check_pblh_f
#endif

    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: z(shcol,nlev)           ! height above surface [m]
    real(rtype), intent(in)  :: ustar(shcol)            ! surface friction velocity [m/s]
    logical(btype), intent(in)      :: check(shcol)            ! True=>chk if Richardson no.>critcal
    !
    ! Output arguments
    !
    real(rtype), intent(inout) :: pblh(shcol)             ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call pblintd_check_pblh_f(shcol,nlev,nlevi,npbl,z,ustar,check,pblh)
      return
   endif
#endif

    !
    ! PBL height must be greater than some minimum mechanical mixing depth
    ! Several investigators have proposed minimum mechanical mixing depth
    ! relationships as a function of the local friction velocity, u*.  We
    ! make use of a linear relationship of the form h = c u* where c=700.
    ! The scaling arguments that give rise to this relationship most often
    ! represent the coefficient c as some constant over the local coriolis
    ! parameter.  Here we make use of the experimental results of Koracin
    ! and Berkowicz (1988) [BLM, Vol 43] for which they recommend 0.07/f
    ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
    ! latitude value for f so that c = 0.07/f = 700.  Also, do not allow
    ! PBL to exceed some maximum (npbl) number of allowable points
    !
    do i=1,shcol
       if (check(i)) pblh(i) = z(i,nlevi-npbl)
       pblh(i) = max(pblh(i),700.0_rtype*ustar(i))
    end do
    return
end subroutine pblintd_check_pblh


subroutine pblintd_cldcheck(      &
                   shcol,nlev,nlevi, &                  ! Input
                   zi,cldn,          &                  ! Input
                   pblh)                                ! InOutput

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: shoc_pblintd_cldcheck_f
#endif

    !------------------------------Arguments--------------------------------
    ! Input arguments
    !
    integer, intent(in) :: shcol                     ! number of atmospheric columns
    integer, intent(in) :: nlev                      ! number of mid-point layers
    integer, intent(in) :: nlevi                     ! number of interface layers

    real(rtype), intent(in)  :: zi(shcol,nlevi)         ! height above surface [m]
    real(rtype), intent(in)  :: cldn(shcol,nlev)        ! new cloud fraction
    !
    ! In/Output arguments
    !
    real(rtype), intent(inout) :: pblh(shcol)             ! boundary-layer height [m]
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i                       ! longitude index
    logical(btype)  :: cldcheck(shcol)      ! True=>if cloud in lowest layer

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call shoc_pblintd_cldcheck_f(shcol, nlev, nlevi, zi, cldn, pblh)
      return
   endif
#endif

    !
    ! Final requirement on PBL heightis that it must be greater than the depth
    ! of the lowest model level if there is any cloud diagnosed in
    ! the lowest model level.  This is to deal with the inadequacies of the
    ! current "dry" formulation of the boundary layer, where this test is
    ! used to identify circumstances where there is marine stratus in the
    ! lowest level, and to provide a weak ventilation of the layer to avoid
    ! a pathology in the cloud scheme (locking in low-level stratiform cloud)
    ! If  any cloud is diagnosed in the lowest level, set pblh to 50 meters
    ! higher than top interface of lowest level
    !
    do i=1,shcol
       cldcheck(i) = .false.
       if (cldn(i,nlev).ge.0.0_rtype) cldcheck(i) = .true.
       if (cldcheck(i)) pblh(i) = max(pblh(i),zi(i,nlev) + 50._rtype)
    end do
    return
end subroutine pblintd_cldcheck

  !==============================================================
  ! Linear interpolation to get values on various grids

subroutine linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)

#ifdef SCREAM_CONFIG_IS_CMAKE
    use shoc_iso_f, only: linear_interp_f
#endif

    implicit none

    integer, intent(in) :: km1, km2
    integer, intent(in) :: ncol
    real(rtype), intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real(rtype), intent(in) :: x2(ncol,km2)
    real(rtype), intent(in) :: minthresh
    real(rtype), intent(out) :: y2(ncol,km2)

    integer :: k1, k2, i

#ifdef SCREAM_CONFIG_IS_CMAKE
   if (use_cxx) then
      call linear_interp_f(x1,x2,y1,y2,km1,km2,ncol,minthresh)
      return
   endif
#endif

#if 1
    !i = check_grid(x1,x2,km1,km2,ncol)
    if (km1 .eq. km2+1) then
       do k2 = 1,km2
          k1 = k2+1
          do i = 1,ncol
             y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
          end do
       end do
    elseif (km2 .eq. km1+1) then
       k2 = 1
       do i = 1,ncol
          y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
       end do
       do k2 = 2, km2-1
          k1 = k2
          do i = 1,ncol
             y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
          end do
       end do
       k2 = km2
       do i = 1,ncol
          y2(i,k2) = y1(i,km1-1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1-1))/(x1(i,km1)-x1(i,km1-1))
       end do
    else
       !print *,km1,km2
    end if
    do k2 = 1,km2
       do i = 1,ncol
          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif
       end do
    end do
#else
    do i=1,ncol
       do k2=1,km2
          if( x2(i,k2) <= x1(i,1) ) then
             y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
          elseif( x2(i,k2) >= x1(i,km1) ) then
             y2(i,k2) = y1(i,km1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1))/(x1(i,km1)-x1(i,km1-1))
          else
             do k1 = 2,km1
                if( (x2(i,k2)>=x1(i,k1-1)).and.(x2(i,k2)<x1(i,k1)) ) then
                   y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
                endif
             enddo ! end k1 loop
          endif

          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif

       enddo ! end k2 loop
    enddo ! i loop
#endif

    return

end subroutine linear_interp

subroutine compute_brunt_shoc_length(nlev,nlevi,shcol,dz_zt,thv,thv_zi,brunt)

  !=========================================================
  !
  ! Computes the brunt_visala frequency

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_brunt_shoc_length_f
#endif

  implicit none
  integer, intent(in) :: nlev, nlevi, shcol
  ! Grid difference centereted on thermo grid [m]
  real(rtype), intent(in) :: dz_zt(shcol,nlev)
  ! virtual potential temperature [K]
  real(rtype), intent(in) :: thv(shcol,nlev)
  ! virtual potential temperature [K] at interface
  real(rtype), intent(in) :: thv_zi(shcol,nlevi)
  ! brunt vaisala frequency [s-1]
  real(rtype), intent(out) :: brunt(shcol, nlev)
  integer k, i

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call compute_brunt_shoc_length_f(nlev,nlevi,shcol,dz_zt,thv,thv_zi,brunt)
    return
  endif
#endif

  do k=1,nlev
    do i=1,shcol
      brunt(i,k) = (ggr/thv(i,k)) * (thv_zi(i,k) - thv_zi(i,k+1))/dz_zt(i,k)
    enddo
  enddo

end subroutine compute_brunt_shoc_length

subroutine compute_l_inf_shoc_length(nlev,shcol,zt_grid,dz_zt,tke,l_inf)

  !=========================================================
  !

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_l_inf_shoc_length_f
#endif

  implicit none
  integer, intent(in) :: nlev, shcol
  real(rtype), intent(in) :: zt_grid(shcol,nlev), dz_zt(shcol,nlev), tke(shcol,nlev)
  real(rtype), intent(inout) :: l_inf(shcol)
  real(rtype) :: tkes, numer(shcol), denom(shcol)
  integer k, i

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call compute_l_inf_shoc_length_f(nlev,shcol,zt_grid,dz_zt,tke,l_inf)
    return
  endif
#endif

  numer(:) = 0._rtype
  denom(:) = 0._rtype
  
  ! MJC [10/31/24]: Add do_edmf to use original code when running just SHOC.
  do k=1,nlev
    do i=1,shcol
      if (do_edmf) then
        if ( tke(i,k) .ge. 0.0005_rtype ) then
           tkes=bfb_sqrt(tke(i,k))
        else
           tkes = 0._rtype    
        endif
        numer(i)=numer(i)+tkes*zt_grid(i,k)*dz_zt(i,k)
        denom(i)=denom(i)+tkes*dz_zt(i,k)
      else
        tkes=bfb_sqrt(tke(i,k))
        numer(i)=numer(i)+tkes*zt_grid(i,k)*dz_zt(i,k)
        denom(i)=denom(i)+tkes*dz_zt(i,k)
      endif  
    enddo
  enddo
  
  do i=1,shcol
    if (do_edmf) then
      if (denom(i) .eq. 0._rtype) then
        l_inf(i)=l_inf_const     
      else
        l_inf(i)=0.1_rtype*(numer(i)/denom(i))
      endif
    else
      l_inf(i)=0.1_rtype*(numer(i)/denom(i))
    endif
  enddo

end subroutine compute_l_inf_shoc_length

subroutine compute_shoc_mix_shoc_length(nlev,shcol,tke,brunt,zt_grid,l_inf,shoc_mix)

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: compute_shoc_mix_shoc_length_f
#endif

  implicit none

  integer, intent(in) :: nlev, shcol
  ! turbulent kinetic energy [m^2/s^2]
  real(rtype), intent(in) :: tke(shcol,nlev)
  ! brunt vaisala frequency [s-1]
  real(rtype), intent(in) :: brunt(shcol,nlev)
  ! heights, for thermo grid [m]
  real(rtype), intent(in) :: zt_grid(shcol,nlev)
  real(rtype), intent(in) :: l_inf(shcol)

  ! Turbulent length scale [m]
  real(rtype), intent(out) :: shoc_mix(shcol,nlev)

  !  LOCAL VARIABLES
  real(rtype) :: brunt2(shcol,nlev)
  integer k, i
  real(rtype) :: tkes

  ! Turnover timescale [s]
  real(rtype), parameter :: tscale = 400._rtype

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call compute_shoc_mix_shoc_length_f(nlev,shcol,tke,brunt,zt_grid,l_inf,& !Input
                                        shoc_mix) ! Ouptut
    return
  endif
#endif

  brunt2(:,:) = 0.0

  do k=1,nlev
    do i=1,shcol

      tkes = sqrt(tke(i,k))

      if(brunt(i,k) .ge. 0) brunt2(i,k) = brunt(i,k)
 
      ! MJC [10/31/24]: Add do_edmf to use original code when running just SHOC.
      if (do_edmf) then
        ! MJC: Original + l_inf_const  
        shoc_mix(i,k)=min(maxlen,(2.8284_rtype*sqrt(1._rtype/((1._rtype/(tscale*tkes*vk*zt_grid(i,k)))&
        +(1._rtype/(tscale*tkes*l_inf_const))+0.01_rtype*(brunt2(i,k)/tke(i,k)))))/length_fac)
      else
        shoc_mix(i,k)=min(maxlen,(2.8284_rtype*sqrt(1._rtype/((1._rtype/(tscale*tkes*vk*zt_grid(i,k)))&
        +(1._rtype/(tscale*tkes*l_inf(i)))+0.01_rtype*(brunt2(i,k)/tke(i,k)))))/length_fac)
      endif
          
    enddo ! end i loop (column loop)
  enddo ! end k loop (vertical loop)

end subroutine compute_shoc_mix_shoc_length

subroutine check_length_scale_shoc_length(nlev,shcol,host_dx,host_dy,shoc_mix)
  ! Do checks on the length scale.  Make sure it is not
  !  larger than the grid mesh of the host model.

#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_iso_f, only: check_length_scale_shoc_length_f
#endif

  implicit none
  integer, intent(in) :: nlev, shcol
  real(rtype), intent(in) :: host_dx(shcol), host_dy(shcol)
  ! Turbulent length scale [m]
  real(rtype), intent(inout) :: shoc_mix(shcol, nlev)
  integer k, i

#ifdef SCREAM_CONFIG_IS_CMAKE
  if (use_cxx) then
    call check_length_scale_shoc_length_f(nlev,shcol,host_dx,host_dy,shoc_mix)
    return
  endif
#endif

  do k=1,nlev
    do i=1,shcol
      shoc_mix(i,k)=min(maxlen,shoc_mix(i,k))
      shoc_mix(i,k)=max(minlen,shoc_mix(i,k))
      shoc_mix(i,k)=min(bfb_sqrt(host_dx(i)*host_dy(i)),shoc_mix(i,k))
    enddo
  enddo

end subroutine check_length_scale_shoc_length

!=========================================================
!=========================================================
! MJC: Poisson code for MF entrainment
!call Poisson( nz, mf_nup, entf, enti, 69._rtype)

  subroutine poisson(nz,nup,lambda,poi,state)
  !**********************************************************************
  ! Set a unique (but reproduceble) seed for the kiss RNG
  ! Call Poisson deviate
  ! By Adam Herrington
  !**********************************************************************
#ifdef SCREAM_CONFIG_IS_CMAKE
   use shoc_eam_host_stubs, only: ShrKissRandGen
#else
   use shr_RandNum_mod, only: ShrKissRandGen
#endif

       integer,                     intent(in)  :: nz,nup
       real(rtype), dimension(1,4),      intent(in)  :: state
       real(rtype), dimension(nz,nup), intent(in)  :: lambda
       integer,  dimension(nz,nup), intent(out) :: poi
       integer,  dimension(1,4)                 :: tmpseed
       integer                                  :: i,j
       type(ShrKissRandGen)                     :: kiss_gen

       ! Compute seed
       tmpseed(1,1) = int((state(1,1) - int(state(1,1))) * 1000000000._rtype)
       tmpseed(1,2) = int((state(1,2) - int(state(1,2))) * 1000000000._rtype)
       tmpseed(1,3) = int((state(1,3) - int(state(1,3))) * 1000000000._rtype)
       tmpseed(1,4) = int((state(1,4) - int(state(1,4))) * 1000000000._rtype)

       ! Set seed
       kiss_gen = ShrKissRandGen(tmpseed)

       do i=1,nz
         do j=1,nup
           call hybridRNG(kiss_gen,lambda(i,j),poi(i,j))
         enddo
       enddo

  end subroutine poisson

  subroutine hybridRNG(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Interface for the two poisson rng subroutines
  ! chooses the appropriate subroutine based on the value of lambda
  !**********************************************************************
#ifdef SCREAM_CONFIG_IS_CMAKE
   use shoc_eam_host_stubs, only: ShrKissRandGen
#else
   use shr_RandNum_mod, only: ShrKissRandGen
#endif

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(rtype),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       if (lambda < 10._rtype) then
          call knuth(kiss_gen,lambda,kout)
       else
          call hormann(kiss_gen,lambda,kout)
       end if

  end subroutine hybridRNG

  subroutine knuth(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson from Knuth 
  ! The Art of Computer Programming, v2, 137-138
  ! By Adam Herrington
  !**********************************************************************
#ifdef SCREAM_CONFIG_IS_CMAKE
   use shoc_eam_host_stubs, only: ShrKissRandGen
#else
   use shr_RandNum_mod, only: ShrKissRandGen
#endif

       type(ShrKissRandGen), intent(inout) :: kiss_gen
       real(rtype),             intent(in)    :: lambda
       integer,              intent(out)   :: kout

       ! Local variables
       real(rtype), dimension(1,1) :: tmpuni
       real(rtype)                 :: puni, explam
       integer                  :: k

       k = 0
       explam = exp(-1._rtype*lambda)
       puni = 1._rtype
       do while (puni > explam)
         k = k + 1
         call kiss_gen%random(tmpuni)
         puni = puni*tmpuni(1,1)
       end do
       kout = k - 1

  end subroutine knuth

  subroutine hormann(kiss_gen,lambda,kout)
  !**********************************************************************
  ! Discrete random poisson
  ! Implements Poisson Transformed Rejection with Squeeze (PTRS) 
  ! from W. Hormann Insurance: Mathematics and Economics 12, 39-45 (1993) 
  ! By Jake Reschke
  !**********************************************************************
#ifdef SCREAM_CONFIG_IS_CMAKE
  use shoc_eam_host_stubs, only: ShrKissRandGen
#else
  use shr_RandNum_mod, only: ShrKissRandGen
#endif

      type(ShrKissRandGen), intent(inout) :: kiss_gen
      real(rtype),             intent(in)    :: lambda
      integer,              intent(out)   :: kout

      ! Local variables
      real(rtype), dimension(1,1) :: U,V
      real(rtype)                 :: a,b,vr,alphinv,us,loggam
      integer                  :: k,i

      b = 0.931_rtype + 2.53_rtype*sqrt(lambda)
      a = -0.059_rtype + 0.02483_rtype*b
      vr = 0.9277_rtype - 3.6224_rtype/(b - 2._rtype)
      alphinv = 1.1239_rtype + 1.1328_rtype/(b - 3.4_rtype)

      do
         call kiss_gen%random(U)
         call kiss_gen%random(V)
         U(1,1) = U(1,1) - 0.5_rtype
         us = 0.5_rtype - abs(U(1,1))
         k = floor( (2._rtype*a/us + b)*U(1,1) + lambda + 0.43_rtype )
         if (us >= 0.07_rtype .and.  V(1,1) <= vr) then
            kout = k
            exit
         end if
         if (k <= 0 .or. (us < 0.013_rtype .and. V(1,1) > us)) then
            cycle
         end if
         ! compute log(k!). If k >=10 use stirling's approximation
         if (k < 10) then
            loggam = 0._rtype
            do i = 1, k
               loggam = loggam + log(1._rtype*i)
            end do
         else
            loggam = log(sqrt(2._rtype*pi)) + (k + 0.5_rtype)*log(1._rtype*k) - k + (1._rtype/12._rtype - 1._rtype/(360._rtype*k*k))/k
         end if
         if (log( V(1,1)*alphinv/(a/(us*us) + b) ) <= -1._rtype*lambda + k*log(lambda) - loggam) then
            kout = k
            exit
         end if
      end do

  end subroutine hormann


end module

!==============================================================
! This is the end of the SHOC parameterization
! We hope you have enjoyed your time here
