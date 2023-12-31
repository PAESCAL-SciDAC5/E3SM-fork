module clubb_intr_types

  use shr_kind_mod, only: r8 => shr_kind_r8 

  implicit none

  public

  type clubb_mean_2d_t

   integer :: ncol

   real(r8), pointer, dimension(:,:) :: um       ! mean east-west wind                          [m/s]
   real(r8), pointer, dimension(:,:) :: vm       ! mean north-south wind                        [m/s]

   real(r8), pointer, dimension(:,:) :: thlm     ! mean liquid water potential temperature      [K]
   real(r8), pointer, dimension(:,:) :: rtm      ! mean total water (vapor+liquid) mixing ratio [kg/kg]

   real(r8), pointer, dimension(:,:) :: rcm      ! mean cloud (liquid) water mixing ratio       [kg/kg]

   real(r8), pointer, dimension(:,:) :: exner_clubb  ! CLUBB's exner function                       [-]
   real(r8), pointer, dimension(:,:) :: pmid         ! air pressure                                 [Pa]
   real(r8), pointer, dimension(:,:) :: pdel         ! layer thickness in pressure coordinate       [Pa]
   real(r8), pointer, dimension(:,:) :: t            ! air temperature                              [K]
   real(r8), pointer, dimension(:,:) :: zi           ! geopotential height at layer interface       [m]
   real(r8), pointer, dimension(:,:) :: zm           ! geopotential height at layer midpoint        [m]
   real(r8), pointer, dimension(:,:) :: omega        ! vertical velocity in pressure coordinate [Pa/s]
   real(r8), pointer, dimension(:)   :: lat          ! latitude  [rad]
   real(r8), pointer, dimension(:)   :: lon          ! longitude [rad]

   real(r8), pointer, dimension(:,:,:) :: q          ! tracer mixing ratio  [-]

   real(r8), pointer, dimension(:,:) :: qrl
   real(r8), pointer, dimension(:,:) :: prer_evap 

  end type clubb_mean_2d_t

  type clubb_mnts_2d_t

   real(r8), pointer, dimension(:,:) :: wp2      ! vertical velocity variance                   [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: wp3      ! third moment of vertical velocity            [m^3/s^3]
   real(r8), pointer, dimension(:,:) :: wpthlp   ! turbulent flux of thetal                     [m/s K]
   real(r8), pointer, dimension(:,:) :: wprtp    ! turbulent flux of moisture                   [m/s kg/kg]
   real(r8), pointer, dimension(:,:) :: rtpthlp  ! covariance of thetal and qt                  [kg/kg K]
   real(r8), pointer, dimension(:,:) :: rtp2     ! moisture variance                            [kg^2/kg^2]
   real(r8), pointer, dimension(:,:) :: thlp2    ! temperature variance                         [K^2]
   real(r8), pointer, dimension(:,:) :: up2      ! east-west wind variance                      [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vp2      ! north-south wind variance                    [m^2/s^2]

   real(r8), pointer, dimension(:,:) :: wpthvp     ! < w'th_v' > (momentum levels)                [m/s K]
   real(r8), pointer, dimension(:,:) :: wp2thvp    ! < w'^2 th_v' > (thermodynamic levels)        [m^2/s^2 K]
   real(r8), pointer, dimension(:,:) :: rtpthvp    ! < r_t'th_v' > (momentum levels)              [kg/kg K]
   real(r8), pointer, dimension(:,:) :: thlpthvp   ! < th_l'th_v' > (momentum levels)             [K^2]

   real(r8), pointer, dimension(:,:) :: upwp     ! east-west momentum flux                      [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vpwp     ! north-south momentum flux                    [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: um_pert  ! perturbed meridional wind                    [m/s]
   real(r8), pointer, dimension(:,:) :: vm_pert  ! perturbed zonal wind                         [m/s]
   real(r8), pointer, dimension(:,:) :: upwp_pert! perturbed meridional wind flux               [m^2/s^2]
   real(r8), pointer, dimension(:,:) :: vpwp_pert! perturbed zonal wind flux                    [m^2/s^2]

  end type clubb_mnts_2d_t

  !----------------------------------
  type clubb_to_host_t

    real(r8), pointer, dimension(:,:) ::  khzt    ! eddy diffusivity on thermo levels [m^2/s]
    real(r8), pointer, dimension(:,:) ::  khzm    ! eddy diffusivity on thermo levels [m^2/s]
    real(r8), pointer, dimension(:,:) ::  qclvar  ! cloud water variance [kg^2/kg^2]
    real(r8), pointer, dimension(:,:) ::  cloud_frac  ! cloud fraction [-]

    real(r8), pointer, dimension(:,:) ::  pdf_zm_w_1
    real(r8), pointer, dimension(:,:) ::  pdf_zm_w_2
    real(r8), pointer, dimension(:,:) ::  pdf_zm_varnce_w_1
    real(r8), pointer, dimension(:,:) ::  pdf_zm_varnce_w_2
    real(r8), pointer, dimension(:,:) ::  pdf_zm_mixt_frac

  end type clubb_to_host_t

end module clubb_intr_types
