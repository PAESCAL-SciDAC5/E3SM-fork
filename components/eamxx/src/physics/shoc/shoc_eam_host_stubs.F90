module shoc_eam_host_stubs
  !---------------------------------------------------------------------------
  ! Stand-ins for the EAM host modules that the SHOC+MF shoc.F90 depends on and
  ! that the EAMxx CMake build of shoc.F90 (libshoc, shoc_in_and_out driver,
  ! BFB tests) does not compile:
  !
  !   wv_saturation   -> qsat, no_ip_hltalt        (MF plume condensation)
  !   time_manager    -> get_nstep, is_first_step,
  !                      is_first_restart_step     (MF cold-pool bookkeeping)
  !   shr_RandNum_mod -> ShrKissRandGen            (Poisson entrainment draws)
  !
  ! shoc.F90 selects this module only under SCREAM_CONFIG_IS_CMAKE; the EAM
  ! build keeps using the real host modules.
  !
  ! Fidelity notes (numbers are EAM's, from shr_const_mod / physconst):
  !  * qsat: EAM looks es up in a 1 K table with linear interpolation
  !    (wv_saturation::estblf) built from the GoffGratch formulas. Here the
  !    GoffGratch formula is evaluated directly (wv_sat_methods::wv_sat_svp_trans
  !    with the water/ice blend over ttrice = 20 K), i.e. the same physics
  !    without the table interpolation error. This only affects the MF
  !    condensation (condensation_mf), never the SHOC PDF, so do_edmf=.false.
  !    stays bit-for-bit with standard SHOC.
  !  * get_nstep: the host step counter is advanced by the C bridge
  !    (shoc_iso_c::shoc_main_c) once per shoc_main call, so a driver that
  !    calls shoc_main once per host step sees nstep = 0,1,2,... exactly as
  !    EAM's time_manager would report.
  !  * ShrKissRandGen: same KISS generator (share/RandNum kissvec.c) and the
  !    same seed handling as shr_RandNum_mod, so the Poisson draws are the
  !    ones EAM would produce for the same seed.
  !---------------------------------------------------------------------------
  use physics_utils, only: rtype
  use kissvec_mod,   only: kissvec
  implicit none
  private

  public :: qsat, no_ip_hltalt
  public :: get_nstep, is_first_step, is_first_restart_step, shoc_host_set_nstep
  public :: ShrKissRandGen

  ! EAM constants (shr_const_mod.F90 / physconst.F90)
  real(rtype), parameter :: tmelt   = 273.15_rtype               ! freezing point of fresh water [K]
  real(rtype), parameter :: h2otrip = 273.16_rtype               ! triple point of fresh water [K]
  real(rtype), parameter :: tboil   = 373.16_rtype               ! boiling point of water at 1 atm [K]
  real(rtype), parameter :: latvap  = 2.501e6_rtype              ! latent heat of evaporation [J/kg]
  real(rtype), parameter :: mwwv    = 18.016_rtype               ! molecular weight of water vapor
  real(rtype), parameter :: mwdair  = 28.966_rtype               ! molecular weight of dry air
  real(rtype), parameter :: epsilo  = mwwv/mwdair                ! ratio of h2o to dry air molecular weights
  real(rtype), parameter :: omeps   = 1.0_rtype - epsilo
  real(rtype), parameter :: ttrice  = 20.00_rtype                ! water->ice es transition range [K]

  ! Host time-step counter (EAM time_manager stand-in)
  integer, save :: nstep_host = 0

  ! KISS random generator with the ShrKissRandGen API subset used by shoc.F90
  ! (constructor from an integer seed(:,4) array, %random(array), %finalize()).
  integer, parameter :: r8 = selected_real_kind(12)
  type ShrKissRandGen
    integer, allocatable, private :: seed(:,:)
  contains
    procedure, non_overridable :: random   => kiss_random
    procedure, non_overridable :: finalize => kiss_finalize
  end type ShrKissRandGen

  interface ShrKissRandGen
    module procedure ShrKissRandGen_constructor
  end interface ShrKissRandGen

contains

  !------------------------------------------------------------------
  ! wv_saturation stand-ins
  !------------------------------------------------------------------
  elemental function GoffGratch_svp_water(t) result(es)
    real(rtype), intent(in) :: t  ! Temperature in Kelvin
    real(rtype) :: es             ! SVP in Pa
    ! uncertain below -70 C  (wv_sat_methods::GoffGratch_svp_water)
    es = 10._rtype**(-7.90298_rtype*(tboil/t-1._rtype)+ &
         5.02808_rtype*log10(tboil/t)- &
         1.3816e-7_rtype*(10._rtype**(11.344_rtype*(1._rtype-t/tboil))-1._rtype)+ &
         8.1328e-3_rtype*(10._rtype**(-3.49149_rtype*(tboil/t-1._rtype))-1._rtype)+ &
         log10(1013.246_rtype))*100._rtype
  end function GoffGratch_svp_water

  elemental function GoffGratch_svp_ice(t) result(es)
    real(rtype), intent(in) :: t  ! Temperature in Kelvin
    real(rtype) :: es             ! SVP in Pa
    ! good down to -100 C  (wv_sat_methods::GoffGratch_svp_ice)
    es = 10._rtype**(-9.09718_rtype*(h2otrip/t-1._rtype)-3.56654_rtype* &
         log10(h2otrip/t)+0.876793_rtype*(1._rtype-t/h2otrip)+ &
         log10(6.1071_rtype))*100._rtype
  end function GoffGratch_svp_ice

  elemental function svp_trans(t) result(es)
    ! wv_sat_methods::wv_sat_svp_trans with the GoffGratch scheme (EAM default)
    real(rtype), intent(in) :: t
    real(rtype) :: es
    real(rtype) :: esice, weight
    if (t >= (tmelt - ttrice)) then
       es = GoffGratch_svp_water(t)
    else
       es = 0.0_rtype
    end if
    if (t < tmelt) then
       esice = GoffGratch_svp_ice(t)
       if ( (tmelt - t) > ttrice ) then
          weight = 1.0_rtype
       else
          weight = (tmelt - t)/ttrice
       end if
       es = weight*esice + (1.0_rtype - weight)*es
    end if
  end function svp_trans

  elemental subroutine qsat(t, p, es, qs)
    ! wv_saturation::qsat (mandatory arguments only): saturation vapor
    ! pressure and saturation specific humidity at (t, p).
    real(rtype), intent(in)  :: t   ! Temperature [K]
    real(rtype), intent(in)  :: p   ! Pressure [Pa]
    real(rtype), intent(out) :: es  ! Saturation vapor pressure [Pa]
    real(rtype), intent(out) :: qs  ! Saturation specific humidity [kg/kg]
    es = svp_trans(t)
    ! wv_sat_methods::wv_sat_svp_to_qsat
    if ( (p - es) <= 0._rtype ) then
       qs = 1.0_rtype
    else
       qs = epsilo*es / (p - omeps*es)
    end if
    ! Ensures returned es is consistent with limiters on qs.
    es = min(es, p)
  end subroutine qsat

  elemental subroutine no_ip_hltalt(t, hltalt)
    ! wv_saturation::no_ip_hltalt: latent heat of vaporization of pure liquid
    ! water at temperature t (constant slope -2369 J/(kg K) = cpv - cw above
    ! freezing).
    real(rtype), intent(in)  :: t
    real(rtype), intent(out) :: hltalt
    hltalt = latvap
    if (t >= tmelt) then
       hltalt = hltalt - 2369.0_rtype*(t-tmelt)
    end if
  end subroutine no_ip_hltalt

  !------------------------------------------------------------------
  ! time_manager stand-ins
  !------------------------------------------------------------------
  subroutine shoc_host_set_nstep(n)
    integer, intent(in) :: n
    nstep_host = n
  end subroutine shoc_host_set_nstep

  function get_nstep() result(n)
    integer :: n
    n = nstep_host
  end function get_nstep

  function is_first_step() result(f)
    logical :: f
    f = (nstep_host == 0)
  end function is_first_step

  function is_first_restart_step() result(f)
    logical :: f
    f = .false.   ! the standalone driver never restarts
  end function is_first_restart_step

  !------------------------------------------------------------------
  ! shr_RandNum_mod::ShrKissRandGen stand-in (identical algorithm/seeding)
  !------------------------------------------------------------------
  function ShrKissRandGen_constructor(seed) result(rand_gen)
    integer, intent(in) :: seed(:,:)
    type(ShrKissRandGen) :: rand_gen
    allocate(rand_gen%seed(size(seed, 1),4))
    rand_gen%seed = seed
  end function ShrKissRandGen_constructor

  subroutine kiss_random(self, array)
    class(ShrKissRandGen), intent(inout) :: self
    real(r8), dimension(:,:), intent(out) :: array
    integer :: nstream, i
    nstream = size(self%seed, 1)
    do i = 1, size(array, 2)
       call kissvec(self%seed(:,1), self%seed(:,2), self%seed(:,3), &
            self%seed(:,4), array(:,i), nstream)
    end do
  end subroutine kiss_random

  subroutine kiss_finalize(self)
    class(ShrKissRandGen), intent(inout) :: self
    if ( allocated(self%seed) ) deallocate(self%seed)
  end subroutine kiss_finalize

end module shoc_eam_host_stubs
