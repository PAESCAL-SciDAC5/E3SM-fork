module lscale_mod

    implicit none

    private
    public :: calculate_lscale
    public :: lscale_init

    contains

subroutine lscale_init()
   !-----------------------------------------------------------------------
   ! Purpose:
   ! Initialization related to the calculation of turbulent length scales
   ! in the physics driver.
   !-----------------------------------------------------------------------
   use cam_history,      only: addfld, add_default
   
   ! Add output variables

   call addfld ('LSCALE',     (/ 'lev' /), 'A', 'm', 'Mixing Length Scale')
   call addfld ('LSCALE_UP',  (/ 'lev' /), 'A', 'm', 'Upward Mixing Length Scale')
   call addfld ('LSCALE_DOWN',(/ 'lev' /), 'A', 'm', 'Downward Mixing Length Scale')

   call add_default('LSCALE',     1, ' ')
   call add_default('LSCALE_UP',  1, ' ')
   call add_default('LSCALE_DOWN',1, ' ')

end subroutine lscale_init

subroutine calculate_lscale(ncol, pcols, pver, pverp,       &
                            temp_in, pmid_in, qv_in, qc_in, &
                            zm_in,   zi_in,   tke_in,       &
                            lscale_out, lscale_up_out, lscale_down_out )
    !------------------------------------------------------------------------------ 
    ! Purpose: 
    ! Calculate mixing length scales using CLUBB's compute_mixing_length subroutine
    ! and pass the results back to the calling subroutine for diagnostic output
    ! and future AMR work.
    !------------------------------------------------------------------------------
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use physconst,          only: gravit, rair, cpair, latvap, zvir, epsilo

    implicit none

    !--------------------
    ! Arguments
    !--------------------
    integer, intent(in) :: ncol, pcols, pver, pverp     ! pverp = pver + 1

    ! Meteorological fields needed to diagnose Lscale.

    real(r8), intent(in)  :: temp_in(pcols,pver)   ! temperature [K] at layer midpoints
    real(r8), intent(in)  :: pmid_in(pcols,pver)   ! pressure [Pa] at layer midpoints
    real(r8), intent(in)  ::   qv_in(pcols,pver)   ! water vapor mixing ratio [kg/kg] at layer midpoints
    real(r8), intent(in)  ::   qc_in(pcols,pver)   ! cloud liquid mixing ratio [kg/kg] at layer midpoints
    real(r8), intent(in)  ::   zm_in(pcols,pver)   ! geopotential height [m] at layer midpoints
    real(r8), intent(in)  ::   zi_in(pcols,pverp)  ! geopotential height [m] at layer interfaces
    real(r8), intent(in)  ::  tke_in(pcols,pverp)  ! turbulent kinetic energy [m2/s2] at layer interfaces

    ! Lscale and its upward and downward components, defined at layer midpoints

    real(r8), intent(out) :: lscale_out     (pcols,pver)   ! turbulent mixing length scale [m]
    real(r8), intent(out) :: lscale_up_out  (pcols,pver)   ! turbulent mixing length scale, upward [m]
    real(r8), intent(out) :: lscale_down_out(pcols,pver)   ! turbulent mixing length scale, downward [m]

    !-----------------------------------------------------------------------------------------
    ! Local arrays. Their shape and direction of indexing follow CLUBB:
    !  - Variables defined on thermodynamic levels (corresponding to E3SM's layer midpoints)
    !    have a extra ghost level below ground/Earth surface;
    !  - Level indices start from the surface and increase upward.
    !-----------------------------------------------------------------------------------------
    ! Input arrays of CLUBB's compute_mixing_length subroutine

    real(r8) :: zt       (pverp)  ! height at thermodynamic levels
    real(r8) :: zm       (pverp)  ! height at momentum levels
    real(r8) :: dzm      (pverp)  ! grid spacing at momentum levels
    real(r8) :: invrs_dzm(pverp)  ! 1/dzm

    real(r8) :: thvm     (pverp)  ! virtual potential temperature
    real(r8) :: thlm     (pverp)  ! liquid water potential temperature
    real(r8) :: rtm      (pverp)  ! total water mixing ratio
    real(r8) :: em       (pverp)  ! turbulent kinetic energy
    real(r8) :: p_in_Pa  (pverp)  ! pressure in Pa
    real(r8) :: exner    (pverp)  ! exner function
    real(r8) :: thv_ds   (pverp)  ! dry static virtual potential temperature

    ! Output arrays of CLUBB's compute_mixing_length subroutine.
    ! These are defined at layer midpoints (CLUBB's thermodynamic levels)
    ! and have a "ghost" level below Earth's surface.

    real(r8) :: lscale_tmp     (pverp)
    real(r8) :: lscale_up_tmp  (pverp) 
    real(r8) :: lscale_down_tmp(pverp)

    !--------------------
    ! Local scalars
    !--------------------
    integer :: ii, kk, kflip    ! array indices

    real(r8) :: temp_k
    real(r8) :: theta_v
    real(r8) :: press_pa
    real(r8) :: qv, qc
    real(r8) :: exner_clubb

    real(r8), parameter :: p0_clubb = 100000._r8      ! reference pressure (Pa)

    !--------------------------------------
    ! Process each column separately
    !--------------------------------------
    do ii = 1, ncol

       !====================================================
       ! 1. Prepare input variables for this column
       !    following clubb_tend_cam in clubb_intr.F90.
       !====================================================
       do kk = 1, pver

          !------------------------------------
          ! 1.1 Thermodynamic-level variables
          !------------------------------------
          kflip = pver - kk + 1   ! for flipping to CLUBB ordering

          temp_k = temp_in(ii, kflip)
          press_pa = pmid_in(ii, kflip)
          exner_clubb = (p0_clubb/press_pa)**(rair/cpair)

          qv = qv_in(ii, kflip)  ! water vapor
          qc = qc_in(ii, kflip)  ! cloud liquid

          ! virtual potential temperature per CLUBB
          theta_v = temp_k*exner_clubb*(1.0_r8 + zvir*qv - qc)

          rtm(kk+1) = qv + qc ! total water mixing ratio per CLUBB
          thlm(kk+1) = (temp_k - (latvap/cpair)*qc)*exner_clubb
          p_in_Pa(kk+1) = press_pa
          exner(kk+1) = 1.0_r8 / exner_clubb

          ! note here thv_ds is assigned actual moist theta_v 
          ! following what is done in clubb_tend_cam in clubb_intr.F90
          thv_ds(kk+1) = theta_v

          ! thvm is then calculated again in advance_clubb_core
          ! as follows:
          thvm(kk+1) = thlm(kk+1) + ((1.0_r8 - epsilo)/epsilo) * thv_ds(kk+1) *rtm(kk+1) &
                        + (latvap/(cpair*exner(kk+1)) - (1.0_r8/epsilo) *thv_ds(kk+1))*qc

          ! thermodynamic-level heights
          zt(kk+1) = zm_in(ii, kflip) - zi_in(ii, pver+1)

          !-------------------------------
          ! 1.2 Momentum-level variables
          !-------------------------------
          em(kk) = tke_in(ii,pverp-kk+1)                    ! TKE
          zm(kk) =  zi_in(ii,pverp-kk+1) - zi_in(ii,pver+1) ! geopot. height above sfc (of layer interfaces) 

       end do

       !------------------------------------------------------------------------------------
       ! 1.3 For variables on zt levels (layer midpoints), fill in values for ghost level.
       !------------------------------------------------------------------------------------
            zt(1) = -1.0_r8 * zt(2)
         exner(1) =   exner(2)
       p_in_Pa(1) = p_in_Pa(2)
           rtm(1) =     rtm(2)
          thlm(1) =    thlm(2)
        thv_ds(1) =  thv_ds(2)
          thvm(1) =    thvm(2)

       !------------------------------------------------------------------------------------
       ! 1.5 For variables on zm levels (layer interfaces), set values at model top.
       !------------------------------------------------------------------------------------
       em(pverp) = tke_in(ii,1)
       zm(pverp) =  zi_in(ii,1) - zi_in(ii,pver+1)

       !------------------------------------------------------------------------------------
       ! 1.6 Calculate dz and 1/dz.
       !------------------------------------------------------------------------------------
       do kk = 1, pver
               dzm(kk) = zt(kk+1) - zt(kk)
         invrs_dzm(kk) = 1.0_r8 / dzm(kk)
       end do
             dzm(pverp) =       dzm(pver)  ! does not matter much what we put here
       invrs_dzm(pverp) = invrs_dzm(pver)

       !========================================================
       ! 2. Call CLUBB's compute_mixing_length for this column
       !========================================================
       call compute_mixing_length_standalone( &
         pverp, zt, zm, dzm, invrs_dzm, & 
         thvm, thlm, rtm, em, p_in_Pa, exner, thv_ds, &
         lscale_tmp, lscale_up_tmp, lscale_down_tmp )

       !========================================================
       ! 3. Flip the vertical indexing for output to host model;
       !    ignore the ghost level below Earth's surface.
       !========================================================
       do kk = 1, pver
          lscale_out     (ii,kk) = lscale_tmp     (pverp-kk+1)
          lscale_up_out  (ii,kk) = lscale_up_tmp  (pverp-kk+1)
          lscale_down_out(ii,kk) = lscale_down_tmp(pverp-kk+1)
       end do

    end do ! ii = 1,ncol

end subroutine calculate_lscale

subroutine compute_mixing_length_standalone( &
    nzmax, zt, zm, dzm, invrs_dzm, &
    thvm, thlm, rtm, em, p_in_Pa, exner, thv_ds, &
    Lscale, Lscale_up, Lscale_down )

    ! Taken from CLUBB, modified to be standalone - H. Xiao
    ! Description:
    !   Larson's 5th moist, nonlocal length scale
    ! 
    ! References:
    !   Section 3b ( /Eddy length formulation/ ) of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    ! 
    ! Notes: 
    ! 
    !   The equation for the rate of change of theta_l and r_t of the parcel with
    !   respect to height, due to entrainment, is:
    ! 
    !           d(thl_par)/dz = - mu * ( thl_parcel - thl_environment );
    ! 
    !           d(rt_par)/dz = - mu * ( rt_parcel - rt_environment );
    ! 
    !   where mu is the entrainment rate,
    !   such that:
    ! 
    !           mu = (1/m)*(dm/dz);
    ! 
    !   where m is the mass of the parcel.  The value of mu is set to be a
    !   constant.
    ! 
    !   The differential equations are solved for given the boundary condition 
    !   and given the fact that the value of thl_environment and rt_environment 
    !   are treated as changing linearly for a parcel of air from one grid level
    !   to the next.
    ! 
    !   For the special case where entrainment rate, mu, is set to 0,
    !   thl_parcel and rt_parcel remain constant
    ! 
    ! 
    !   The equation for Lscale_up is:
    ! 
    !       INT(z_i:z_i+Lscale_up) g * ( thv_par - thvm ) / thvm dz = -em(z_i);
    ! 
    !   and for Lscale_down
    !   
    !       INT(z_i-Lscale_down:z_i) g * ( thv_par - thvm ) / thvm dz = em(z_i);
    ! 
    !   where thv_par is theta_v of the parcel, thvm is the mean
    !   environmental value of theta_v, z_i is the altitude that the parcel
    !   started from, and em is the mean value of TKE at
    !   altitude z_i (which gives the parcel its initial boost)
    ! 
    !   The increment of CAPE (convective air potential energy) for any two 
    !   successive vertical levels is:
    !   
    !       Upwards:
    !           CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
    !
    !       Downwards: 
    !           CAPE_incr = INT(z_(-1):z_0) g * ( thv_par - thvm ) / thvm dz
    ! 
    !   Thus, the derivative of CAPE with respect to height is:
    ! 
    !           dCAPE/dz = g * ( thv_par - thvm ) / thvm.
    ! 
    !   A purely trapezoidal rule is used between levels, and is considered
    !   to vary linearly at all altitudes.  Thus, dCAPE/dz is considered to be 
    !   of the form:  A * (z-zo) + dCAPE/dz|_(z_0),
    !   where A = ( dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) ) / ( z_1 - z_0 )
    ! 
    !   The integral is evaluated to find the CAPE increment between two
    !   successive vertical levels.  The result either adds to or depletes
    !   from the total amount of energy that keeps the parcel ascending/descending.
    ! 
    ! 
    ! IMPORTANT NOTE:
    !   This subroutine has been optimized by adding precalculations, rearranging
    !   equations to avoid divides, and modifying the algorithm entirely.
    !       -Gunther Huebler
    ! 
    !   The algorithm previously used looped over every grid level, following a
    !   a parcel up from its initial grid level to its max. The very nature of this 
    !   algorithm is an N^2 
    !--------------------------------------------------------------------------------

    ! mu = (1/M) dM/dz > 0.  mu=0 for no entrainment.
    ! Siebesma recommends mu=2e-3, although most schemes use mu=1e-4
    ! When mu was fixed, we used the value mu = 6.e-4

    use shr_kind_mod,  only: r8 => shr_kind_r8
    use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
                             shr_const_latice, shr_const_latsub, shr_const_rgas, &
                             shr_const_mwwv, shr_const_stebol, shr_const_tkfrz, &
                             shr_const_mwdair, shr_const_g, shr_const_karman, &
                             shr_const_rhofw
    
    implicit none

    ! constants
    ! copied from constants_clubb module
    real(r8), parameter :: &
        Cp = shr_const_cpdair, & ! Dry air specific heat at constant p [J/kg/K]
        Lv = shr_const_latvap, & ! Latent heat of vaporization         [J/kg]
        Rd = shr_const_rdair, & ! Dry air gas constant                [J/kg/K]
        ep  = shr_const_mwwv/shr_const_mwdair, & ! ep  = 0.622  [-]
        ep1 = (1.0_r8-ep)/ep, & ! ep1 = 0.61   [-]
        ep2 = 1.0_r8/ep, & ! ep2 = 1.61   [-]
        grav = shr_const_g, & ! Gravitational acceleration     [m/s^2]
        zero_threshold = 0.0_r8, &
        eps = 1.0e-10_r8  ! small number to avoid division by zero

    ! External
    intrinsic :: min, max, sqrt

    ! Almost-Constant Parameters
    real(r8), parameter ::  & 
      ! don't know why there are two mininum Lscale values
      ! but both are used later in the code, so kept as it is.
      zlmin = 0.1_r8, & ! Minimum value for Lscale [m]
      Lscale_sfclyr_depth = 500._r8, & ! [m]
      em_min = 1.5_r8*(2.0e-2_r8)**2, & ! Minimum TKE [m^2/s^2]
      ! this is done following setup_parameters in parameters_tunable.F90.
      lmin = 40.0_r8 * 0.1_r8 ! Minimum value for Lscale                         [m]

    ! tunable parameters
    real(r8), parameter ::  & 
      ! Fractional extrainment rate per unit altitude [1/m] 
      ! following CLUBB (in parameters_tunable.F90)
    !   mu = 1.0e-3_r8, &
      mu = 5.0e-4_r8, &
      ! Maximum allowable value for Lscale [m]
      ! CLUBB sets it to a quarter of horizontal grid spacing
      ! assuming a 25 km grid spacing for now
      Lscale_max = 0.25_r8 * 25000.0_r8

    ! Input Variables
    integer, intent(in) :: nzmax
    real(r8), dimension(nzmax), intent(in) ::  & 
      thvm,    & ! Virtual potential temp. on themodynamic level  [K]
      thlm,    & ! Liquid potential temp. on themodynamic level   [K]
      rtm,     & ! Total water mixing ratio on themodynamic level [kg/kg]
      em,      & ! em = 3/2 * w'^2; on momentum level             [m^2/s^2]
      exner,   & ! Exner function on thermodynamic level          [-]
      p_in_Pa, & ! Pressure on thermodynamic level                [Pa]
      ! Note:  thv_ds used as a reference theta_l here
      ! (don't know what this means. - H. Xiao)
      thv_ds     ! Dry, base-state theta_v on thermodynamic level [K]

    real(r8), dimension(nzmax), intent(in) ::  & 
      zt,        & ! Height on thermodynamic level                 [m]
      zm,        & ! Height on momentum level                    [m]
      dzm,        & ! Grid spacing on momentum level               [m]
      invrs_dzm   ! 1/dzm                                      [1/m]

    real(r8), dimension(nzmax), intent(out) ::  & 
      Lscale,    & ! Mixing length      [m]
      Lscale_up, & ! Mixing length up   [m]
      Lscale_down  ! Mixing length down [m]

    ! Local Variables

    integer :: i, j

    real(r8) :: tke, CAPE_incr

    real(r8) :: dCAPE_dz_j, dCAPE_dz_j_minus_1, dCAPE_dz_j_plus_1

    ! Temporary arrays to store calculations to speed runtime
    real(r8), dimension(nzmax) :: &
        exp_mu_dzm, &
        invrs_dzm_on_mu, &
        grav_on_thvm, &
        thl_par_j_precalc, &
        rt_par_j_precalc, &
        tke_i, &
        thl_par_1, &
        rt_par_1, &
        tl_par_1, &
        rsatl_par_1, &
        s_par_1, &
        rc_par_1, &
        thv_par_1, &
        dCAPE_dz_1, &
        CAPE_incr_1, &
        Lv_coef, &
        entrain_coef

    ! Minimum value for Lscale that will taper off with height
    real(r8) :: lminh

    ! Parcel quantities at grid level j
    real(r8) :: thl_par_j, rt_par_j, rc_par_j, thv_par_j

    ! Used in latent heating calculation
    real(r8) :: tl_par_j, rsatl_par_j, s_par_j

    ! Variables to make L nonlocal
    real(r8) :: Lscale_up_max_alt, Lscale_down_min_alt

    ! Variables used to precalculate values
    real(r8) :: &
        Lv2_coef, &
        tl_par_j_sqd, &
        invrs_dCAPE_diff, &
        invrs_Lscale_sfclyr_depth

    ! ---- Begin Code ----

    !---------- Mixing length computation ----------------------------------

   !  if( abs(mu) < eps ) then
   !      write(fstderr,*) "Entrainment rate mu cannot be 0"
   !      stop "Fatal error in subroutine compute_mixing_length"
   !  end if

    ! Initialize arrays and precalculate values for computational efficiency
    do i = 1, nzmax

        ! Initialize up and down arrays
        Lscale_up(i) = zlmin
        Lscale_down(i) = zlmin

        ! Precalculate values to avoid unnecessary calculations later
        exp_mu_dzm(i) = exp( -mu * dzm(i) )
        invrs_dzm_on_mu(i) = ( invrs_dzm(i) ) / mu
        grav_on_thvm(i) = grav / thvm(i)
        Lv_coef(i) = Lv / ( exner(i) * cp ) - ep2 * thv_ds(i)
        entrain_coef(i) = ( 1.0_r8 - exp_mu_dzm(i) ) * invrs_dzm_on_mu(i)

    end do

    ! Avoid uninitialized memory (these values are not used in Lscale)
    Lscale_up(1)   = 0.0_r8
    Lscale_down(1) = 0.0_r8

    ! Precalculations of single values to avoid unnecessary calculations later
    Lv2_coef = ep * Lv**2 / ( Rd * cp )
    invrs_Lscale_sfclyr_depth = 1.0_r8 / Lscale_sfclyr_depth

    ! Calculate initial turbulent kinetic energy for zt levels
    ! using linear interpolation following CLUBB
    do i = 2, nzmax
      tke_i (i) = em(i) * ( zt(i) - zm(i-1) ) + &
                  em(i-1) * ( zm(i) - zt(i) )  
      tke_i (i) = max( tke_i(i) / (zm(i) - zm(i-1)), em_min )
    end do
      tke_i (1) = (zt(1) - zm(1)) * (em(2) - em(1)) / (zm(2) - zm(1)) + em(1) 
      tke_i (1) = max( tke_i(1) , em_min )
    

    ! ---------------- Upwards Length Scale Calculation ----------------

    ! Precalculate values for upward Lscale, these are useful only if a parcel can rise
    ! more than one level. They are used in the equations that calculate thl and rt 
    ! recursively for a parcel as it ascends
    do j = 2, nzmax-1

        thl_par_j_precalc(j) = thlm(j) - thlm(j-1) * exp_mu_dzm(j-1)  &
                               - ( thlm(j) - thlm(j-1) ) * entrain_coef(j-1)

        rt_par_j_precalc(j) = rtm(j) - rtm(j-1) * exp_mu_dzm(j-1)  &
                              - ( rtm(j) - rtm(j-1) ) * entrain_coef(j-1)
    end do


    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain 
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations. 

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    do j = 3, nzmax

        thl_par_1(j) = thlm(j) - ( thlm(j) - thlm(j-1) ) * entrain_coef(j-1)

        tl_par_1(j) = thl_par_1(j) * exner(j)

        rt_par_1(j) = rtm(j) - ( rtm(j) - rtm(j-1) ) * entrain_coef(j-1)
       
    end do


    ! Caclculate initial rsatl for parcels at each grid level, this function is elemental
    rsatl_par_1(3:) = sat_mixrat_liq_clubb( p_in_Pa(3:), tl_par_1(3:) )


    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    do j = 3, nzmax

        tl_par_j_sqd = tl_par_1(j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )  
        ! 
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )  
        s_par_1(j) = ( rt_par_1(j) - rsatl_par_1(j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(j) )

        rc_par_1(j) = max( s_par_1(j), zero_threshold )

        ! theta_v of entraining parcel at grid level j 
        thv_par_1(j) = thl_par_1(j) + ep1 * thv_ds(j) * rt_par_1(j) + Lv_coef(j) * rc_par_1(j)

        
        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(j) = grav_on_thvm(j) * ( thv_par_1(j) - thvm(j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(j) = 0.5_r8 * dCAPE_dz_1(j) * dzm(j-1)

    end do

    
    ! Calculate Lscale_up for each grid level. If the TKE from a parcel has not been completely
    ! exhausted by the initial change then continue the exhaustion calculations here for a single
    ! grid level at a time until the TKE is exhausted. 

    Lscale_up_max_alt = 0.0_r8    ! Set initial max value for Lscale_up to 0
    do i = 2, nzmax-2

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i) + CAPE_incr_1(i+1) > 0.0_r8 ) then

            ! Calculate new TKE for parcel
            tke = tke_i(i) + CAPE_incr_1(i+1)

            ! Set j to 2 levels above current Lscale_up level, this is because we've already
            ! determined that the parcel can rise at least 1 full level
            j = i + 2

            ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
            thl_par_j = thl_par_1(i+1)
            rt_par_j  = rt_par_1(i+1)
            dCAPE_dz_j_minus_1 = dCAPE_dz_1(i+1)

            ! Continue change in TKE calculations until it is exhausted or the max grid
            ! level has been reached. j is the next grid level above the level that can
            ! be reached for a parcel starting at level i. If TKE is exhausted in this loop
            ! that means the parcel starting at i cannot reach level j, but has reached j-1
            do while ( j < nzmax )
           
                ! thl, rt of parcel are conserved except for entrainment
                ! 
                ! The values of thl_env and rt_env are treated as changing linearly for a parcel 
                ! of air ascending from level j-1 to level j

                ! theta_l of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(thl_par)/dz = - mu * ( thl_par - thl_env )
                thl_par_j = thl_par_j_precalc(j) + thl_par_j * exp_mu_dzm(j-1)


                ! r_t of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(rt_par)/dz = - mu * ( rt_par - rt_env )
                rt_par_j = rt_par_j_precalc(j) + rt_par_j * exp_mu_dzm(j-1)
                

                ! Include effects of latent heating on Lscale_up 6/12/00
                ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
                ! Probably should use properties of bump 1 in Gaussian, not mean!!!

                tl_par_j = thl_par_j*exner(j)

                rsatl_par_j = sat_mixrat_liq_clubb( p_in_Pa(j), tl_par_j )

                tl_par_j_sqd = tl_par_j**2

                ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
                !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
                ! and SD's beta (eqn. 8),
                !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
                ! 
                ! Simplified by multiplying top and bottom by tl^2 to avoid a 
                ! divide and precalculating ep * Lv**2 / ( Rd * cp )  
                s_par_j = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
                          / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

                rc_par_j = max( s_par_j, zero_threshold )

                ! theta_v of entraining parcel at grid level j 
                thv_par_j = thl_par_j + ep1 * thv_ds(j) * rt_par_j  &
                            + Lv_coef(j) * rc_par_j

                ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
                dCAPE_dz_j = grav_on_thvm(j) * ( thv_par_j - thvm(j) )

                ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
                ! Trapezoidal estimate between grid levels j and j-1
                CAPE_incr = 0.5_r8 * ( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * dzm(j-1)

                ! Exit loop early if tke has been exhaused between level j and j+1
                if ( tke + CAPE_incr <= 0.0_r8 ) then
                    exit
                end if

                ! Save previous dCAPE value for next cycle
                dCAPE_dz_j_minus_1 = dCAPE_dz_j

                ! Caclulate new TKE and increment j
                tke = tke + CAPE_incr
                j = j + 1

            enddo


            ! Add full grid level thickness for each grid level that was passed without the TKE
            ! being exhausted, difference between starting level (i) and last level passed (j-1)
            Lscale_up(i) = Lscale_up(i) + zt(j-1) - zt(i)


            if ( j < nzmax ) then

                ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
                ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.
        
                if ( abs( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) * 2 <= &
                     abs( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * eps ) then

                    ! Special case where dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) = 0
                    ! Find the remaining distance z - z_0 that it takes to 
                    ! exhaust the remaining TKE

                    Lscale_up(i) = Lscale_up(i) + ( - tke / dCAPE_dz_j )

                else

                    ! Case used for most scenarios where dCAPE/dz|_(z_1) /= dCAPE/dz|_(z_0)
                    ! Find the remaining distance z - z_0 that it takes to exhaust the
                    ! remaining TKE (tke_i), using the quadratic formula (only the
                    ! negative (-) root works in this scenario).
                    invrs_dCAPE_diff = 1.0_r8 / ( dCAPE_dz_j - dCAPE_dz_j_minus_1 )

                    Lscale_up(i) = Lscale_up(i) &
                                   - dCAPE_dz_j_minus_1 * invrs_dCAPE_diff * dzm(j-1)  &
                                   - sqrt( dCAPE_dz_j_minus_1**2 &
                                            - 2.0_r8 * tke * invrs_dzm(j-1) & 
                                              * ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) ) &
                                     * invrs_dCAPE_diff  * dzm(j-1)
                endif

            end if

        else    ! TKE for parcel at level (i) was exhaused before one full grid level
            
            ! Find the remaining distance z - z_0 that it takes to exhaust the
            ! remaining TKE (tke_i), using the quadratic formula. Simplified 
            ! since dCAPE_dz_j_minus_1 = 0.0
            Lscale_up(i) = Lscale_up(i) - sqrt( - 2.0_r8 * tke_i(i) &
                                                  * dzm(i) * dCAPE_dz_1(i+1) ) &
                                          / dCAPE_dz_1(i+1)  
        endif


        ! If a parcel at a previous grid level can rise past the parcel at this grid level
        ! then this one should also be able to rise up to that height. This feature insures 
        ! that the profile of Lscale_up will be smooth, thus reducing numerical instability.
        if ( zt(i) + Lscale_up(i) < Lscale_up_max_alt ) then

            ! A lower starting parcel can ascend higher than this one, set height to the max
            ! that any lower starting parcel can ascend to
            Lscale_up(i) = Lscale_up_max_alt - zt(i)
        else 
            
            ! This parcel can ascend higher than any below it, save final height
            Lscale_up_max_alt = Lscale_up(i) + zt(i)
        end if


    end do


    ! ---------------- Downwards Length Scale Calculation ----------------

    ! Precalculate values for downward Lscale, these are useful only if a parcel can descend
    ! more than one level. They are used in the equations that calculate thl and rt 
    ! recursively for a parcel as it descends
    do j = 2, nzmax-1

        thl_par_j_precalc(j) = thlm(j) - thlm(j+1) * exp_mu_dzm(j)  &
                               - ( thlm(j) - thlm(j+1) ) * entrain_coef(j)

        rt_par_j_precalc(j) = rtm(j) - rtm(j+1) * exp_mu_dzm(j)  &
                              - ( rtm(j) - rtm(j+1) ) * entrain_coef(j)
    end do


    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain 
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations. 

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    do j = 2, nzmax-1

        thl_par_1(j) = thlm(j) - ( thlm(j) - thlm(j+1) )  * entrain_coef(j)

        tl_par_1(j) = thl_par_1(j) * exner(j)

        rt_par_1(j) = rtm(j) - ( rtm(j) - rtm(j+1) ) * entrain_coef(j)
       
    end do


    ! Caclculate initial rsatl for parcels at each grid level, this function is elemental
    rsatl_par_1(2:) = sat_mixrat_liq_clubb( p_in_Pa(2:), tl_par_1(2:) )


    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    do j = 2, nzmax-1

        tl_par_j_sqd = tl_par_1(j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )  
        ! 
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )  
        s_par_1(j) = ( rt_par_1(j) - rsatl_par_1(j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(j) )

        rc_par_1(j) = max( s_par_1(j), zero_threshold )

        ! theta_v of entraining parcel at grid level j 
        thv_par_1(j) = thl_par_1(j) + ep1 * thv_ds(j) * rt_par_1(j) + Lv_coef(j) * rc_par_1(j)

        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(j) = grav_on_thvm(j) * ( thv_par_1(j) - thvm(j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(j) = 0.5_r8 * dCAPE_dz_1(j) * dzm(j)

    end do


    ! Calculate Lscale_down for each grid level. If the TKE from a parcel has not been completely
    ! exhausted by the initial change then continue the exhaustion calculations here for a single
    ! grid level at a time until the TKE is exhausted. 

    Lscale_down_min_alt = zt(nzmax)  ! Set initial min value for Lscale_down to max zt
    do i = nzmax, 3, -1

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i) - CAPE_incr_1(i-1) > 0.0_r8 ) then

            ! Calculate new TKE for parcel
            tke = tke_i(i) - CAPE_incr_1(i-1)

            ! Set j to 2 levels below current Lscale_down level, this is because we've already
            ! determined that the parcel can descend at least 1 full level
            j = i - 2

            ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
            thl_par_j = thl_par_1(i-1)
            rt_par_j = rt_par_1(i-1)
            dCAPE_dz_j_plus_1 = dCAPE_dz_1(i-1)


            ! Continue change in TKE calculations until it is exhausted or the min grid
            ! level has been reached. j is the next grid level below the level that can
            ! be reached for a parcel starting at level i. If TKE is exhausted in this loop
            ! that means the parcel starting at i cannot sink to level j, but can sink to j+1
            do while ( j >= 2 )

                ! thl, rt of parcel are conserved except for entrainment
                ! 
                ! The values of thl_env and rt_env are treated as changing linearly for a parcel 
                ! of air descending from level j to level j-1

                ! theta_l of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(thl_par)/dz = - mu * ( thl_par - thl_env )
                thl_par_j = thl_par_j_precalc(j) + thl_par_j * exp_mu_dzm(j)


                ! r_t of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(rt_par)/dz = - mu * ( rt_par - rt_env )
                rt_par_j = rt_par_j_precalc(j) + rt_par_j * exp_mu_dzm(j)


                ! Include effects of latent heating on Lscale_up 6/12/00
                ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
                ! Probably should use properties of bump 1 in Gaussian, not mean!!!

                tl_par_j = thl_par_j*exner(j)

                rsatl_par_j = sat_mixrat_liq_clubb( p_in_Pa(j), tl_par_j )

                tl_par_j_sqd = tl_par_j**2

                ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
                !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
                ! and SD's beta (eqn. 8),
                !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
                ! 
                ! Simplified by multiplying top and bottom by tl^2 to avoid a 
                ! divide and precalculating ep * Lv**2 / ( Rd * cp )  
                s_par_j = (rt_par_j - rsatl_par_j) * tl_par_j_sqd &
                          / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

                rc_par_j = max( s_par_j, zero_threshold )

                ! theta_v of entraining parcel at grid level j 
                thv_par_j = thl_par_j + ep1 * thv_ds(j) * rt_par_j + Lv_coef(j) * rc_par_j

                ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
                dCAPE_dz_j = grav_on_thvm(j) * ( thv_par_j - thvm(j) )

                ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
                ! Trapezoidal estimate between grid levels j+1 and j
                CAPE_incr = 0.5_r8 * ( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * dzm(j)

                ! Exit loop early if tke has been exhaused between level j+1 and j
                if ( tke - CAPE_incr <= 0.0_r8 ) then
                    exit
                endif

                ! Save previous dCAPE value for next cycle
                dCAPE_dz_j_plus_1 = dCAPE_dz_j

                ! Caclulate new TKE and increment j
                tke = tke - CAPE_incr
                j = j - 1

            enddo

            ! Add full grid level thickness for each grid level that was passed without the TKE
            ! being exhausted, difference between starting level (i) and last level passed (j+1)
            Lscale_down(i) = Lscale_down(i) + zt(i) - zt(j+1)

            
            if ( j >= 2 ) then
        
                ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
                ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.

                if ( abs( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) * 2 <= &
                     abs( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * eps ) then

                    ! Special case where dCAPE/dz|_(z_(-1)) - dCAPE/dz|_(z_0) = 0
                    ! Find the remaining distance z_0 - z that it takes to 
                    ! exhaust the remaining TKE

                    Lscale_down(i) = Lscale_down(i) + ( tke / dCAPE_dz_j )

                else

                    ! Case used for most scenarios where dCAPE/dz|_(z_(-1)) /= dCAPE/dz|_(z_0)
                    ! Find the remaining distance z_0 - z that it takes to exhaust the
                    ! remaining TKE (tke_i), using the quadratic formula (only the
                    ! negative (-) root works in this scenario) -- however, the
                    ! negative (-) root is divided by another negative (-) factor,
                    ! which results in an overall plus (+) sign in front of the
                    ! square root term in the equation below).
                    invrs_dCAPE_diff = 1.0_r8 / ( dCAPE_dz_j - dCAPE_dz_j_plus_1 )

                    Lscale_down(i) = Lscale_down(i) &
                                     - dCAPE_dz_j_plus_1 * invrs_dCAPE_diff * dzm(j)  &
                                     + sqrt( dCAPE_dz_j_plus_1**2 &
                                             + 2.0_r8 * tke * invrs_dzm(j)  &
                                               * ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) )  &
                                       * invrs_dCAPE_diff * dzm(j)
                endif

            end if

        else    ! TKE for parcel at level (i) was exhaused before one full grid level

            ! Find the remaining distance z_0 - z that it takes to exhaust the
            ! remaining TKE (tke_i), using the quadratic formula. Simplified 
            ! since dCAPE_dz_j_plus_1 = 0.0
            Lscale_down(i) = Lscale_down(i) + sqrt( 2.0_r8 * tke_i(i) &
                                                    * dzm(i-1) * dCAPE_dz_1(i-1) ) &
                                              / dCAPE_dz_1(i-1) 
        endif

        ! If a parcel at a previous grid level can descend past the parcel at this grid level
        ! then this one should also be able to descend down to that height. This feature insures
        ! that the profile of Lscale_down will be smooth, thus reducing numerical instability.
        if ( zt(i) - Lscale_down(i) > Lscale_down_min_alt ) then
            Lscale_down(i) = zt(i) - Lscale_down_min_alt
        else 
            Lscale_down_min_alt = zt(i) - Lscale_down(i)
        end if

    end do


    ! ---------------- Final Lscale Calculation ----------------

    do i = 2, nzmax, 1

        ! Make lminh a linear function starting at value lmin at the bottom
        ! and going to zero at 500 meters in altitude.
        ! Within a host model, increase mixing length in 500 m layer above *ground*
        lminh = max( zero_threshold, Lscale_sfclyr_depth - ( zt(i) - zm(1) ) ) &
                    * lmin * invrs_Lscale_sfclyr_depth

        Lscale_up(i)    = max( lminh, Lscale_up(i) )
        Lscale_down(i)  = max( lminh, Lscale_down(i) )

        ! When L is large, turbulence is weakly damped
        ! When L is small, turbulence is strongly damped
        ! Use a geometric mean to determine final Lscale so that L tends to become small 
        ! if either Lscale_up or Lscale_down becomes small. 
        Lscale(i) = sqrt( Lscale_up(i) * Lscale_down(i) )

    enddo

    ! Set the value of Lscale at the upper and lower boundaries.
    Lscale(1) = Lscale(2)
    Lscale(nzmax) = Lscale(nzmax-1)

    ! Vince Larson limited Lscale to allow host
    ! model to take over deep convection.  13 Feb 2008.
    Lscale = min( Lscale, Lscale_max )

 
   ! Ensure that no Lscale values are NaN
   !  if ( clubb_at_least_debug_level( 1 ) ) then

   !      call length_check( Lscale, Lscale_up, Lscale_down )

   !      if ( err_code == clubb_fatal_error ) then

   !        write(fstderr,*) "Errors in compute_mixing_length subroutine"

   !        write(fstderr,*) "Intent(in)"

   !        write(fstderr,*) "thvm = ", thvm
   !        write(fstderr,*) "thlm = ", thlm
   !        write(fstderr,*) "rtm = ", rtm
   !        write(fstderr,*) "em = ", em
   !        write(fstderr,*) "exner = ", exner
   !        write(fstderr,*) "p_in_Pa = ", p_in_Pa
   !        write(fstderr,*) "thv_ds = ", thv_ds

   !        write(fstderr,*) "Intent(out)"

   !        write(fstderr,*) "Lscale = ", Lscale
   !        write(fstderr,*) "Lscale_up = ", Lscale_up
   !        write(fstderr,*) "Lscale_down = ", Lscale_down
     
   !      endif ! Fatal error

   !  end if

    return

end subroutine compute_mixing_length_standalone

elemental function sat_mixrat_liq_clubb(p_in_Pa, T_in_K) result(rsat)

   !-------------------------------------------------------------------------
   ! Description:
   !   Used to compute the saturation mixing ratio of liquid water.
   !   Taken from CLUBB.
   !
   ! References:
   !   Formula from Emanuel 1994, 4.4.14
   !-------------------------------------------------------------------------

    use shr_kind_mod,    only: r8 => shr_kind_r8

   implicit none

   ! Input Variables
   real(r8), intent(in) :: &
      p_in_Pa, & ! Pressure    [Pa]
      T_in_K     ! Temperature [K]

   real(r8), parameter :: &
      ep = 0.622_r8  ! Ratio of gas constants, Rd/Rv
   
   real(r8) :: rsat

   ! Local Variables
   real(r8) :: esatv

   ! Calculate the SVP for water vapor.
   esatv = sat_vapor_press_liq_flatau_clubb(T_in_K)

   ! If esatv exceeds the air pressure, then assume esatv~=0.5*pressure 
   ! and set rsat = ep = 0.622
   if (p_in_Pa-esatv < 1.0_r8) then
      rsat = ep
   else
      ! Formula for Saturation Mixing Ratio:
      !
      ! rs = (epsilon) * [ esat / ( p - esat ) ];
      ! where epsilon = R_d / R_v
      rsat = ep * esatv / (p_in_Pa - esatv)
   end if

   return
end function sat_mixrat_liq_clubb

elemental function sat_vapor_press_liq_flatau_clubb(T_in_K) result(esat)

    ! Description:
    !   Computes SVP (Saturation Vapor Pressure) for water vapor.
    !
    ! References:
    !   "Polynomial Fits to Saturation Vapor Pressure" Flatau, Walko,
    !   and Cotton. (1992) Journal of Applied Meteorology, Vol. 31,
    !   pp. 1507--1513
    !------------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8

   implicit none

   ! Constant parameters
   !
   ! Relative error norm expansion (-50 to 50 deg_C) from
   ! Table 3 of pp. 1510 of Flatau et al. 1992 (Water Vapor)
   ! (The 100 coefficient converts from mb to Pa)
   ! real, dimension(7), parameter :: a = & 
   ! 100.* (/ 6.11176750,      0.443986062,     0.143053301E-01, & 
   !          0.265027242E-03, 0.302246994E-05, 0.203886313E-07, & 
   !          0.638780966E-10 /)

   ! Relative error norm expansion (-85 to 70 deg_C) from
   ! Table 4 of pp. 1511 of Flatau et al.
   ! real(kind=r8), dimension(9), parameter :: a = &
   ! 100._r8 * &
   ! Commented out because the form has been redone, causing these number to no longer be needed,
   ! leaving them in for now for reference.
   !       (/ 6.11583699_r8,      0.444606896_r8,     0.143177157E-01_r8, &
   !       0.264224321E-03_r8, 0.299291081E-05_r8, 0.203154182E-07_r8, & 
   !       0.702620698E-10_r8, 0.379534310E-13_r8,-0.321582393E-15_r8 /)

   real(r8), parameter :: min_T_in_C = -85._r8, &  ! [deg_C]
                         T_freeze_K = 273.15_r8   ! [K]

   ! Input Variables
   real(r8), intent(in) :: T_in_K   ! Temperature [K]

   ! Output Variables
   real(r8) :: esat  ! Saturation vapor pressure over water [Pa]

   ! Local Variables
   real(r8) :: T_in_C, T_in_C_sqd
   ! integer :: i  ! Loop index

   ! ---- Begin Code ----

   ! Determine deg K - 273.15
   T_in_C = T_in_K - T_freeze_K

   ! Since this approximation is only good out to -85 degrees Celsius we
   ! truncate the result here (Flatau, et al. 1992)
   T_in_C = max(T_in_C, min_T_in_C)

   ! Polynomial approx. (Flatau, et al. 1992)

   ! This is the generalized formula but is not computationally efficient.
   ! Based on Wexler's expressions(2.1)-(2.4) (See Flatau et al. p 1508)
   ! e_{sat} = a_1 + a_2 (T - T_0) + ... + a_{n+1} (T - T_0)^n

   ! esat = a(1)
   ! do i = 2, size(a), 1
   !   esat = esat + a(i) * (T_in_C)**(i-1)
   ! end do

   ! The 8th order polynomial fit. When running deep 
   ! convective cases I noticed that absolute temperature often dips below
   ! -50 deg_C at higher altitudes, where the 6th order approximation is
   ! not accurate. -dschanen 20 Nov 2008
   !
   ! esat = a(1) + T_in_C*(a(2) + T_in_C*(a(3) + T_in_C*(a(4) + T_in_C &
   ! *(a(5) + T_in_C*(a(6) + T_in_C*(a(7) + T_in_C*(a(8) + T_in_C*(a(9))))))))

   ! Factoring the polynomial above and changing it into this form allows the cpu
   ! to complete the calculations out of order. This is because modern cpus can complete
   ! multiple instructions at once if they do not depend on each other. In the above case
   ! each instruction relies on the result of the last. In this version however, the terms
   ! in the parentheses could potentially be calculated in parallel by different execution
   ! units in the cpu, then only when those terms are being multiplied together do the 
   ! instructions need to be done one at a time. See clubb issue 834 for more info.
   !   - Gunther Huebler, Aug 2018
   
   T_in_C_sqd = T_in_C**2

   esat = -3.21582393e-14_r8 * (T_in_C - 646.5835252598777_r8) &
                        * (T_in_C + 90.72381630364440_r8) &
                        * (T_in_C_sqd + 111.0976961559954_r8 * T_in_C + 6459.629194243118_r8) &
                        * (T_in_C_sqd + 152.3131930092453_r8 * T_in_C + 6499.774954705265_r8) &
                        * (T_in_C_sqd + 174.4279584934021_r8 * T_in_C + 7721.679732114084_r8)

   return
end function sat_vapor_press_liq_flatau_clubb

end module lscale_mod
