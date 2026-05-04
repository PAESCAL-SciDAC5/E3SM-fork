module shoc_lscale_mod

    !---------------------------------------------------------------------------
    ! Purpose:
    !   SHOC-internal diagnostic Larson nonlocal moist length scale.
    !
    !   This is a SHOC-side parallel of the standalone routine
    !   `calculate_lscale` in [components/eam/src/physics/cam/lscale_mod.F90],
    !   which is used by the physics driver (tphysbc) to compute LSCALE,
    !   LSCALE_UP, LSCALE_DOWN once per macmic substep BEFORE CLUBB or SHOC
    !   runs.
    !
    !   This module exists because in turb_standalone (in-and-out) simulations
    !   with SHOC, the macmic loop runs only once per host time step, while
    !   SHOC subcycles internally for `nadv` substeps inside `shoc_main`.
    !   The driver-level lscale therefore captures only the initial state,
    !   not the per-substep evolution that the AMR/Part-1 paper requires.
    !
    !   Calling this module from inside SHOC's nadv loop produces lscale
    !   values consistent with each substep's evolved state. The placement
    !   inside `shoc_main` is right after `pblintd` and before `shoc_length`,
    !   structurally analogous to where CLUBB calls its own
    !   `compute_mixing_length` inside `advance_clubb_core`.
    !
    !   The algorithm is a direct port of the workhorse
    !   `compute_mixing_length_standalone` (Larson nonlocal moist length scale,
    !   Golaz et al. 2002 JAS Vol. 59, pp. 3540-3551, Section 3b "Eddy length
    !   formulation"), with two changes:
    !
    !   (a) Precision: `rtype` (configurable via SHOC's physics_utils) instead
    !       of hardcoded `r8`. SHOC may run in single precision under SCREAM
    !       CMAKE builds; the wrapper code stays consistent with the rest of
    !       SHOC.
    !
    !   (b) Inputs: takes SHOC's already-evolved thermodynamic state directly
    !       (thetal -> thlm, qw -> rtm, thv -> thvm and thv_ds, inv_exner is
    !       inverted to give CLUBB-convention exner). The standalone version
    !       in lscale_mod.F90 recomputes these from T, p, qv, qc using CLUBB
    !       formulas; here we skip that step because SHOC has the same
    !       quantities already computed and self-consistent with the rest of
    !       the SHOC closure.
    !
    !   Validation strategy: compare LSCALE_SHOC (from this module) against
    !   LSCALE (from lscale_mod, via tphysbc) at the first SHOC substep of
    !   the first macmic substep. They should match closely but not bit-for-
    !   bit (TKE source differs slightly: physpkg uses pbuf-stored interface-grid
    !   TKE from the previous host step; whereas this routine builds interface TKE
    !   on the fly from SHOC's prognostic midpoint-grid TKE via linear average).
    !
    ! Provenance:
    !   Original CLUBB algorithm: Larson, Golaz et al.
    !   Standalone E3SM port    : H. Xiao
    !   SHOC-side parallel port : M. Chinita, 2026
    !---------------------------------------------------------------------------

    use physics_utils, only: rtype

    implicit none

    private
    public :: shoc_compute_lscale

    contains

!=============================================================================
subroutine shoc_compute_lscale( &
                shcol, nlev, nlevi,                              &  ! in: dims
                thetal, qw, tke, thv,                            &  ! in: SHOC state on midpoints (top->bottom, no ghost)
                pres, inv_exner,                                 &  ! in: pres on midpoints, 1/Exner on midpoints
                zt_grid, zi_grid,                                &  ! in: heights AGL on midpoints/interfaces
                lscale_out, lscale_up_out, lscale_down_out)         ! out: midpoint, top->bottom, no ghost

    !-----------------------------------------------------------------------
    ! Wrapper: per-column translation between SHOC's grid convention
    ! (top-to-bottom, no ghost) and the CLUBB-internal convention required
    ! by `compute_mixing_length_shoc` (bottom-to-top, with one ghost level
    ! below the surface for thermodynamic-grid variables).
    !-----------------------------------------------------------------------

    implicit none

    !-----------------
    ! Arguments
    !-----------------
    integer,     intent(in) :: shcol, nlev, nlevi  ! nlevi = nlev + 1

    real(rtype), intent(in)  :: thetal   (shcol, nlev)    ! liq water pot. temp [K]
    real(rtype), intent(in)  :: qw       (shcol, nlev)    ! total water mixing ratio [kg/kg]
    real(rtype), intent(in)  :: tke      (shcol, nlev)    ! turbulent kinetic energy [m^2/s^2]
    real(rtype), intent(in)  :: thv      (shcol, nlev)    ! virtual potential temperature [K]
    real(rtype), intent(in)  :: pres     (shcol, nlev)    ! pressure [Pa]
    real(rtype), intent(in)  :: inv_exner(shcol, nlev)    ! 1/Exner = (p0/p)^(Rd/cp) [-]
    real(rtype), intent(in)  :: zt_grid  (shcol, nlev)    ! height AGL on midpoints [m]
    real(rtype), intent(in)  :: zi_grid  (shcol, nlevi)   ! height AGL on interfaces [m]

    real(rtype), intent(out) :: lscale_out     (shcol, nlev)
    real(rtype), intent(out) :: lscale_up_out  (shcol, nlev)
    real(rtype), intent(out) :: lscale_down_out(shcol, nlev)

    !-----------------------------------------------------------------------
    ! Local arrays. Their shape and direction of indexing follow CLUBB:
    !  - Variables defined on thermodynamic levels (corresponding to SHOC's
    !    layer midpoints) have an extra ghost level below the surface;
    !  - Level indices start from the surface and increase upward.
    !-----------------------------------------------------------------------
    real(rtype) :: zt       (nlevi)
    real(rtype) :: zm       (nlevi)
    real(rtype) :: dzm      (nlevi)
    real(rtype) :: invrs_dzm(nlevi)

    real(rtype) :: thvm   (nlevi)
    real(rtype) :: thlm   (nlevi)
    real(rtype) :: rtm    (nlevi)
    real(rtype) :: em     (nlevi)
    real(rtype) :: p_in_Pa(nlevi)
    real(rtype) :: exner  (nlevi)        ! NB: this stores (p/p0)^(Rd/cp) =
                                          ! 1/inv_exner; this matches the
                                          ! convention that
                                          ! compute_mixing_length_shoc expects
                                          ! (T = thetal * exner)
    real(rtype) :: thv_ds (nlevi)

    real(rtype) :: lscale_tmp     (nlevi)
    real(rtype) :: lscale_up_tmp  (nlevi)
    real(rtype) :: lscale_down_tmp(nlevi)

    integer :: ii, kk, kflip, k_shoc

    !-----------------------------------------------------------------------
    ! Process each column independently.
    !-----------------------------------------------------------------------
    do ii = 1, shcol

       !====================================================
       ! 1. Prepare input variables for this column,
       !    flipping vertical indexing to CLUBB's convention.
       !====================================================
       do kk = 1, nlev

          !--------------------------------------------------
          ! 1.1 Thermodynamic-level variables
          !     CLUBB index kk+1 (slot 1 reserved for ghost),
          !     SHOC index kflip = nlev - kk + 1.
          !--------------------------------------------------
          kflip = nlev - kk + 1

          thlm   (kk+1) = thetal   (ii, kflip)
          rtm    (kk+1) = qw       (ii, kflip)
          thvm   (kk+1) = thv      (ii, kflip)
          ! thv_ds is set to the moist theta_v, following the convention used
          ! in clubb_intr.F90 (clubb_tend_cam: thv_ds_zt = moist thv) and
          ! mirrored in calculate_lscale (lscale_mod.F90:142). Inside
          ! compute_mixing_length, this is treated as a "reference theta"
          ! used to recompute thv from (thl, rt, qc) for the parcel.
          thv_ds (kk+1) = thv      (ii, kflip)
          p_in_Pa(kk+1) = pres     (ii, kflip)
          ! SHOC's inv_exner = (p0/p)^(R/cp). The algorithm needs the
          ! "standard" Exner = (p/p0)^(R/cp) = 1/inv_exner.
          exner  (kk+1) = 1.0_rtype / inv_exner(ii, kflip)
          ! zt_grid is already height above ground (subtraction of surface
          ! interface height done in shoc_intr.F90 before shoc_main is called).
          zt     (kk+1) = zt_grid  (ii, kflip)

       end do

       !-----------------------------------------------------------
       ! 1.2 Fill ghost level (slot 1, below ground) for zt-grid
       !     variables. Mirror zt; zero-gradient for the rest.
       !-----------------------------------------------------------
       zt     (1) = -1.0_rtype * zt(2)
       exner  (1) = exner  (2)
       p_in_Pa(1) = p_in_Pa(2)
       rtm    (1) = rtm    (2)
       thlm   (1) = thlm   (2)
       thv_ds (1) = thv_ds (2)
       thvm   (1) = thvm   (2)

       !-----------------------------------------------------------
       ! 1.3 Momentum-level variables (interfaces). NO ghost level
       !     on this grid -- all `nlevi` slots are real interfaces,
       !     reordered surface-up.
       !-----------------------------------------------------------
       ! Heights: zi_grid is already AGL.
       do kk = 1, nlev
          zm(kk) = zi_grid(ii, nlevi - kk + 1)
       end do
       zm(nlevi) = zi_grid(ii, 1)   ! TOM interface

       !-----------------------------------------------------------
       ! 1.4 TKE on momentum (interface) levels.
       !
       !     Inside shoc_main, only midpoint TKE is available (the
       !     interpolation to interface grid happens in shoc_intr,
       !     via linear_interp, AFTER shoc_main returns and writes
       !     to pbuf 'tke').
       !
       !     We replicate that same linear-in-z interpolation here
       !     so LSCALE_SHOC and the driver-level LSCALE see
       !     consistently-interpolated interface TKE. The remaining
       !     LSCALE_SHOC vs LSCALE difference is then the
       !     time-state difference only (LSCALE uses pbuf 'tke'
       !     populated at the previous host step; LSCALE_SHOC uses
       !     within-substep midpoint TKE interpolated on the fly).
       !
       !     For interior interfaces (k_shoc = 2..nlev, between two
       !     midpoints), the formula matches linear_interp's
       !     interior branch [shoc.F90:4795-4798]. The surface and
       !     TOM interfaces use linear extrapolation from the two
       !     nearest midpoints, also matching linear_interp's
       !     boundary branches [4791-4793, 4801-4803].
       !-----------------------------------------------------------

       ! Interior interfaces. CLUBB kk=2..nlev maps to SHOC
       ! interface index k_shoc = nlevi-kk+1 in 2..nlev.
       do kk = 2, nlev
          k_shoc = nlevi - kk + 1
          ! Surrounding SHOC midpoints: k_shoc-1 (above), k_shoc (below).
          em(kk) = tke(ii, k_shoc-1) &
                   + ( tke(ii, k_shoc) - tke(ii, k_shoc-1) ) &
                   * ( zi_grid(ii, k_shoc) - zt_grid(ii, k_shoc-1) ) &
                   / ( zt_grid(ii, k_shoc) - zt_grid(ii, k_shoc-1) )
       end do

       ! Surface interface (CLUBB kk=1, SHOC k_shoc=nlevi):
       ! linear extrapolation using the two lowest midpoints.
       em(1) = tke(ii, nlev-1) &
               + ( tke(ii, nlev) - tke(ii, nlev-1) ) &
               * ( zi_grid(ii, nlevi) - zt_grid(ii, nlev-1) ) &
               / ( zt_grid(ii, nlev)  - zt_grid(ii, nlev-1) )

       ! Top-of-model interface (CLUBB kk=nlevi, SHOC k_shoc=1):
       ! linear extrapolation using the two highest midpoints.
       em(nlevi) = tke(ii, 1) &
                   + ( tke(ii, 2) - tke(ii, 1) ) &
                   * ( zi_grid(ii, 1) - zt_grid(ii, 1) ) &
                   / ( zt_grid(ii, 2) - zt_grid(ii, 1) )

       !-----------------------------------------------------------
       ! 1.5 Grid spacing and inverse spacing on momentum grid.
       !-----------------------------------------------------------
       do kk = 1, nlev
          dzm      (kk) = zt(kk+1) - zt(kk)
          invrs_dzm(kk) = 1.0_rtype / dzm(kk)
       end do
       dzm      (nlevi) = dzm      (nlev)   ! never used inside the algorithm
       invrs_dzm(nlevi) = invrs_dzm(nlev)

       !====================================================
       ! 2. Run the algorithm on this column.
       !====================================================
       call compute_mixing_length_shoc( &
            nlevi, zt, zm, dzm, invrs_dzm, &
            thvm, thlm, rtm, em, p_in_Pa, exner, thv_ds, &
            lscale_tmp, lscale_up_tmp, lscale_down_tmp )

       !====================================================
       ! 3. Flip outputs back to SHOC convention (top-to-
       !    bottom, no ghost). The ghost level (CLUBB index 1)
       !    is dropped on the way out.
       !====================================================
       do kk = 1, nlev
          lscale_out     (ii, kk) = lscale_tmp     (nlevi - kk + 1)
          lscale_up_out  (ii, kk) = lscale_up_tmp  (nlevi - kk + 1)
          lscale_down_out(ii, kk) = lscale_down_tmp(nlevi - kk + 1)
       end do

    end do  ! ii = 1, shcol

end subroutine shoc_compute_lscale

!=============================================================================
subroutine compute_mixing_length_shoc( &
    nzmax, zt, zm, dzm, invrs_dzm, &
    thvm, thlm, rtm, em, p_in_Pa, exner, thv_ds, &
    Lscale, Lscale_up, Lscale_down )

    !---------------------------------------------------------------------
    ! Larson's 5th moist, nonlocal length scale.
    ! Direct port of `compute_mixing_length_standalone` from
    ! [components/eam/src/physics/cam/lscale_mod.F90:209-978] by H. Xiao,
    ! with r8 -> rtype substitution. The algorithm is unchanged.
    !
    ! References:
    !   Section 3b ("Eddy length formulation") of
    !   "A PDF-Based Model for Boundary Layer Clouds. Part I:
    !    Method and Model Description", Golaz, et al. (2002),
    !   JAS, Vol. 59, pp. 3540--3551.
    !
    ! See lscale_mod.F90 for full algorithmic documentation; the comments
    ! in this routine are kept terse to avoid duplication.
    !---------------------------------------------------------------------

    use shr_const_mod, only: shr_const_rdair, shr_const_cpdair, shr_const_latvap, &
                             shr_const_mwwv, shr_const_mwdair, shr_const_g

    implicit none

    ! Constants. Use the same host-model constants the lscale_mod standalone
    ! uses, so LSCALE_SHOC and LSCALE share thermodynamic constants and the
    ! validation comparison isn't biased by constant choice.
    real(rtype), parameter :: &
        Cp   = shr_const_cpdair, &
        Lv   = shr_const_latvap, &
        Rd   = shr_const_rdair,  &
        ep   = shr_const_mwwv / shr_const_mwdair, &  ! ~ 0.622
        ep1  = (1.0_rtype - ep) / ep, &              ! ~ 0.61
        ep2  = 1.0_rtype / ep, &                     ! ~ 1.61
        grav = shr_const_g, &
        zero_threshold = 0.0_rtype, &
        eps_div = 1.0e-10_rtype  ! small number to guard divides

    intrinsic :: min, max, sqrt

    ! Almost-Constant Parameters (matching lscale_mod.F90 exactly so
    ! validation against LSCALE is meaningful).
    real(rtype), parameter ::  &
      zlmin               = 0.1_rtype,       &           ! [m]
      Lscale_sfclyr_depth = 500.0_rtype,     &           ! [m]
      em_min              = 1.5_rtype * (2.0e-2_rtype)**2, &  ! [m^2/s^2]
      lmin                = 40.0_rtype * 0.1_rtype       ! [m]

    ! Tunable parameters (hardcoded to match lscale_mod.F90 / clubb_mu = 5e-4).
    ! TODO: expose `mu` via SHOC namelist if SHOC tuning warrants.
    ! TODO: replace `Lscale_max` hardcode with 0.25 * min(host_dx, host_dy)
    !       once we plumb the host-grid spacing into this module.
    real(rtype), parameter ::  &
      mu         = 5.0e-4_rtype, &
      Lscale_max = 0.25_rtype * 25000.0_rtype

    ! Input/Output
    integer,     intent(in)  :: nzmax
    real(rtype), dimension(nzmax), intent(in) :: &
      thvm, thlm, rtm, em, exner, p_in_Pa, thv_ds, &
      zt, zm, dzm, invrs_dzm
    real(rtype), dimension(nzmax), intent(out) :: &
      Lscale, Lscale_up, Lscale_down

    ! Local Variables
    integer :: i, j

    real(rtype) :: tke, CAPE_incr
    real(rtype) :: dCAPE_dz_j, dCAPE_dz_j_minus_1, dCAPE_dz_j_plus_1

    real(rtype), dimension(nzmax) :: &
        exp_mu_dzm, invrs_dzm_on_mu, grav_on_thvm, &
        thl_par_j_precalc, rt_par_j_precalc, &
        tke_i, &
        thl_par_1, rt_par_1, tl_par_1, &
        rsatl_par_1, s_par_1, rc_par_1, thv_par_1, &
        dCAPE_dz_1, CAPE_incr_1, &
        Lv_coef, entrain_coef

    real(rtype) :: lminh
    real(rtype) :: thl_par_j, rt_par_j, rc_par_j, thv_par_j
    real(rtype) :: tl_par_j, rsatl_par_j, s_par_j
    real(rtype) :: Lscale_up_max_alt, Lscale_down_min_alt
    real(rtype) :: Lv2_coef, tl_par_j_sqd, invrs_dCAPE_diff, &
                    invrs_Lscale_sfclyr_depth

    ! ---- Begin Code ----

    ! Initialize arrays and precalculate values for computational efficiency.
    do i = 1, nzmax
        Lscale_up(i)        = zlmin
        Lscale_down(i)      = zlmin
        exp_mu_dzm(i)       = exp( -mu * dzm(i) )
        invrs_dzm_on_mu(i)  = ( invrs_dzm(i) ) / mu
        grav_on_thvm(i)     = grav / thvm(i)
        Lv_coef(i)          = Lv / ( exner(i) * Cp ) - ep2 * thv_ds(i)
        entrain_coef(i)     = ( 1.0_rtype - exp_mu_dzm(i) ) * invrs_dzm_on_mu(i)
    end do

    ! Avoid uninitialized memory at the ghost level (these are dropped on output).
    Lscale_up  (1) = 0.0_rtype
    Lscale_down(1) = 0.0_rtype

    Lv2_coef                  = ep * Lv**2 / ( Rd * Cp )
    invrs_Lscale_sfclyr_depth = 1.0_rtype / Lscale_sfclyr_depth

    ! Calculate initial TKE for zt levels by linear interpolation, following CLUBB.
    do i = 2, nzmax
      tke_i(i) = em(i)   * ( zt(i) - zm(i-1) ) + &
                 em(i-1) * ( zm(i) - zt(i)   )
      tke_i(i) = max( tke_i(i) / (zm(i) - zm(i-1)), em_min )
    end do
    tke_i(1) = (zt(1) - zm(1)) * (em(2) - em(1)) / (zm(2) - zm(1)) + em(1)
    tke_i(1) = max( tke_i(1), em_min )

    ! ============== Upwards Length Scale Calculation ==============

    ! Precalculate values for upward Lscale.
    do j = 2, nzmax-1
        thl_par_j_precalc(j) = thlm(j) - thlm(j-1) * exp_mu_dzm(j-1)  &
                               - ( thlm(j) - thlm(j-1) ) * entrain_coef(j-1)
        rt_par_j_precalc(j)  = rtm(j) - rtm(j-1) * exp_mu_dzm(j-1)  &
                               - ( rtm(j) - rtm(j-1) ) * entrain_coef(j-1)
    end do

    ! Initial parcel properties at each grid level (vectorized first level).
    do j = 3, nzmax
        thl_par_1(j) = thlm(j) - ( thlm(j) - thlm(j-1) ) * entrain_coef(j-1)
        tl_par_1 (j) = thl_par_1(j) * exner(j)
        rt_par_1 (j) = rtm(j)  - ( rtm(j)  - rtm(j-1)  ) * entrain_coef(j-1)
    end do

    rsatl_par_1(3:) = sat_mixrat_liq_shoc( p_in_Pa(3:), tl_par_1(3:) )

    do j = 3, nzmax
        tl_par_j_sqd = tl_par_1(j)**2

        ! s = (rt - rsatl) / (1 + beta*rsatl); simplified by *tl^2.
        s_par_1 (j) = ( rt_par_1(j) - rsatl_par_1(j) ) * tl_par_j_sqd &
                      / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(j) )

        rc_par_1(j) = max( s_par_1(j), zero_threshold )
        thv_par_1(j) = thl_par_1(j) + ep1 * thv_ds(j) * rt_par_1(j) + Lv_coef(j) * rc_par_1(j)

        dCAPE_dz_1 (j) = grav_on_thvm(j) * ( thv_par_1(j) - thvm(j) )
        ! Trapezoid; dCAPE at z_0 = 0 for this initial calculation.
        CAPE_incr_1(j) = 0.5_rtype * dCAPE_dz_1(j) * dzm(j-1)
    end do

    Lscale_up_max_alt = 0.0_rtype
    do i = 2, nzmax-2

        if ( tke_i(i) + CAPE_incr_1(i+1) > 0.0_rtype ) then

            tke = tke_i(i) + CAPE_incr_1(i+1)
            j   = i + 2

            thl_par_j          = thl_par_1(i+1)
            rt_par_j           = rt_par_1(i+1)
            dCAPE_dz_j_minus_1 = dCAPE_dz_1(i+1)

            do while ( j < nzmax )

                thl_par_j = thl_par_j_precalc(j) + thl_par_j * exp_mu_dzm(j-1)
                rt_par_j  = rt_par_j_precalc (j) + rt_par_j  * exp_mu_dzm(j-1)

                tl_par_j     = thl_par_j * exner(j)
                rsatl_par_j  = sat_mixrat_liq_shoc( p_in_Pa(j), tl_par_j )
                tl_par_j_sqd = tl_par_j**2

                s_par_j  = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
                           / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )
                rc_par_j = max( s_par_j, zero_threshold )

                thv_par_j  = thl_par_j + ep1 * thv_ds(j) * rt_par_j  &
                             + Lv_coef(j) * rc_par_j

                dCAPE_dz_j = grav_on_thvm(j) * ( thv_par_j - thvm(j) )
                CAPE_incr  = 0.5_rtype * ( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * dzm(j-1)

                if ( tke + CAPE_incr <= 0.0_rtype ) exit

                dCAPE_dz_j_minus_1 = dCAPE_dz_j
                tke                = tke + CAPE_incr
                j                  = j + 1
            end do

            Lscale_up(i) = Lscale_up(i) + zt(j-1) - zt(i)

            if ( j < nzmax ) then
                if ( abs( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) * 2 <= &
                     abs( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * eps_div ) then
                    ! dCAPE/dz constant -> linear-in-z exhaustion.
                    Lscale_up(i) = Lscale_up(i) + ( - tke / dCAPE_dz_j )
                else
                    ! Generic case: quadratic in z.
                    invrs_dCAPE_diff = 1.0_rtype / ( dCAPE_dz_j - dCAPE_dz_j_minus_1 )
                    Lscale_up(i) = Lscale_up(i) &
                                   - dCAPE_dz_j_minus_1 * invrs_dCAPE_diff * dzm(j-1) &
                                   - sqrt( dCAPE_dz_j_minus_1**2 &
                                            - 2.0_rtype * tke * invrs_dzm(j-1) &
                                              * ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) ) &
                                     * invrs_dCAPE_diff * dzm(j-1)
                end if
            end if

        else  ! TKE exhausted within first grid level
            Lscale_up(i) = Lscale_up(i) - sqrt( - 2.0_rtype * tke_i(i) &
                                                 * dzm(i) * dCAPE_dz_1(i+1) ) &
                                          / dCAPE_dz_1(i+1)
        end if

        ! Monotonicity (smoothing) constraint for nonlocality.
        if ( zt(i) + Lscale_up(i) < Lscale_up_max_alt ) then
            Lscale_up(i) = Lscale_up_max_alt - zt(i)
        else
            Lscale_up_max_alt = Lscale_up(i) + zt(i)
        end if

    end do

    ! ============== Downwards Length Scale Calculation ==============

    do j = 2, nzmax-1
        thl_par_j_precalc(j) = thlm(j) - thlm(j+1) * exp_mu_dzm(j)  &
                               - ( thlm(j) - thlm(j+1) ) * entrain_coef(j)
        rt_par_j_precalc(j)  = rtm(j) - rtm(j+1) * exp_mu_dzm(j)  &
                               - ( rtm(j) - rtm(j+1) ) * entrain_coef(j)
    end do

    do j = 2, nzmax-1
        thl_par_1(j) = thlm(j) - ( thlm(j) - thlm(j+1) )  * entrain_coef(j)
        tl_par_1 (j) = thl_par_1(j) * exner(j)
        rt_par_1 (j) = rtm(j)  - ( rtm(j)  - rtm(j+1)  ) * entrain_coef(j)
    end do

    rsatl_par_1(2:) = sat_mixrat_liq_shoc( p_in_Pa(2:), tl_par_1(2:) )

    do j = 2, nzmax-1
        tl_par_j_sqd = tl_par_1(j)**2
        s_par_1 (j) = ( rt_par_1(j) - rsatl_par_1(j) ) * tl_par_j_sqd &
                      / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(j) )
        rc_par_1(j) = max( s_par_1(j), zero_threshold )
        thv_par_1(j) = thl_par_1(j) + ep1 * thv_ds(j) * rt_par_1(j) + Lv_coef(j) * rc_par_1(j)

        dCAPE_dz_1 (j) = grav_on_thvm(j) * ( thv_par_1(j) - thvm(j) )
        CAPE_incr_1(j) = 0.5_rtype * dCAPE_dz_1(j) * dzm(j)
    end do

    Lscale_down_min_alt = zt(nzmax)
    do i = nzmax, 3, -1

        if ( tke_i(i) - CAPE_incr_1(i-1) > 0.0_rtype ) then

            tke = tke_i(i) - CAPE_incr_1(i-1)
            j   = i - 2

            thl_par_j         = thl_par_1(i-1)
            rt_par_j          = rt_par_1(i-1)
            dCAPE_dz_j_plus_1 = dCAPE_dz_1(i-1)

            do while ( j >= 2 )

                thl_par_j = thl_par_j_precalc(j) + thl_par_j * exp_mu_dzm(j)
                rt_par_j  = rt_par_j_precalc (j) + rt_par_j  * exp_mu_dzm(j)

                tl_par_j     = thl_par_j * exner(j)
                rsatl_par_j  = sat_mixrat_liq_shoc( p_in_Pa(j), tl_par_j )
                tl_par_j_sqd = tl_par_j**2

                s_par_j  = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
                           / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )
                rc_par_j = max( s_par_j, zero_threshold )

                thv_par_j = thl_par_j + ep1 * thv_ds(j) * rt_par_j + Lv_coef(j) * rc_par_j

                dCAPE_dz_j = grav_on_thvm(j) * ( thv_par_j - thvm(j) )
                CAPE_incr  = 0.5_rtype * ( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * dzm(j)

                if ( tke - CAPE_incr <= 0.0_rtype ) exit

                dCAPE_dz_j_plus_1 = dCAPE_dz_j
                tke               = tke - CAPE_incr
                j                 = j - 1
            end do

            Lscale_down(i) = Lscale_down(i) + zt(i) - zt(j+1)

            if ( j >= 2 ) then
                if ( abs( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) * 2 <= &
                     abs( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * eps_div ) then
                    Lscale_down(i) = Lscale_down(i) + ( tke / dCAPE_dz_j )
                else
                    invrs_dCAPE_diff = 1.0_rtype / ( dCAPE_dz_j - dCAPE_dz_j_plus_1 )
                    Lscale_down(i) = Lscale_down(i) &
                                     - dCAPE_dz_j_plus_1 * invrs_dCAPE_diff * dzm(j) &
                                     + sqrt( dCAPE_dz_j_plus_1**2 &
                                             + 2.0_rtype * tke * invrs_dzm(j)  &
                                               * ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) ) &
                                       * invrs_dCAPE_diff * dzm(j)
                end if
            end if

        else
            Lscale_down(i) = Lscale_down(i) + sqrt( 2.0_rtype * tke_i(i) &
                                                    * dzm(i-1) * dCAPE_dz_1(i-1) ) &
                                              / dCAPE_dz_1(i-1)
        end if

        if ( zt(i) - Lscale_down(i) > Lscale_down_min_alt ) then
            Lscale_down(i) = zt(i) - Lscale_down_min_alt
        else
            Lscale_down_min_alt = zt(i) - Lscale_down(i)
        end if

    end do

    ! ============== Final Lscale Calculation ==============

    do i = 2, nzmax, 1
        ! Linear taper of the surface-layer minimum from `lmin` at the
        ! surface to zero at Lscale_sfclyr_depth above the ground.
        lminh = max( zero_threshold, Lscale_sfclyr_depth - ( zt(i) - zm(1) ) ) &
                    * lmin * invrs_Lscale_sfclyr_depth
        Lscale_up  (i) = max( lminh, Lscale_up  (i) )
        Lscale_down(i) = max( lminh, Lscale_down(i) )
        Lscale     (i) = sqrt( Lscale_up(i) * Lscale_down(i) )
    end do

    ! Boundary values (top and ghost). The ghost (i=1) is dropped on output.
    Lscale(1)     = Lscale(2)
    Lscale(nzmax) = Lscale(nzmax-1)

    ! Cap at maximum (intended to let host model take over deep convection).
    Lscale = min( Lscale, Lscale_max )

    return
end subroutine compute_mixing_length_shoc

!=============================================================================
elemental function sat_mixrat_liq_shoc(p_in_Pa, T_in_K) result(rsat)
   !---------------------------------------------------------------------
   ! Saturation mixing ratio of liquid water (Emanuel 1994, eqn. 4.4.14).
   ! Direct port of `sat_mixrat_liq_clubb` from lscale_mod.F90:980-1024.
   !---------------------------------------------------------------------

   implicit none

   real(rtype), intent(in) :: p_in_Pa, T_in_K

   real(rtype), parameter :: ep = 0.622_rtype  ! Rd/Rv

   real(rtype) :: rsat
   real(rtype) :: esatv

   esatv = sat_vapor_press_liq_flatau_shoc(T_in_K)

   ! Guard against pressure approaching saturation vapor pressure.
   if (p_in_Pa - esatv < 1.0_rtype) then
      rsat = ep
   else
      rsat = ep * esatv / (p_in_Pa - esatv)
   end if

   return
end function sat_mixrat_liq_shoc

!=============================================================================
elemental function sat_vapor_press_liq_flatau_shoc(T_in_K) result(esat)
   !---------------------------------------------------------------------
   ! Saturation vapor pressure over water via the Flatau et al. (1992)
   ! 8th-order polynomial fit (Table 4 of Flatau, Walko & Cotton 1992,
   ! J. Appl. Met. 31, 1507-1513), in the factored form computed by
   ! G. Huebler (CLUBB issue 834) for cpu pipelining.
   ! Direct port of `sat_vapor_press_liq_flatau_clubb` from
   ! lscale_mod.F90:1026-1120.
   !---------------------------------------------------------------------

   implicit none

   real(rtype), parameter :: min_T_in_C = -85.0_rtype, &
                              T_freeze_K = 273.15_rtype

   real(rtype), intent(in) :: T_in_K

   real(rtype) :: esat
   real(rtype) :: T_in_C, T_in_C_sqd

   T_in_C = T_in_K - T_freeze_K
   ! The polynomial is only valid out to -85 deg_C.
   T_in_C = max(T_in_C, min_T_in_C)

   T_in_C_sqd = T_in_C**2

   esat = -3.21582393e-14_rtype * (T_in_C - 646.5835252598777_rtype) &
                               * (T_in_C +  90.72381630364440_rtype) &
                               * (T_in_C_sqd + 111.0976961559954_rtype * T_in_C + 6459.629194243118_rtype) &
                               * (T_in_C_sqd + 152.3131930092453_rtype * T_in_C + 6499.774954705265_rtype) &
                               * (T_in_C_sqd + 174.4279584934021_rtype * T_in_C + 7721.679732114084_rtype)

   return
end function sat_vapor_press_liq_flatau_shoc

end module shoc_lscale_mod
