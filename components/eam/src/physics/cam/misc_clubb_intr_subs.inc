  !--------------------------------------------------------------------------------------------
  ! Code separated from clubb_tend_cam by Hui Wan, 2023-12
  !--------------------------------------------------------------------------------------------
  subroutine update_mnts_for_rain_evap( rtp2, thlp2, wprtp, wpthlp,                 &! inout
                                        rtm, rcm, cloud_frac, thlm, wm_zt, exner,   &! in
                                        pre_in, pdf_params, landfrac, rnevap_effic, &! in
                                        nz, dtime, clubb_do_deep                    )! in

    use clubb_api_module, only: update_xp2_mc_api

    real(core_rknd),intent(inout)  :: rtp2  (nz)
    real(core_rknd),intent(inout)  :: thlp2 (nz)
    real(core_rknd),intent(inout)  :: wprtp (nz)
    real(core_rknd),intent(inout)  :: wpthlp(nz)

    real(core_rknd),intent(in)     :: rtm       (nz)
    real(core_rknd),intent(in)     :: rcm       (nz)
    real(core_rknd),intent(in)     :: cloud_frac(nz)
    real(core_rknd),intent(in)     :: thlm      (nz)
    real(core_rknd),intent(in)     :: wm_zt     (nz)
    real(core_rknd),intent(in)     :: exner     (nz)
    real(core_rknd),intent(in)     :: pre_in    (nz)

    type(pdf_parameter),intent(in) :: pdf_params

    real(core_rknd),intent(in)     :: landfrac
    real(core_rknd),intent(in)     :: rnevap_effic
    integer,        intent(in)     :: nz
    real(core_rknd),intent(in)     :: dtime
    logical,        intent(in)     :: clubb_do_deep

    real(core_rknd) ::         rvm_in(nz)
    real(core_rknd) ::    rtp2_mc_out(nz)
    real(core_rknd) ::   thlp2_mc_out(nz)
    real(core_rknd) ::   wprtp_mc_out(nz)
    real(core_rknd) ::  wpthlp_mc_out(nz)
    real(core_rknd) :: rtpthlp_mc_out(nz)

    real(core_rknd) :: dum_core_rknd

    !---------------------
    rvm_in = rtm - rcm
    call update_xp2_mc_api( nz, dtime, cloud_frac, rcm, rvm_in, thlm, wm_zt, exner, pre_in, pdf_params, &
                            rtp2_mc_out, thlp2_mc_out, wprtp_mc_out, wpthlp_mc_out, rtpthlp_mc_out)

    if (clubb_do_deep) then
       dum_core_rknd = 1._core_rknd
    else
       dum_core_rknd = 1._core_rknd - landfrac
    end if

    ! update turbulent moments based on rain evaporation

    rtp2  = rtp2  + rnevap_effic * dum_core_rknd *  rtp2_mc_out * dtime
    thlp2 = thlp2 + rnevap_effic * dum_core_rknd * thlp2_mc_out * dtime

    if (.not. clubb_do_deep) then
       wprtp  = wprtp  + rnevap_effic * dum_core_rknd *  wprtp_mc_out * dtime
       wpthlp = wpthlp + rnevap_effic * dum_core_rknd * wpthlp_mc_out * dtime
    endif

  end subroutine update_mnts_for_rain_evap

  !--------------------------------------------------------------------------------------------
  ! Code separated from clubb_tend_cam by Hui Wan, 2023-12
  !--------------------------------------------------------------------------------------------
  subroutine enchance_thlp2_for_cloud_top_cooling( thlp2, rcm_zt, qrl_zt, thlprcp, nz, dtime )

    use clubb_api_module, only: zt2zm_api, calculate_thlp2_rad_api

    real(core_rknd),intent(inout) :: thlp2  (nz)
    real(core_rknd),intent(in)    :: rcm_zt (nz)
    real(core_rknd),intent(in)    :: qrl_zt (nz)
    real(core_rknd),intent(in)    :: thlprcp(nz)
    integer,        intent(in)    :: nz
    real(core_rknd),intent(in)    :: dtime

    real(core_rknd) :: rcm_zm       (nz)
    real(core_rknd) :: qrl_zm       (nz)
    real(core_rknd) :: thlp2_rad_out(nz)

    !----------------------------------------------------------
    ! Interpolate arrays from thermo levels to momentum levels
    !----------------------------------------------------------
    rcm_zm = zt2zm_api( rcm_zt )
    qrl_zm = zt2zm_api( qrl_zt )

    !-------------------------------------------
    ! Calculate thlp2 tendency due to radiation
    !-------------------------------------------
    thlp2_rad_out(:) = 0._core_rknd
    call calculate_thlp2_rad_api( nz, rcm_zm, thlprcp, qrl_zm, thlp2_rad_out)

    !--------------
    ! Update thlp2
    !--------------
    thlp2(:) = thlp2(:) + thlp2_rad_out(:) * dtime
    thlp2(:) = max( thl_tol**2, thlp2(:) )

  end subroutine enchance_thlp2_for_cloud_top_cooling
