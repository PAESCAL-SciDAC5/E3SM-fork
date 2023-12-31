  subroutine enchance_thlp2_for_cloud_top_cooling( rcm_zt, qrl_zt, thlprcp, thlp2, nz, dtime )

    use clubb_api_module, only: zt2zm_api, calculate_thlp2_rad_api

    real(core_rknd),intent(in)    :: rcm_zt (nz)
    real(core_rknd),intent(in)    :: qrl_zt (nz)
    real(core_rknd),intent(in)    :: thlprcp(nz)
    real(core_rknd),intent(inout) :: thlp2  (nz)

    integer,        intent(in) :: nz
    real(core_rknd),intent(in) :: dtime

    real(core_rknd) :: rcm_zm (nz)
    real(core_rknd) :: qrl_zm (nz)
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
