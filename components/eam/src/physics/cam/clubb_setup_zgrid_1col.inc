
  subroutine setup_clubb_zgrid_for_1column( pver, pverp, host_zi, host_zm, &
                                            zt_bot, sfc_elevation,         &
                                            zt_out_1col, zi_out_1col )

      use clubb_api_module, only: read_parameters_api, nparams
      use clubb_api_module, only: setup_grid_heights_api, setup_parameters_api

      integer,        intent(in)  :: pver
      integer,        intent(in)  :: pverp
      real(r8),       intent(in)  :: host_zi(:)
      real(r8),       intent(in)  :: host_zm(:)
      real(core_rknd),intent(out) :: zt_bot
      real(core_rknd),intent(out) :: sfc_elevation
      real(r8),       intent(out) :: zt_out_1col(pverp)  ! height of CLUBB's thermo   levels, for history output [m]
      real(r8),       intent(out) :: zi_out_1col(pverp)  ! height of CLUBB's momentum levels, for history output [m]

      real(core_rknd) :: zi_g(pverp)
      real(core_rknd) :: zt_g(pverp)

      integer :: k
      integer :: err_code

      logical,parameter :: l_implemented = .true.
      real(core_rknd), dimension(nparams)  :: clubb_params ! Adjustable CLUBB parameters (C1, C2 ...)

      !--------
      !  Define the CLUBB momentum grid (in height, units of m)
      do k=1,pverp
         zi_g(k) = real( host_zi(pverp-k+1)-host_zi(pver+1), kind=core_rknd )
      enddo

      !  Define the CLUBB thermodynamic grid (in units of m)
      do k=1,pver
         zt_g(k+1) = real( host_zm(pver-k+1)-host_zi(pver+1), kind=core_rknd )
      enddo
      zt_g(1) = -1._core_rknd*zt_g(2)  ! ghost point is below surface

      !  Read in parameters for CLUBB. Pack the default and updated (via nml) tunable parameters into clubb_params
      call read_parameters_api( -99, "", clubb_params )

      !  Setup CLUBB's z-grid because heights are given by the host model and hence can change

      call setup_grid_heights_api(l_implemented, grid_type, zi_g(2), zi_g(1), zi_g, zt_g)
      call setup_parameters_api(zi_g(2), clubb_params, pverp, grid_type, zi_g, zt_g, err_code)

      !----------------------------------
      !----------------------------------
      zt_bot = zt_g(2)  ! lowst thermo level above ground. Used for specifying ustar in single column runs
      sfc_elevation = real( host_zi(pver+1), kind=core_rknd )  ! the elevation of the surface. Input to clubb_core

      !----------------------------------
      do k=1,pverp
         zt_out_1col(k) = real( zt_g(pverp-k+1), kind=r8 )
         zi_out_1col(k) = real( zi_g(pverp-k+1), kind=r8 )
      enddo

      zi_out_1col(1) = 0._r8

  end subroutine setup_clubb_zgrid_for_1column
