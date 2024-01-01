

  subroutine advance_clubb_core_api_eam( nz, dtime, sfc_elevation, fcoriolis, &! in
                                         core_forcing, core_sfc,     &! in
                                         clubb_misc,   core_auxil,   &! in
                                         core_prog,    core_diag,    &! inout or out
                                         pdf_params,   pdf_params_zm,&! inout or out
                                         !
                                         upwp_sfc_pert, vpwp_sfc_pert,      &! intent(in)
                                         um_pert_col, vm_pert_col,          &! 
                                         upwp_pert_col, vpwp_pert_col       &! intent(inout)
                                         )

     use clubb_intr_core_types, only: core_auxil_t, core_prog_t, core_diag_t, core_forcing_t, core_sfc_t, clubb_misc_t
     use clubb_api_module,      only: pdf_parameter
     use clubb_api_module,      only: advance_clubb_core_api

     use clubb_api_module,      only: clubb_fatal_error, fstderr
     use shr_log_mod ,          only: errMsg => shr_log_errMsg
     use cam_abortutils,        only: endrun
     use phys_grid,             only: get_gcol_p
     use time_manager,          only: get_nstep

     integer, intent(in) :: nz

     real(core_rknd),intent(in) :: dtime 
     real(core_rknd),intent(in) :: sfc_elevation
     real(core_rknd),intent(in) :: fcoriolis

     type(core_forcing_t) :: core_forcing
     type(core_sfc_t)     :: core_sfc
     type(core_auxil_t)   :: core_auxil
     type(core_prog_t)    :: core_prog
     type(core_diag_t)    :: core_diag
     type(clubb_misc_t)   :: clubb_misc

     type(pdf_parameter) :: pdf_params    ! PDF parameters (thermo. levs.)   [units vary]
     type(pdf_parameter) :: pdf_params_zm ! PDF parameters on momentum levs. [units vary]

     !---------------------------
     ! For linearize_pbl_winds

     real(core_rknd), pointer :: upwp_sfc_pert    ! u'w' at surface                               [m^2/s^2]
     real(core_rknd), pointer :: vpwp_sfc_pert    ! v'w' at surface                               [m^2/s^2]
     real(core_rknd), pointer, dimension(:) :: um_pert_col ! Pointer to a particular column of um
     real(core_rknd), pointer, dimension(:) :: vm_pert_col
     real(core_rknd), pointer, dimension(:) :: upwp_pert_col
     real(core_rknd), pointer, dimension(:) :: vpwp_pert_col

     !----------------------------------------------------------------------------------------------
     ! l_implemented = .true. is declared and set in clubb_intr.F90
     real(core_rknd), parameter :: host_dx = 100000._core_rknd   ! Host model deltax [m]
     real(core_rknd), parameter :: host_dy = 100000._core_rknd   ! Host model deltay [m]

     !----------------------------------------------------------------------------------------------
     ! The following local variables need to be declared as they are parts of
     ! advance_clubb_core_api's argument list, 
     ! but the input variables are set zero and the output varialbes are not used by the host model.
     !----------------------------------------------------------------------------------------------
     real(core_rknd) ::  rtp3_in(nz)         ! thermodynamic levels (r_t'^3 )               [(kg/kg)^3]
     real(core_rknd) :: thlp3_in(nz)         ! thermodynamic levels (th_l'^3)               [K^3]

     ! hydromet_dim = 0 is declared and set in clubb_intr.F90

     real(core_rknd) :: sclrm        (nz,sclr_dim)  ! Passive scalar mean (thermo. levels)          [units vary]
     real(core_rknd) :: wpsclrp      (nz,sclr_dim)  ! w'sclr' (momentum levels)                     [{units vary} m/s]
     real(core_rknd) :: sclrp2       (nz,sclr_dim)  ! sclr'^2 (momentum levels)                     [{units vary}^2]
     real(core_rknd) :: sclrprtp     (nz,sclr_dim)  ! sclr'rt' (momentum levels)                    [{units vary} (kg/kg)]
     real(core_rknd) :: sclrpthlp    (nz,sclr_dim)  ! sclr'thlp' (momentum levels)                  [{units vary} (K)]
     real(core_rknd) :: sclrpthvp    (nz,sclr_dim)  ! sclr'th_v' (momentum levels)                  [units vary]
     real(core_rknd) :: sclrm_forcing(nz,sclr_dim)  ! Passive scalar forcing                        [{units vary}/s]
     real(core_rknd) :: wpsclrp_sfc  (   sclr_dim)  ! Scalar flux at surface                        [{units vary} m/s]

     real(core_rknd) :: hydromet   (nz,hydromet_dim)
     real(core_rknd) :: wphydrometp(nz,hydromet_dim)
     real(core_rknd) :: wp2hmp     (nz,hydromet_dim)
     real(core_rknd) :: rtphmp_zt  (nz,hydromet_dim)
     real(core_rknd) :: thlphmp_zt (nz,hydromet_dim)

     real(core_rknd) :: radf(nz)               ! input to clubb_core, set to zero 
     real(core_rknd) :: ice_supersat_frac(nz)  ! output from clubb_core, not used.

     integer :: err_code

     !-------------------------------------------------------------------------
     ! rtp3_in and thlp3_in are not currently used in CLUBB's default code.

     rtp3_in(:)  = 0.0_core_rknd
     thlp3_in(:) = 0.0_core_rknd

     ! Higher order scalar inouts, set to zero

               sclrm(:,:) = 0._core_rknd
             wpsclrp(:,:) = 0._core_rknd
              sclrp2(:,:) = 0._core_rknd
            sclrprtp(:,:) = 0._core_rknd
           sclrpthlp(:,:) = 0._core_rknd
           sclrpthvp(:,:) = 0._core_rknd
       sclrm_forcing(:,:) = 0._core_rknd
         wpsclrp_sfc(  :) = 0._core_rknd

     ! Hydrometeors

     hydromet(:,:)       = 0._core_rknd
     wphydrometp(:,:)    = 0._core_rknd
     wp2hmp(:,:)         = 0._core_rknd
     rtphmp_zt(:,:)      = 0._core_rknd
     thlphmp_zt(:,:)     = 0._core_rknd

     radf(:) = 0._core_rknd

         call advance_clubb_core_api &
              ( l_implemented, dtime, fcoriolis, sfc_elevation, &! intent(in)
              hydromet_dim,      & ! intent(in)
              !
              core_forcing%thlm, &
              core_forcing%rtm,  &
              core_forcing%um,   &
              core_forcing%vm,   &
              sclrm_forcing,     &
              core_forcing%edsclr,  &
              core_forcing%wprtp,   &
              core_forcing%wpthlp,  &
              core_forcing%rtp2,    &
              core_forcing%thlp2,   &
              core_forcing%rtpthlp, &
              !
              core_auxil%wm_zm, &! intent(in)
              core_auxil%wm_zt, &! intent(in)
              !
              core_sfc%wpthlp,    &
              core_sfc%wprtp,     &
              core_sfc%upwp,      &
              core_sfc%vpwp,      &! intent(in)
              wpsclrp_sfc,        &
              core_sfc%wpedsclrp, &! intent(in)
              !
              core_auxil%p_in_Pa,         &! intent(in)
              core_auxil%rho_zm,          &
              core_auxil%rho_zt,          &
              core_auxil%exner,           &
              core_auxil%rho_ds_zm,       &
              core_auxil%rho_ds_zt,       &
              core_auxil%invrs_rho_ds_zm, &
              core_auxil%invrs_rho_ds_zt, &
              core_auxil%thv_ds_zm,       &
              core_auxil%thv_ds_zt,       &
              !
              hydromet,         &! intent(in)
              clubb_misc%rfrzm, &! intent(in)
              radf,             &! intent(in)
#ifdef CLUBBND_CAM
              clubb_misc%varmu, &! intent(in)
#endif
              wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &! intent(in)
              host_dx, host_dy, &                          ! intent(in)
              !
              core_prog%um,     &
              core_prog%vm,     &
              core_prog%upwp,   &! intent(inout)
              core_prog%vpwp,   &
              core_prog%up2,    &
              core_prog%vp2,    &! intent(inout)
              core_prog%thlm,   &
              core_prog%rtm,    &
              core_prog%wprtp,  &
              core_prog%wpthlp, &! intent(inout)
              core_prog%wp2,    &
              core_prog%wp3,    &
              core_prog%rtp2,   &! intent(inout)
              rtp3_in,          &
              core_prog%thlp2,  &
              thlp3_in,         &
              core_prog%rtpthlp,&! intent(inout)
              sclrm,                       &! intent(inout)
              sclrp2, sclrprtp, sclrpthlp, &! intent(inout)
              wpsclrp, &
              core_prog%edsclr, &
              err_code,         &! intent(inout)
              !
              core_diag%rcm,        &
              core_diag%cloud_frac, &! intent(inout)
              !
              core_diag%wpthvp,  &
              core_diag%wp2thvp, &
              core_diag%rtpthvp, &
              core_diag%thlpthvp,&! intent(inout)
              !
              sclrpthvp,         &! intent(inout)
              pdf_params,        &! intent(inout)
              pdf_params_zm,     &! intent(inout)
              !
              core_diag%khzm,    &
              core_diag%khzt,    &
              core_diag%qclvar,  &
              core_diag%thlprcp, &! intent(out)
              core_diag%wprcp,   &
              !
              ice_supersat_frac,      &! intent(out)
              core_diag%rcm_in_layer, &! intent(out)
              core_diag%cloud_cover,  &! intent(out)
              !
              upwp_sfc_pert, vpwp_sfc_pert,      &! intent(in)
              um_pert_col, vm_pert_col,          &! 
              upwp_pert_col, vpwp_pert_col)       ! intent(inout)

         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Fatal error in CLUBB: at timestep ", get_nstep(), &
                 "LAT (radians): ",            core_auxil%lat, &
                 "LON (radians): ",            core_auxil%lon, &
                 "LAT (degrees): ", rad_to_deg*core_auxil%lat, &
                 "LON (degrees): ", rad_to_deg*core_auxil%lon, &
                 "Global Column Number: ", get_gcol_p( core_auxil%lchnk,core_auxil%icol )
            call endrun('clubb_tend_cam:  Fatal error in CLUBB library'//errmsg(__FILE__,__LINE__))
         end if

  end subroutine advance_clubb_core_api_eam

