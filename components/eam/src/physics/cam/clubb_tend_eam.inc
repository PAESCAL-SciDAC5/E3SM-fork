  subroutine clubb_tend_eam( ptend_loc, state1, pbuf, cam_in, lq, hdtime,   &! in
                             itim_old, cld_macmic_num_steps, macmic_it,     &! in
                             ixq, ixcldliq, ixcldice,                       &! in
                             do_tms, ksrftms, linearize_pbl_winds,          &! in
                             apply_const,  dz_g, wprcp                      &! out
                            )

   use scamMOD,        only: single_column,scm_clubb_iop_name
   use perf_mod,       only: t_startf, t_stopf
   use physics_types,  only: physics_state
   use physics_types,  only: physics_ptend, physics_ptend_init
   use camsrfexch,     only: cam_in_t
   use physics_buffer, only: physics_buffer_desc

   use clubb_intr_host_types, only: clubb_mean_2d_t, clubb_mnts_2d_t, clubb_to_host_t
   use clubb_intr_core_types
   use clubb_api_module, only:  pdf_parameter
   use clubb_api_module, only:  stats_tsamp, stats_tout, l_stats, stats_begin_timestep_api

   !------------------
   ! Arguments
   !------------------
   type(physics_ptend),intent(inout)  :: ptend_loc
   type(physics_state),intent(inout)  :: state1
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(cam_in_t),      intent(in)    :: cam_in

   real(r8),intent(in) :: hdtime                          ! Host model macmic coupling timestep [s]
   integer, intent(in) :: itim_old
   integer, intent(in) :: cld_macmic_num_steps, macmic_it
   integer, intent(in) :: ixq, ixcldliq, ixcldice

   logical, intent(in) :: lq(:)

   logical, intent(in) :: do_tms
   real(r8),intent(in) :: ksrftms(:)

   logical, intent(in) :: linearize_pbl_winds

   real(r8),intent(out) :: apply_const
   real(r8),intent(out) :: dz_g(:)

   real(r8),intent(out) :: wprcp(pcols,pverp)   ! CLUBB liquid water flux                       [m/s kg/kg]

   real(r8) ::       qclvar(pcols,pverp)        ! cloud water variance [kg^2/kg^2]
   real(r8) :: rcm_in_layer(pcols,pverp)        ! CLUBB in-cloud liquid water mixing ratio      [kg/kg]
   real(r8) ::  cloud_cover(pcols,pverp)        ! CLUBB in-cloud cloud fraction                 [fraction]

   !------------------
   ! Local variables
   !------------------
   integer :: i,ncol, lchnk
   real(r8) :: te_b
   real(r8) :: clubb_s(pcols,pver)

   real(r8) :: varmu(pcols)
   real(r8) :: edsclr_out(pcols,pverp,edsclr_dim)     ! Scalars to be diffused through CLUBB [units vary]

   type(clubb_mean_2d_t) :: host_mean
   type(clubb_mnts_2d_t) :: host_mnts
   type(clubb_to_host_t) :: c2h

   type(core_prog_t )   :: core_prog
   type(core_diag_t )   :: core_diag
   type(core_auxil_t)   :: core_auxil
   type(core_forcing_t) :: core_forcing
   type(core_sfc_t)     :: core_sfc
   type(clubb_misc_t)   :: clubb_misc

   type(pdf_parameter), pointer :: pdf_params    ! PDF parameters (thermo. levs.) [units vary]
   type(pdf_parameter), pointer :: pdf_params_zm ! PDF parameters on momentum levs. [units vary]

   integer :: t   ! loop index of CLUBB's timestep loop
   integer :: nz

   real(core_rknd) :: hdtime_core_rknd          ! host model's cloud macmic timestep in core_rknd [s]
   real(core_rknd) :: core_dtime                ! CLUBB time step [s]
   integer         :: n_clubb_core_step

   real(core_rknd) :: core_rknd_landfrac
   real(core_rknd) :: core_rknd_rnevap_effic

   real(core_rknd) :: zt_bot
   real(core_rknd) :: fcoriolis           ! Coriolis forcing [s^-1]
   real(core_rknd) :: sfc_elevation       ! Elevation of ground [m AMSL]

   integer :: time_elapsed                ! time keep track of stats [s]

   logical :: l_pbl_winds_diag
   real(r8) :: sfc_v_diff_tau(pcols)      ! Response to tau perturbation, m/s
   real(r8), parameter :: pert_tau = 0.1_r8 ! tau perturbation, Pa

   real(r8) :: zt_out(pcols,pverp)              ! output for the thermo CLUBB grid              [m]
   real(r8) :: zi_out(pcols,pverp)              ! output for momentum CLUBB grid                [m]

   !------------------
   ncol  = state1%ncol
   lchnk = state1%lchnk

   !-----------------------------
   ! Initialize ptend for CLUBB
   !-----------------------------
   call physics_ptend_init(ptend_loc, state1%psetcols, 'clubb', ls=.true., lu=.true., lv=.true., lq=lq)

   !----------------------------------------------------------------------------------
   ! Load input (initial and boundary conditions) from host model's data structures
   !----------------------------------------------------------------------------------
   call gather_host_fields( cld_macmic_num_steps, macmic_it, itim_old, state1, &! in
                            pbuf, cam_in, host_mean, host_mnts, apply_const    )! inout, out

   call setup_clubb_to_host_fields( pbuf, c2h, itim_old, linearize_pbl_winds )

   !===========================================================================================
   !  Loop over all columns in lchnk to advance CLUBB core
   !===========================================================================================
   call t_startf('adv_clubb_core_col_loop')

   do i=1,ncol   ! loop over columns

      !-----------------------------------------------------------------
      ! Calculate column total energy before CLUBB. The value is saved
      ! for the energy fixer at the end of the column loop.
      !-----------------------------------------------------------------
      call column_total_energy( state1, cam_in, i, hdtime, &! in
                                te_b                       )! out, scalar

      !------------------------------------------------------------------------------------
      ! CLUBB can be substepped with respect to the host model's timestep hdtime.
      ! Determine the substepsize and number of substeps for CLUBB.
      ! Note that in the standard EAM which uses the same vertical grid for all processes,
      ! this call only needs to be done once at the beginning of an integration.
      ! But we place the call here (meaning repeating the calculation every timestep and
      ! in each grid column), as we are preparing for the implementation of AMR.
      !------------------------------------------------------------------------------------
      call determine_clubb_dtime( hdtime, clubb_timestep, hdtime_core_rknd, core_dtime, n_clubb_core_step ) ! 2xin, 3xout

      ! Set time_elapsed to host model time step, this is for CLUBB's budget stats
      time_elapsed = hdtime_core_rknd

      ! Set stats output and increment equal to CLUBB and host dt
      stats_tsamp = core_dtime
      stats_tout  = hdtime_core_rknd

      !------------------------------------------------------------------------------------
      ! Setup CLUBB's z-grid for this column of the host model
      !------------------------------------------------------------------------------------
      nz = pverp
      call setup_clubb_zgrid_for_1column( nz, pverp, host_mean%zi(i,:), host_mean%zm(i,:),   &
                                          zt_bot, sfc_elevation, zt_out(i,:), zi_out(i,:)    )

      !--------------------------------------
      ! Prepare input to advance_clubb_core
      !--------------------------------------
      ! Allocate memory for derived-type variables

      call clubb_core_fld_alloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc, &! inout
                                 nz, edsclr_dim, linearize_pbl_winds )! in

      ! Associate pointers for CLUBB's PDF parameters

      pdf_params    => pdf_params_chnk(i,lchnk)
      pdf_params_zm => pdf_params_zm_chnk(i,lchnk)

      ! Now fill those data structures with values

      call map_1column_of_info_from_host_to_clubb( i, host_mean, host_mnts, c2h,         &! in
                                                   pverp, nz, edsclr_dim, lq,            &! in
                                                   ixq, ixcldliq, ixcldice,              &! in
                                                   macmic_it, clubb_do_adv, do_expldiff, &! in
                                                   linearize_pbl_winds, pert_tau,        &! in
                                                   core_auxil, core_prog, core_diag,     &! inout
                                                   core_forcing, core_sfc, clubb_misc,   &! inout
                                                   pdf_params_zm,                        &! inout
                                                   dz_g, varmu(i)                        )! out


      !**************************************************************************************
      ! CLUBB's time-stepping sub-cycles
      !**************************************************************************************
      call t_startf('adv_clubb_core_ts_loop')
      do t=1,n_clubb_core_step    ! do needed number of "sub" timesteps for each macmic step

         !  Increment the statistics then being stats timestep
         if (l_stats) then
            time_elapsed = time_elapsed+core_dtime
            call stats_begin_timestep_api(time_elapsed, 1, 1)
         endif

         !--------------------------------------
         ! Set forcing terms and surface values
         !--------------------------------------
         call set_clubb_forcing_and_sfc( i, pverp, nz, host_mean, core_auxil, &! in
                                         do_tms, ksrftms,                     &! in
                                         core_forcing, core_sfc, fcoriolis    )! inout or out

         if (single_column) call set_clubb_sfc_for_scm( scm_clubb_iop_name, &! in
                                                        i, pver, host_mean, &! in
                                                        zt_bot, theta0,     &! in
                                                        core_sfc            )! inout

         !--------------------------------------
         ! Advance CLUBB CORE for one timestep
         !--------------------------------------
         call t_startf('advance_clubb_core')
         call advance_clubb_core_api_eam( nz, core_dtime, sfc_elevation, fcoriolis, &! in
                                          core_forcing, core_sfc,     &! in
                                          clubb_misc,   core_auxil,   &! in
                                          core_prog,    core_diag,    &! inout or out
                                          pdf_params,   pdf_params_zm &! inout or out
                                          )
         call t_stopf('advance_clubb_core')

         !------------------------------------------------------------
         ! Consider impact of rain evaporation on turbulent moments 
         !------------------------------------------------------------
         if (do_rainturb) then

            core_rknd_landfrac     = real( cam_in%landfrac(i), kind=core_rknd )
            core_rknd_rnevap_effic = real( clubb_rnevap_effic, kind=core_rknd )

            call update_mnts_for_rain_evap( core_prog%rtp2,       core_prog%thlp2,        &! inout
                                            core_prog%wprtp,      core_prog%wpthlp,       &! inout
                                            core_prog%rtm,        core_diag%rcm,          &! in
                                            core_diag%cloud_frac, core_prog%thlm,         &! in
                                            core_auxil%wm_zt,     core_auxil%exner,       &! in
                                            clubb_misc%prer_evap, pdf_params,             &! in
                                            core_rknd_landfrac,   core_rknd_rnevap_effic, &! in
                                            nz, core_dtime, clubb_do_deep                 )! in
         endif

         !-----------------------------------------------
         ! Consider impact of cloud-top cooling on thlp2
         !-----------------------------------------------
         if (do_cldcool) then
            call enchance_thlp2_for_cloud_top_cooling( core_prog%thlp2,                  &! inout
                                                       core_diag%rcm, clubb_misc%qrl_zt, &! in
                                                       core_diag%thlprcp,                &! in
                                                       nz, core_dtime                    )
         endif

         !------------------------------------------------------------------- 
         !  Check to see if stats should be output, here stats are read into
         !  output arrays to make them conformable to CAM output
         !------------------------------------------------------------------- 
         if (l_stats) call stats_end_timestep_clubb(lchnk,i,out_zt,out_zm, out_radzt,out_radzm,out_sfc)

      enddo  ! end time loop
      call t_stopf('adv_clubb_core_ts_loop')
      !**************************************************************************************
      ! End of CLUBB's time-stepping sub-cycles
      !**************************************************************************************

      !-------------------------------------------------------------------------------
      ! Transfer information from local column arrays to host model's data structures
      !-------------------------------------------------------------------------------
      call map_1column_of_info_from_clubb_to_host( i, pverp, nz, edsclr_dim,                      &! in
                                                   macmic_it, cld_macmic_num_steps, clubb_do_adv, &! in
                                                   linearize_pbl_winds,                           &! in 
                                                   core_prog, core_diag, pdf_params_zm, clubb_misc, &! inout
                                                   host_mean, host_mnts, c2h, edsclr_out,         &! inout
                                                   wprcp, rcm_in_layer, cloud_cover, qclvar,      &! inout
                                                   sfc_v_diff_tau                                 )! inout
      ! Clean up

      call clubb_core_fld_dealloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc, linearize_pbl_winds )

      !---------------------------------------------------------------------------------------------------------
      ! After advancing CLUBB in time, apply a total energy fixer by adjusting the ending s (dry static energy)
      !---------------------------------------------------------------------------------------------------------
      call column_total_energy_fixer( te_b, host_mean%thlm, host_mean%rtm, host_mean%rcm, host_mean%exner_clubb, &! in
                                      host_mean%um, host_mean%vm, host_mnts%wp2, state1, i,                      &! in
                                      clubb_s(i,:)                                                               )! out

   enddo  ! end column loop
   call t_stopf('adv_clubb_core_col_loop')
   !===========================================================================================
   ! End of column loop
   !===========================================================================================

   !-------------------------------------------------------------
   ! Calculate ptend for CLUBB
   !-------------------------------------------------------------
   call derive_ptend( host_mean, host_mnts, state1, clubb_s, edsclr_out, hdtime,    &! in
                      ixq, ixcldliq, clubb_do_adv, macmic_it, cld_macmic_num_steps, &! in
                      ptend_loc )! inout

   deallocate( host_mean%exner_clubb )

   !-------------------------------------------------------------
   ! Diagnostics and output
   !-------------------------------------------------------------
   l_pbl_winds_diag = linearize_pbl_winds .and. (macmic_it == cld_macmic_num_steps)
   call clubb_related_diags_on_host_grid( pbuf, ptend_loc, host_mean, host_mnts, c2h,      &! in
                                          qclvar, deep_scheme, clubb_do_deep, varmu,       &! in
                                          l_pbl_winds_diag, sfc_v_diff_tau, pert_tau,      &! in
                                          apply_const, zi_out, zt_out, rcm_in_layer, wprcp )

   !  Output CLUBB's statistics
   call clubb_stats_output(pcols,lchnk)

  end subroutine clubb_tend_eam
