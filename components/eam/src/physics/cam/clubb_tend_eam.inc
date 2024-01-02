
   !-----------------------------
   ! Initialize ptend for CLUBB
   !-----------------------------
   call physics_ptend_init(ptend_loc, state1%psetcols, 'clubb', ls=.true., lu=.true., lv=.true., lq=lq)

   !-------------------------------------------------------------
   ! Load input (initial and boundary conditions) from host model's data structures
   !-------------------------------------------------------------
   call gather_host_fields( cld_macmic_num_steps, macmic_it, itim_old, state1, &! in
                            pbuf, cam_in, host_mean, host_mnts, apply_const    )! inout, out

   call setup_clubb_to_host_fields( pbuf, c2h, itim_old )

#include "clubb_misc_inout_fields.inc"

   call t_startf('adv_clubb_core_col_loop')
   !----------------------------------------------------------
   !  Loop over all columns in lchnk to advance CLUBB core
   !----------------------------------------------------------
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
      call determine_clubb_dtime( hdtime, clubb_timestep, hdtime_core_rknd, dtime, n_clubb_core_step ) ! 2xin, 3xout

      ! Set time_elapsed to host model time step, this is for CLUBB's budget stats
      time_elapsed = hdtime_core_rknd

      ! Set stats output and increment equal to CLUBB and host dt
      stats_tsamp = dtime
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

      call clubb_core_fld_alloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc, nz, edsclr_dim )

      ! Associate pointers for CLUBB's PDF parameters

      pdf_params    => pdf_params_chnk(i,lchnk)
      pdf_params_zm => pdf_params_zm_chnk(i,lchnk)

      ! Now fill those data structures with values

      call map_1column_of_info_from_host_to_clubb( i, host_mean, host_mnts, c2h,         &! in
                                                   pverp, nz, edsclr_dim, lq,            &! in
                                                   ixq, ixcldliq, ixcldice,              &! in
                                                   macmic_it, clubb_do_adv, do_expldiff, &! in
                                                   core_auxil, core_prog, core_diag,     &! inout
                                                   core_forcing, core_sfc, clubb_misc,   &! inout
                                                   pdf_params_zm,                        &! inout
                                                   dz_g, varmu(i)                        )! out
#include "more_info_to_clubb.inc"


      !**************************************************************************************
      ! CLUBB's time-stepping sub-cycles
      !**************************************************************************************
      call t_startf('adv_clubb_core_ts_loop')
      do t=1,n_clubb_core_step    ! do needed number of "sub" timesteps for each CAM step

         !  Increment the statistics then being stats timestep
         if (l_stats) then
            time_elapsed = time_elapsed+dtime
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
         call advance_clubb_core_api_eam( nz, dtime, sfc_elevation, fcoriolis, &! in
                                          core_forcing, core_sfc,     &! in
                                          clubb_misc,   core_auxil,   &! in
                                          core_prog,    core_diag,    &! inout or out
                                          pdf_params,   pdf_params_zm,&! inout or out
                                          !
                                          upwp_sfc_pert, vpwp_sfc_pert,      &! intent(in)
                                          um_pert_col, vm_pert_col,          &!
                                          upwp_pert_col, vpwp_pert_col       &! intent(inout)
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
                                            nz, dtime, clubb_do_deep                      )! in
         endif

         !-----------------------------------------------
         ! Consider impact of cloud-top cooling on thlp2
         !-----------------------------------------------
         if (do_cldcool) then
            call enchance_thlp2_for_cloud_top_cooling( core_prog%thlp2,                  &! inout
                                                       core_diag%rcm, clubb_misc%qrl_zt, &! in
                                                       core_diag%thlprcp,                &! in
                                                       nz, dtime                         )
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

#include "info_from_clubb_to_host.inc"

      ! Clean up
      call clubb_core_fld_dealloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc )

      !---------------------------------------------------------------------------------------------------------
      ! After advancing CLUBB in time, apply a total energy fixer by adjusting the ending s (dry static energy)
      !---------------------------------------------------------------------------------------------------------
      call column_total_energy_fixer( te_b, host_mean%thlm, host_mean%rtm, host_mean%rcm, host_mean%exner_clubb, &! in
                                      host_mean%um, host_mean%vm, host_mnts%wp2, state1, i,                      &! in
                                      clubb_s(i,:)                                                               )! out

   enddo  ! end column loop
   call t_stopf('adv_clubb_core_col_loop')

   !-------------------------------------------------------------
   ! Calculate ptend for CLUBB
   !-------------------------------------------------------------
   call derive_ptend( host_mean, host_mnts, state1, clubb_s, edsclr_out, hdtime,    &! in
                      ixq, ixcldliq, clubb_do_adv, macmic_it, cld_macmic_num_steps, &! in
                      ptend_loc )! inout

   deallocate( host_mean%exner_clubb )

   !===================
   ! History output
   !===================
   call outfld( 'RVMTEND_CLUBB', ptend_loc%q(:,:,ixq), pcols, lchnk)
   call outfld( 'RCMTEND_CLUBB', ptend_loc%q(:,:,ixcldliq), pcols, lchnk)
   call outfld( 'RIMTEND_CLUBB', ptend_loc%q(:,:,ixcldice), pcols, lchnk)
   call outfld( 'TTEND_CLUBB',   ptend_loc%s/cpair,pcols, lchnk)
   call outfld( 'UTEND_CLUBB',   ptend_loc%u,pcols, lchnk)
   call outfld( 'VTEND_CLUBB',   ptend_loc%v,pcols, lchnk)

   call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)
   cmeliq(:,:) = ptend_loc%q(:,:,ixcldliq)
   call outfld( 'CMELIQ',        cmeliq, pcols, lchnk)

   if (clubb_do_deep) call outfld( 'MU_CLUBB',      varmu      ,pcols, lchnk)

   ! Linearization of PBL winds: diagnose wsresp and tau_est for history output.
   ! The outfld calls are in cam_diagnostics.F90

   if (linearize_pbl_winds .and. macmic_it == cld_macmic_num_steps) then
      do i = 1, ncol
         wsresp(i) = sfc_v_diff_tau(i) / pert_tau
         ! Estimated tau in balance with wind is the tau we just used.
         if (cam_in%wsx(i) == 0._r8 .or. cam_in%wsy(i) == 0._r8) then
            ! Work around an odd FPE issue with intel compiler.
            tau_est(i) = abs(cam_in%wsx(i)) + abs(cam_in%wsy(i))
         else
            tau_est(i) = hypot(cam_in%wsx(i), cam_in%wsy(i))
         end if
      end do
   end if

#include "clubb_stats_output.inc"

