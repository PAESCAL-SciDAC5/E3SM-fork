
      zi_out(i,1) = 0._r8

      !======================================
      ! Prepare input to advance_clubb_core
      !======================================
      ! Allocate memory for derived-type variables

      call clubb_core_fld_alloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc, pverp, edsclr_dim )

      ! Associate pointers for CLUBB's PDF parameters

      pdf_params    => pdf_params_chnk(i,lchnk)
      pdf_params_zm => pdf_params_zm_chnk(i,lchnk)

      ! Now fill those data structures with values

      call map_1column_of_info_from_host_to_clubb( i, host_mean, host_mnts, c2h,         &! in
                                                   pver, pverp, edsclr_dim, lq,          &! in
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

         !  Advance CLUBB CORE one timestep in the future
         call t_startf('advance_clubb_core')

#include "clubb_timestep_init_forcing_and_sfc.inc"

#include "advance_clubb_core_api_eam.inc"

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
                                            pverp, dtime, clubb_do_deep                   )! in
         endif

         !-----------------------------------------------
         ! Consider impact of cloud-top cooling on thlp2
         !-----------------------------------------------
         if (do_cldcool) then
            call enchance_thlp2_for_cloud_top_cooling( core_prog%thlp2,                  &! inout
                                                       core_diag%rcm, clubb_misc%qrl_zt, &! in
                                                       core_diag%thlprcp,                &! in
                                                       pverp, dtime                      )
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


      !=================================================
      call clubb_core_fld_dealloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc )

