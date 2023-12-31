
      zi_out(i,1) = 0._r8

      call clubb_core_fld_alloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, clubb_misc, pverp, edsclr_dim )

      !============================================
      ! Prepare input to advance_clubb_core
      !============================================
      call map_1column_of_info_from_host_to_clubb( i, host_mean, host_mnts, c2h,               &
                                                   pver, pverp, edsclr_dim, lq, ixq, ixcldliq, &
                                                   macmic_it, clubb_do_adv, do_expldiff,       &
                                                   core_auxil, core_prog, core_diag,           &
                                                   core_forcing, core_sfc, clubb_misc, dz_g    )

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

         !============================================
         pdf_zm_w_1_inout        = pdf_params_zm%w_1
         pdf_zm_w_2_inout        = pdf_params_zm%w_2
         pdf_zm_varnce_w_1_inout = pdf_params_zm%varnce_w_1
         pdf_zm_varnce_w_2_inout = pdf_params_zm%varnce_w_2
         pdf_zm_mixt_frac_inout  = pdf_params_zm%mixt_frac

         if (do_rainturb) then

            do k=2,pverp
               pre_in(k)    = real(prer_evap(i,pverp-k+1), kind = core_rknd)
            enddo
            pre_in(1) = pre_in(2)

            rvm_in = core_prog%rtm - core_diag%rcm
            call update_xp2_mc_api(pverp, dtime, core_diag%cloud_frac, &
                                   core_diag%rcm, rvm_in, core_prog%thlm, &
                                   core_auxil%wm_zt, core_auxil%exner, pre_in, pdf_params, &
            rtp2_mc_out, thlp2_mc_out, &
            wprtp_mc_out, wpthlp_mc_out, &
            rtpthlp_mc_out)

            if (clubb_do_deep) then
               dum_core_rknd = 1._core_rknd
            else
               dum_core_rknd = (1._core_rknd - real(cam_in%landfrac(i), kind = core_rknd))
            end if

            ! update turbulent moments based on rain evaporation
            core_prog%rtp2  = core_prog%rtp2  + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd *  rtp2_mc_out * dtime
            core_prog%thlp2 = core_prog%thlp2 + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * thlp2_mc_out * dtime
            if (.not. clubb_do_deep) then
               core_prog%wprtp  = core_prog%wprtp  + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd *  wprtp_mc_out * dtime
               core_prog%wpthlp = core_prog%wpthlp + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * wpthlp_mc_out * dtime
            endif
!                     rtpthlp_in = rtpthlp_in + rtpthlp_mc_out * dtime

         endif

         !----------------------
         if (do_cldcool) then
            call enchance_thlp2_for_cloud_top_cooling( core_diag%rcm, clubb_misc%qrl_zt, core_diag%thlprcp, &
                                                       core_prog%thlp2, pverp, dtime )
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

