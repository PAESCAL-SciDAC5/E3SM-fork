
      zi_out(i,1) = 0._r8

      call clubb_core_state_alloc( core_state, pverp )

      !============================================
      ! Prepare input to advance_clubb_core
      !============================================
      invrs_gravit = 1._r8 / gravit

      do k=1,pver
         dz_g(k) = host_mean%zi(i,k)-host_mean%zi(i,k+1)  ! compute thickness
      enddo

      ! k = pver
      ! dz_g_bot = host_mean%zi(i,k)-host_mean%zi(i,k+1)  ! compute thickness

      !  Compute thermodynamic stuff needed for CLUBB on thermo levels.
      !  Inputs for the momentum levels are set below setup_clubb core
      do k=1,pver
         kkhost = pver-k+1

         core_state%p_in_Pa(k+1)   = real( host_mean%pmid(i,kkhost), kind=core_rknd )
         core_state%exner(k+1)     = real( 1._r8/host_mean%exner_clubb(i,kkhost), kind=core_rknd )

         thvtmp =  host_mean%t(i,kkhost) &
                  *host_mean%exner_clubb(i,kkhost) &
                  *( 1._r8+zvir*host_mean%q(i,kkhost,ixq)-host_mean%q(i,kkhost,ixcldliq) )
         core_state%thv_ds_zt(k+1) = real(thvtmp, kind = core_rknd)                       ! thetav on thermo
      enddo

      core_state%p_in_Pa  (1) = core_state%p_in_Pa(2)
      core_state%exner    (1) = core_state%exner  (2)
      core_state%thv_ds_zt(1) = core_state%thv_ds_zt(2)

      !----------
      do k=1,pver
         kkhost = pver-k+1
         core_state%rho_ds_zt(k+1) = real( invrs_gravit*host_mean%pdel(i,kkhost)/dz_g(kkhost), kind=core_rknd )
      enddo
      core_state%rho_ds_zt(1) = core_state%rho_ds_zt(2)


      !----------
      core_state%rho_zt(:)  = core_state%rho_ds_zt(:)

      core_state%invrs_rho_ds_zt(2:pverp) = 1._core_rknd/core_state%rho_ds_zt(2:pverp) ! Inverse ds rho at thermo
      core_state%invrs_rho_ds_zt(1)       = core_state%invrs_rho_ds_zt(2)

      !  Compute some inputs from the thermodynamic grid to the momentum grid
      !  Note: this can only be done after "call setup_grid_heights_api"

      core_state%rho_ds_zm       = zt2zm_api( core_state%rho_ds_zt )
      core_state%rho_zm          = zt2zm_api( core_state%rho_zt    )
      core_state%invrs_rho_ds_zm = zt2zm_api( core_state%invrs_rho_ds_zt)
      core_state%thv_ds_zm       = zt2zm_api( core_state%thv_ds_zt )

      !----------
      !  Compute mean w wind on thermo grid, convert from omega to w

      core_state%wm_zt(1) = 0._core_rknd
      do k=1,pver
        kkhost = pver-k+1
        dum1 = -1._r8*host_mean%omega(i,kkhost) * real(core_state%invrs_rho_ds_zt(k+1),kind=r8) *invrs_gravit
        core_state%wm_zt(k+1) = real(dum1,kind=core_rknd)
      enddo

      core_state%wm_zm = zt2zm_api( core_state%wm_zt )

      !----------
      do k=1,pver
         kkhost = pver-k+1
         rfrzm(k+1)     = real( host_mean%q(i,kkhost,ixcldice),                 kind=core_rknd )
         radf(k+1)      = real( radf_clubb(i,kkhost),                           kind=core_rknd )
      enddo
      radf(1)      = radf(2)
      rfrzm(1)     = rfrzm(2)
      !----------

      !  Need to flip arrays around for CLUBB core
      do k=1,pverp

         core_state%thlm(k) = real(host_mean%thlm(i,pverp-k+1), kind = core_rknd)
         core_state%rtm (k) = real(host_mean%rtm(i,pverp-k+1), kind = core_rknd)
         core_state%rcm (k) = real(host_mean%rcm(i,pverp-k+1), kind = core_rknd)

         core_state%um(k)      = real(host_mean%um(i,pverp-k+1), kind = core_rknd)
         core_state%vm(k)      = real(host_mean%vm(i,pverp-k+1), kind = core_rknd)

         upwp_in(k)    = real(host_mnts%upwp(i,pverp-k+1), kind = core_rknd)
         vpwp_in(k)    = real(host_mnts%vpwp(i,pverp-k+1), kind = core_rknd)
         up2_in(k)     = real(host_mnts%up2(i,pverp-k+1), kind = core_rknd)
         vp2_in(k)     = real(host_mnts%vp2(i,pverp-k+1), kind = core_rknd)
         wp2_in(k)     = real(host_mnts%wp2(i,pverp-k+1), kind = core_rknd)
         wp3_in(k)     = real(host_mnts%wp3(i,pverp-k+1), kind = core_rknd)
         rtp2_in(k)    = real(host_mnts%rtp2(i,pverp-k+1), kind = core_rknd)
         thlp2_in(k)   = real(host_mnts%thlp2(i,pverp-k+1), kind = core_rknd)

         wprtp_in(k)   = real(host_mnts%wprtp(i,pverp-k+1), kind = core_rknd)
         wpthlp_in(k)  = real(host_mnts%wpthlp(i,pverp-k+1), kind = core_rknd)
         rtpthlp_in(k) = real(host_mnts%rtpthlp(i,pverp-k+1), kind = core_rknd)

         wpthvp_inout(k)     = real(host_mnts%wpthvp(i,pverp-k+1), kind = core_rknd)
         wp2thvp_inout(k)    = real(host_mnts%wp2thvp(i,pverp-k+1), kind = core_rknd)
         rtpthvp_inout(k)    = real(host_mnts%rtpthvp(i,pverp-k+1), kind = core_rknd)
         thlpthvp_inout(k)   = real(host_mnts%thlpthvp(i,pverp-k+1), kind = core_rknd)
         cloud_frac_inout(k) = real(cloud_frac(i,pverp-k+1), kind = core_rknd)

         ! also flip the arrays for the pdf_params_zm variables to have
         ! consistent vertical orientation for all variables in restart file
         ! also need to flip back after calling advance_clubb_core

         pdf_zm_w_1_inout(k) = pdf_zm_w_1(i,pverp-k+1)
         pdf_zm_w_2_inout(k) = pdf_zm_w_2(i,pverp-k+1)
         pdf_zm_varnce_w_1_inout(k) = pdf_zm_varnce_w_1(i,pverp-k+1)
         pdf_zm_varnce_w_2_inout(k) = pdf_zm_varnce_w_2(i,pverp-k+1)
         pdf_zm_mixt_frac_inout(k) =  pdf_zm_mixt_frac(i,pverp-k+1)
      enddo


      if (linearize_pbl_winds) then
         ! Each host model time step, reset the perturbed variables to be equal to
         ! the unperturbed values.
         if (macmic_it == 1) then
            um_pert(i,:) = core_state%um
            vm_pert(i,:) = core_state%vm
            upwp_pert(i,:) = upwp_in
            vpwp_pert(i,:) = vpwp_in
         end if
         allocate(um_pert_col(pverp))
         allocate(vm_pert_col(pverp))
         allocate(upwp_pert_col(pverp))
         allocate(vpwp_pert_col(pverp))
         um_pert_col = um_pert(i,:)
         vm_pert_col = vm_pert(i,:)
         upwp_pert_col = upwp_pert(i,:)
         vpwp_pert_col = vpwp_pert(i,:)

         allocate(upwp_sfc_pert)
         allocate(vpwp_sfc_pert)
         ! Prefer to perturb wind/stress in the direction of the existing stress.
         ! However, if there's no existing surface stress, just perturb zonal
         ! wind/stress.
         if (abs(cam_in%wsx(i)) < 1.e-12 .and. abs(cam_in%wsy(i)) < 1.e-12) then
            upwp_sfc_pert = upwp_sfc + pert_tau / core_state%rho_ds_zm(1)
            vpwp_sfc_pert = vpwp_sfc
         else
            upwp_sfc_pert = upwp_sfc + cam_in%wsx(i) * &
                 (pert_tau / (core_state%rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
            vpwp_sfc_pert = vpwp_sfc + cam_in%wsy(i) * &
                 (pert_tau / (core_state%rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
         end if
      else
         nullify(upwp_sfc_pert)
         nullify(vpwp_sfc_pert)
         nullify(um_pert_col)
         nullify(vm_pert_col)
         nullify(upwp_pert_col)
         nullify(vpwp_pert_col)
      end if


      if (clubb_do_adv) then
        if (macmic_it .eq. 1) then
          wp2_in=zt2zm_api(wp2_in)
          wpthlp_in=zt2zm_api(wpthlp_in)
          wprtp_in=zt2zm_api(wprtp_in)
          up2_in=zt2zm_api(up2_in)
          vp2_in=zt2zm_api(vp2_in)
          thlp2_in=zt2zm_api(thlp2_in)
          rtp2_in=zt2zm_api(rtp2_in)
          rtpthlp_in=zt2zm_api(rtpthlp_in)

          do k=1,pverp
            thlp2_in(k)=max(thl_tol**2,thlp2_in(k))
            rtp2_in(k)=max(rt_tol**2,rtp2_in(k))
            wp2_in(k)=max(w_tol_sqd,wp2_in(k))
            up2_in(k)=max(w_tol_sqd,up2_in(k))
            vp2_in(k)=max(w_tol_sqd,vp2_in(k))
          enddo
        endif
      endif

      !-----------------------------------------------
      ! Tracers 
      !-----------------------------------------------
      !  Initialize these to prevent crashing behavior

      edsclr_in(:,:)      = 0._core_rknd
      edsclr_out(:,:)     = 0._r8

      !  Get input values for tracers
 
      icnt=0
      do ixind=1,pcnst
         if (lq(ixind))  then
            icnt=icnt+1
            do k=1,pver
               edsclr_in(k+1,icnt) = real(host_mean%q(i,pver-k+1,ixind), kind = core_rknd)
            enddo
            edsclr_in(1,icnt) = edsclr_in(2,icnt)
         end if
      enddo

      if (do_expldiff) then
        do k=1,pver
          edsclr_in(k+1,icnt+1) = real(host_mean%thlm(i,pver-k+1), kind = core_rknd)
          edsclr_in(k+1,icnt+2) = real(host_mean%rtm(i,pver-k+1), kind = core_rknd)
        enddo

        edsclr_in(1,icnt+1) = edsclr_in(2,icnt+1)
        edsclr_in(1,icnt+2) = edsclr_in(2,icnt+2)
      endif

      ! --------------------------------------------------------- !
      ! Compute cloud-top radiative cooling contribution to CLUBB !
      ! --------------------------------------------------------- !
      ! Sandbox version of code to take into account meso organization

      if (clubb_do_deep) then
         orgparam = 0._r8
         delpavg = 0._r8

         do k = 1, pver
           if (abs(prer_evap(i,k)) .gt. 0._r8) then
             orgparam = orgparam + (abs(prer_evap(i,k)) * 1000._r8 * 1000._r8 * 2._r8 ) * host_mean%pdel(i,k)
             delpavg = delpavg + host_mean%pdel(i,k)
           endif
         enddo

         if (delpavg .gt. 0._r8) then
           orgparam = orgparam/delpavg
         endif

         ! Now compute new entrainment rate based on organization
         varmu(i) = mu / (1._r8 + orgparam * 100._r8)
         varmu2   = real(varmu(i), kind = core_rknd)

      endif
      ! --------------------------------------------------------- !
      ! End cloud-top radiative cooling contribution to CLUBB     !
      ! --------------------------------------------------------- !

      pdf_params    => pdf_params_chnk(i,lchnk)
      pdf_params_zm => pdf_params_zm_chnk(i,lchnk)

      if ( is_first_restart_step() .and. ipdf_call_placement .eq. ipdf_post_advance_fields ) then
         ! assign the values read back from restart file
         ! This is necessary when ipdf_call_placement = 2
         pdf_params_zm%w_1 = pdf_zm_w_1_inout
         pdf_params_zm%w_2 = pdf_zm_w_2_inout
         pdf_params_zm%varnce_w_1 = pdf_zm_varnce_w_1_inout
         pdf_params_zm%varnce_w_2 = pdf_zm_varnce_w_2_inout
         pdf_params_zm%mixt_frac = pdf_zm_mixt_frac_inout
      end if
      !============================================


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

#include "clubb_timestep_init_forcing_and_boundary.inc"

         !  Advance CLUBB CORE one timestep in the future
         call t_startf('advance_clubb_core')
         !Balli- to do: check whether initent-ins and intent-inouts are actually what they say

         call advance_clubb_core_api &
              ( l_implemented, dtime, fcoriolis, sfc_elevation, hydromet_dim, & ! intent(in)
              thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &         ! intent(in)
              sclrm_forcing, edsclrm_forcing, wprtp_forcing, &             ! intent(in)
              wpthlp_forcing, rtp2_forcing, thlp2_forcing, &               ! intent(in)
              rtpthlp_forcing,  &
              !
              core_state%wm_zm, &
              core_state%wm_zt, &                             ! intent(in)
              !
              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &                 ! intent(in)
              wpsclrp_sfc, wpedsclrp_sfc, &                                ! intent(in)
              !
              core_state%p_in_Pa,         &
              core_state%rho_zm,          &
              core_state%rho_zt,          &
              core_state%exner,           &                          ! intent(in)
              core_state%rho_ds_zm,       &
              core_state%rho_ds_zt,       &
              core_state%invrs_rho_ds_zm, &                          ! intent(in)
              core_state%invrs_rho_ds_zt, &
              core_state%thv_ds_zm,       &
              core_state%thv_ds_zt,       &
              !
              hydromet, &           ! intent(in)
              rfrzm, radf, &                                               ! intent(in)
#ifdef CLUBBND_CAM
              varmu2, &                                                    ! intent(in)
#endif
              wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &                ! intent(in)
              host_dx, host_dy, &                                          ! intent(in)
              core_state%um, &
              core_state%vm, &
              upwp_in, &                                     ! intent(inout)
              vpwp_in, up2_in, vp2_in, &                                   ! intent(inout)
              core_state%thlm, &
              core_state%rtm,  &
              wprtp_in, wpthlp_in, &                      ! intent(inout)
              wp2_in, wp3_in, rtp2_in, &                                   ! intent(inout)
              rtp3_in, thlp2_in, thlp3_in, rtpthlp_in, &                   ! intent(inout)
              sclrm,   &                                                   ! intent(inout)
              sclrp2, sclrprtp, sclrpthlp, &                               ! intent(inout)
              wpsclrp, edsclr_in, err_code, &                              ! intent(inout)
              !
              core_state%rcm, &
              !
              cloud_frac_inout, &                               ! intent(inout)
              wpthvp_inout, wp2thvp_inout, rtpthvp_inout, thlpthvp_inout, & ! intent(inout)
              sclrpthvp_inout, &                                            ! intent(inout)
              pdf_params, pdf_params_zm, &                                 ! intent(inout)
              khzm_out, khzt_out, qclvar_out, thlprcp_out, &               ! intent(out)
              wprcp_out, ice_supersat_frac, &                              ! intent(out)
              rcm_in_layer_out, cloud_cover_out, &                         ! intent(out)
              upwp_sfc_pert, vpwp_sfc_pert, &                              ! intent(in)
              um_pert_col, vm_pert_col, upwp_pert_col, vpwp_pert_col)      ! intent(inout)
         call t_stopf('advance_clubb_core')

         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Fatal error in CLUBB: at timestep ", get_nstep(), &
                 "LAT (radians): ", host_mean%lat(i), &
                 "LON (radians): ", host_mean%lon(i), &
                 "LAT (degrees): ", rad_to_deg*host_mean%lat(i), &
                 "LON (degrees): ", rad_to_deg*host_mean%lon(i), &
                 "Global Column Number: ", get_gcol_p(lchnk,i)
            call endrun('clubb_tend_cam:  Fatal error in CLUBB library'//errmsg(__FILE__,__LINE__))
         end if

         !============================================
         pdf_zm_w_1_inout = pdf_params_zm%w_1
         pdf_zm_w_2_inout = pdf_params_zm%w_2
         pdf_zm_varnce_w_1_inout = pdf_params_zm%varnce_w_1
         pdf_zm_varnce_w_2_inout = pdf_params_zm%varnce_w_2
         pdf_zm_mixt_frac_inout = pdf_params_zm%mixt_frac

         if (do_rainturb) then

            do k=2,pverp
               pre_in(k)    = real(prer_evap(i,pverp-k+1), kind = core_rknd)
            enddo
            pre_in(1) = pre_in(2)

            rvm_in = core_state%rtm - core_state%rcm
            call update_xp2_mc_api(pverp, dtime, cloud_frac_inout, &
                                   core_state%rcm, rvm_in, core_state%thlm, &
                                   core_state%wm_zt, core_state%exner, pre_in, pdf_params, &
            rtp2_mc_out, thlp2_mc_out, &
            wprtp_mc_out, wpthlp_mc_out, &
            rtpthlp_mc_out)

            if (clubb_do_deep) then
               dum_core_rknd = 1._core_rknd
            else
               dum_core_rknd = (1._core_rknd - real(cam_in%landfrac(i), kind = core_rknd))
            end if

            ! update turbulent moments based on rain evaporation
            rtp2_in  = rtp2_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * rtp2_mc_out * dtime
            thlp2_in = thlp2_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * thlp2_mc_out * dtime
            if (.not. clubb_do_deep) then
               wprtp_in = wprtp_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * wprtp_mc_out * dtime
               wpthlp_in = wpthlp_in + real(clubb_rnevap_effic, kind = core_rknd) * dum_core_rknd * wpthlp_mc_out * dtime
            endif
!                     rtpthlp_in = rtpthlp_in + rtpthlp_mc_out * dtime

         endif

         !----------------------
         if (do_cldcool) then
            do k=1,pver
               kkhost = pver-k+1
               qrl_clubb(k+1) = real( qrl(i,kkhost)/(cpair*host_mean%pdel(i,kkhost)), kind=core_rknd )
            enddo
            qrl_clubb(1) = qrl_clubb(2)

            rcm_out_zm = zt2zm_api( core_state%rcm )
            qrl_zm = zt2zm_api(qrl_clubb)
            thlp2_rad_out(:) = 0._r8
            call calculate_thlp2_rad_api(pverp, rcm_out_zm, thlprcp_out, qrl_zm, thlp2_rad_out)
            thlp2_in = thlp2_in + thlp2_rad_out * dtime
            thlp2_in = max(thl_tol**2,thlp2_in)
          endif

          !  Check to see if stats should be output, here stats are read into
          !  output arrays to make them conformable to CAM output
          if (l_stats) call stats_end_timestep_clubb(lchnk,i,out_zt,out_zm,&
                                                     out_radzt,out_radzm,out_sfc)

      enddo  ! end time loop
      call t_stopf('adv_clubb_core_ts_loop')
      !**************************************************************************************
      ! End of CLUBB's time-stepping sub-cycles
      !**************************************************************************************

      !=====================================================
      ! Put values back to host model's data structures
      !=====================================================
      if (clubb_do_adv) then
         if (macmic_it .eq. cld_macmic_num_steps) then
            wp2_in=zm2zt_api(wp2_in)
            wpthlp_in=zm2zt_api(wpthlp_in)
            wprtp_in=zm2zt_api(wprtp_in)
            up2_in=zm2zt_api(up2_in)
            vp2_in=zm2zt_api(vp2_in)
            thlp2_in=zm2zt_api(thlp2_in)
            rtp2_in=zm2zt_api(rtp2_in)
            rtpthlp_in=zm2zt_api(rtpthlp_in)

            do k=1,pverp
               thlp2_in(k)=max(thl_tol**2,thlp2_in(k))
               rtp2_in(k)=max(rt_tol**2,rtp2_in(k))
               wp2_in(k)=max(w_tol_sqd,wp2_in(k))
               up2_in(k)=max(w_tol_sqd,up2_in(k))
               vp2_in(k)=max(w_tol_sqd,vp2_in(k))
            enddo
         endif
      endif

      if (linearize_pbl_winds) then
         ! Copy column variables back to pbuf arrays.
         um_pert(i,:) = um_pert_col
         vm_pert(i,:) = vm_pert_col
         upwp_pert(i,:) = upwp_pert_col
         vpwp_pert(i,:) = vpwp_pert_col
         deallocate(um_pert_col)
         deallocate(vm_pert_col)
         deallocate(upwp_pert_col)
         deallocate(vpwp_pert_col)
         deallocate(upwp_sfc_pert)
         deallocate(vpwp_sfc_pert)
         if (abs(cam_in%wsx(i)) < 1.e-12 .and. abs(cam_in%wsy(i)) < 1.e-12) then
            sfc_v_diff_tau(i) = um_pert(i,2) - core_state%um(2)
         else
            sfc_v_diff_tau(i) = ((um_pert(i,2) - core_state%um(2))*cam_in%wsx(i) &
                 + (vm_pert(i,2) - core_state%vm(2))*cam_in%wsy(i)) &
                 / hypot(cam_in%wsx(i), cam_in%wsy(i))
         end if
      end if

      !  Arrays need to be "flipped" to CAM grid
      do k=1,pverp

          host_mean% um(i,k) = real(core_state%um(pverp-k+1), kind = r8)
          host_mean% vm(i,k) = real(core_state%vm(pverp-k+1), kind = r8)

          host_mean% thlm(i,k) = real( core_state%thlm(pverp-k+1), kind=r8 )
          host_mean% rtm (i,k) = real( core_state% rtm(pverp-k+1), kind=r8 )
          host_mean% rcm (i,k) = real( core_state% rcm(pverp-k+1), kind=r8 )

          host_mnts% upwp(i,k)         = real(upwp_in(pverp-k+1), kind = r8)
          host_mnts% vpwp(i,k)         = real(vpwp_in(pverp-k+1), kind = r8)
          host_mnts% wpthvp(i,k)       = real(wpthvp_inout(pverp-k+1), kind = r8)
          host_mnts% wp2thvp(i,k)      = real(wp2thvp_inout(pverp-k+1), kind = r8)
          host_mnts% rtpthvp(i,k)      = real(rtpthvp_inout(pverp-k+1), kind = r8)
          host_mnts% thlpthvp(i,k)     = real(thlpthvp_inout(pverp-k+1), kind = r8)
          host_mnts% up2(i,k)          = real(up2_in(pverp-k+1), kind = r8)
          host_mnts% vp2(i,k)          = real(vp2_in(pverp-k+1), kind = r8)

          host_mnts% wprtp(i,k)        = real(wprtp_in(pverp-k+1), kind = r8)
          host_mnts% wpthlp(i,k)       = real(wpthlp_in(pverp-k+1), kind = r8)
          host_mnts% wp2(i,k)          = real(wp2_in(pverp-k+1), kind = r8)
          host_mnts% wp3(i,k)          = real(wp3_in(pverp-k+1), kind = r8)
          host_mnts% rtp2(i,k)         = real(rtp2_in(pverp-k+1), kind = r8)
          host_mnts% thlp2(i,k)        = real(thlp2_in(pverp-k+1), kind = r8)
          host_mnts% rtpthlp(i,k)      = real(rtpthlp_in(pverp-k+1), kind = r8)

          wprcp(i,k)        = real(wprcp_out(pverp-k+1), kind = r8)
          cloud_frac(i,k)   = min(real(cloud_frac_inout(pverp-k+1), kind = r8),1._r8)
          rcm_in_layer(i,k) = real(rcm_in_layer_out(pverp-k+1), kind = r8)
          cloud_cover(i,k)  = min(real(cloud_cover_out(pverp-k+1), kind = r8),1._r8)

          c2h%khzm(i,k)     = real(khzm_out(pverp-k+1), kind = r8)
          c2h%khzt(i,k)     = real(khzt_out(pverp-k+1), kind = r8)

          qclvar(i,k)       = min(1._r8,real(qclvar_out(pverp-k+1), kind = r8))
          sclrpthvp(i,k,:)  = real(sclrpthvp_inout(pverp-k+1,:), kind = r8)
          pdf_zm_w_1(i,k) = pdf_zm_w_1_inout(pverp-k+1)
          pdf_zm_w_2(i,k) = pdf_zm_w_2_inout(pverp-k+1)
          pdf_zm_varnce_w_1(i,k) = pdf_zm_varnce_w_1_inout(pverp-k+1)
          pdf_zm_varnce_w_2(i,k) = pdf_zm_varnce_w_2_inout(pverp-k+1)
          pdf_zm_mixt_frac(i,k) =  pdf_zm_mixt_frac_inout(pverp-k+1)

          do ixind=1,edsclr_dim
              edsclr_out(k,ixind) = real(edsclr_in(pverp-k+1,ixind), kind = r8)
          enddo

      enddo
      !=================================================
      call clubb_core_state_dealloc( core_state )
