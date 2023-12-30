
      zi_out(i,1) = 0._r8

      call clubb_core_fld_alloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc, pverp, edsclr_dim )

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

         core_auxil%p_in_Pa(k+1)   = real( host_mean%pmid(i,kkhost), kind=core_rknd )
         core_auxil%exner(k+1)     = real( 1._r8/host_mean%exner_clubb(i,kkhost), kind=core_rknd )

         thvtmp =  host_mean%t(i,kkhost) &
                  *host_mean%exner_clubb(i,kkhost) &
                  *( 1._r8+zvir*host_mean%q(i,kkhost,ixq)-host_mean%q(i,kkhost,ixcldliq) )
         core_auxil%thv_ds_zt(k+1) = real(thvtmp, kind = core_rknd)                       ! thetav on thermo
      enddo

      core_auxil%p_in_Pa  (1) = core_auxil%p_in_Pa(2)
      core_auxil%exner    (1) = core_auxil%exner  (2)
      core_auxil%thv_ds_zt(1) = core_auxil%thv_ds_zt(2)

      !----------
      do k=1,pver
         kkhost = pver-k+1
         core_auxil%rho_ds_zt(k+1) = real( invrs_gravit*host_mean%pdel(i,kkhost)/dz_g(kkhost), kind=core_rknd )
      enddo
      core_auxil%rho_ds_zt(1) = core_auxil%rho_ds_zt(2)


      !----------
      core_auxil%rho_zt(:)  = core_auxil%rho_ds_zt(:)

      core_auxil%invrs_rho_ds_zt(2:pverp) = 1._core_rknd/core_auxil%rho_ds_zt(2:pverp) ! Inverse ds rho at thermo
      core_auxil%invrs_rho_ds_zt(1)       = core_auxil%invrs_rho_ds_zt(2)

      !  Compute some inputs from the thermodynamic grid to the momentum grid
      !  Note: this can only be done after "call setup_grid_heights_api"

      core_auxil%rho_ds_zm       = zt2zm_api( core_auxil%rho_ds_zt )
      core_auxil%rho_zm          = zt2zm_api( core_auxil%rho_zt    )
      core_auxil%invrs_rho_ds_zm = zt2zm_api( core_auxil%invrs_rho_ds_zt)
      core_auxil%thv_ds_zm       = zt2zm_api( core_auxil%thv_ds_zt )

      !----------
      !  Compute mean w wind on thermo grid, convert from omega to w

      core_auxil%wm_zt(1) = 0._core_rknd
      do k=1,pver
        kkhost = pver-k+1
        dum1 = -1._r8*host_mean%omega(i,kkhost) * real(core_auxil%invrs_rho_ds_zt(k+1),kind=r8) *invrs_gravit
        core_auxil%wm_zt(k+1) = real(dum1,kind=core_rknd)
      enddo

      core_auxil%wm_zm = zt2zm_api( core_auxil%wm_zt )

      !----------

      !  Need to flip arrays around for CLUBB core
      do k=1,pverp

         core_prog%thlm(k) = real(host_mean%thlm(i,pverp-k+1), kind = core_rknd)
         core_prog%rtm (k) = real(host_mean%rtm(i,pverp-k+1), kind = core_rknd)

         core_prog%um(k)      = real(host_mean%um(i,pverp-k+1), kind = core_rknd)
         core_prog%vm(k)      = real(host_mean%vm(i,pverp-k+1), kind = core_rknd)

         core_prog%up2(k)     = real(host_mnts%up2(i,pverp-k+1), kind = core_rknd)
         core_prog%vp2(k)     = real(host_mnts%vp2(i,pverp-k+1), kind = core_rknd)
         core_prog%wp2(k)     = real(host_mnts%wp2(i,pverp-k+1), kind = core_rknd)
         core_prog%wp3(k)     = real(host_mnts%wp3(i,pverp-k+1), kind = core_rknd)

         core_prog%upwp(k)    = real(host_mnts%upwp(i,pverp-k+1), kind = core_rknd)
         core_prog%vpwp(k)    = real(host_mnts%vpwp(i,pverp-k+1), kind = core_rknd)

         core_prog%rtp2   (k) = real(host_mnts%rtp2(i,pverp-k+1), kind = core_rknd)
         core_prog%thlp2  (k) = real(host_mnts%thlp2(i,pverp-k+1), kind = core_rknd)
         core_prog%wprtp  (k) = real(host_mnts%wprtp(i,pverp-k+1), kind = core_rknd)
         core_prog%wpthlp (k) = real(host_mnts%wpthlp(i,pverp-k+1), kind = core_rknd)
         core_prog%rtpthlp(k) = real(host_mnts%rtpthlp(i,pverp-k+1), kind = core_rknd)

         core_diag%wpthvp  (k) = real(host_mnts%wpthvp  (i,pverp-k+1), kind = core_rknd)
         core_diag%wp2thvp (k) = real(host_mnts%wp2thvp (i,pverp-k+1), kind = core_rknd)
         core_diag%rtpthvp (k) = real(host_mnts%rtpthvp (i,pverp-k+1), kind = core_rknd)
         core_diag%thlpthvp(k) = real(host_mnts%thlpthvp(i,pverp-k+1), kind = core_rknd)

         core_diag%rcm       (k) = real(host_mean%rcm(i,pverp-k+1), kind = core_rknd)
         core_diag%cloud_frac(k) = real( cloud_frac(i,pverp-k+1), kind = core_rknd)
      enddo

         ! also flip the arrays for the pdf_params_zm variables to have
         ! consistent vertical orientation for all variables in restart file
         ! also need to flip back after calling advance_clubb_core

      do k=1,pverp
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
            um_pert(i,:) = core_prog%um
            vm_pert(i,:) = core_prog%vm
            upwp_pert(i,:) = core_prog%upwp
            vpwp_pert(i,:) = core_prog%vpwp
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
            upwp_sfc_pert = core_sfc%upwp + pert_tau / core_auxil%rho_ds_zm(1)
            vpwp_sfc_pert = core_sfc%vpwp
         else
            upwp_sfc_pert = core_sfc%upwp + cam_in%wsx(i) * &
                 (pert_tau / (core_auxil%rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
            vpwp_sfc_pert = core_sfc%vpwp + cam_in%wsy(i) * &
                 (pert_tau / (core_auxil%rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
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

          core_prog%up2 = zt2zm_api( core_prog%up2 )
          core_prog%vp2 = zt2zm_api( core_prog%vp2 )
          core_prog%wp2 = zt2zm_api( core_prog%wp2 )

          core_prog%wpthlp  = zt2zm_api( core_prog%wpthlp  )
          core_prog%wprtp   = zt2zm_api( core_prog%wprtp   )
          core_prog%thlp2   = zt2zm_api( core_prog%thlp2   )
          core_prog%rtp2    = zt2zm_api( core_prog%rtp2    )
          core_prog%rtpthlp = zt2zm_api( core_prog%rtpthlp )

          do k=1,pverp
            core_prog%thlp2(k)=max(thl_tol**2,core_prog%thlp2(k))
            core_prog%rtp2 (k)=max( rt_tol**2,core_prog%rtp2 (k))

            core_prog%up2(k)=max(w_tol_sqd,core_prog%up2(k))
            core_prog%vp2(k)=max(w_tol_sqd,core_prog%vp2(k))
            core_prog%wp2(k)=max(w_tol_sqd,core_prog%wp2(k))
          enddo
        endif
      endif

      !-----------------------------------------------
      ! Get input values for passive tracers
      !-----------------------------------------------
      core_prog%edsclr(:,:) = 0._core_rknd
 
      icnt=0
      do ixind=1,pcnst
         if (lq(ixind))  then
            icnt=icnt+1
            do k=1,pver
               core_prog%edsclr(k+1,icnt) = real(host_mean%q(i,pver-k+1,ixind), kind = core_rknd)
            enddo
            core_prog%edsclr(1,icnt) = core_prog%edsclr(2,icnt)
         end if
      enddo

      if (do_expldiff) then
        do k=1,pver
          core_prog%edsclr(k+1,icnt+1) = real(host_mean%thlm(i,pver-k+1), kind = core_rknd)
          core_prog%edsclr(k+1,icnt+2) = real(host_mean%rtm (i,pver-k+1), kind = core_rknd)
        enddo

        core_prog%edsclr(1,icnt+1) = core_prog%edsclr(2,icnt+1)
        core_prog%edsclr(1,icnt+2) = core_prog%edsclr(2,icnt+2)
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

      !----------
      do k=1,pver
         kkhost = pver-k+1
         rfrzm(k+1)     = real( host_mean%q(i,kkhost,ixcldice),                 kind=core_rknd )
         radf(k+1)      = real( radf_clubb(i,kkhost),                           kind=core_rknd )
      enddo
      radf(1)      = radf(2)
      rfrzm(1)     = rfrzm(2)

      ice_supersat_frac(1:pverp) = 0._core_rknd
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

#include "clubb_timestep_init_forcing_and_sfc.inc"

#include "clubb_terms_not_in_eqns.inc"

         !  Advance CLUBB CORE one timestep in the future
         call t_startf('advance_clubb_core')
         !Balli- to do: check whether initent-ins and intent-inouts are actually what they say

         call advance_clubb_core_api &
              ( l_implemented, dtime, fcoriolis, sfc_elevation, hydromet_dim, & ! intent(in)
              !
              core_forcing%thlm, &
              core_forcing%rtm,  &
              core_forcing%um,   &
              core_forcing%vm,   &! intent(in)
              sclrm_forcing,     &
              core_forcing%edsclr,  &
              core_forcing%wprtp,   &! intent(in)
              core_forcing%wpthlp,  &
              core_forcing%rtp2,    &
              core_forcing%thlp2,   &! intent(in)
              core_forcing%rtpthlp, &
              !
              core_auxil%wm_zm, &
              core_auxil%wm_zt, &! intent(in)
              !
              core_sfc%wpthlp,    &
              core_sfc%wprtp,     &
              core_sfc%upwp,      &
              core_sfc%vpwp,      &! intent(in)
              wpsclrp_sfc,        &
              core_sfc%wpedsclrp, &! intent(in)
              !
              core_auxil%p_in_Pa,         &
              core_auxil%rho_zm,          &
              core_auxil%rho_zt,          &
              core_auxil%exner,           &! intent(in)
              core_auxil%rho_ds_zm,       &
              core_auxil%rho_ds_zt,       &
              core_auxil%invrs_rho_ds_zm, &! intent(in)
              core_auxil%invrs_rho_ds_zt, &
              core_auxil%thv_ds_zm,       &
              core_auxil%thv_ds_zt,       &
              !
              hydromet,    &! intent(in)
              rfrzm, radf, &! intent(in)
#ifdef CLUBBND_CAM
              varmu2,      &! intent(in)
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
              sclrpthvp_inout, &! intent(inout)
              pdf_params, pdf_params_zm, &! intent(inout)
              !
              core_diag%khzm,    &
              core_diag%khzt,    &
              core_diag%qclvar,  &
              core_diag%thlprcp, &! intent(out)
              core_diag%wprcp,   &
              !
              ice_supersat_frac,                 &! intent(out)
              rcm_in_layer_out, cloud_cover_out, &! intent(out)
              upwp_sfc_pert, vpwp_sfc_pert,      &! intent(in)
              um_pert_col, vm_pert_col,          &! 
              upwp_pert_col, vpwp_pert_col)       ! intent(inout)
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
            do k=1,pver
               kkhost = pver-k+1
               qrl_clubb(k+1) = real( qrl(i,kkhost)/(cpair*host_mean%pdel(i,kkhost)), kind=core_rknd )
            enddo
            qrl_clubb(1) = qrl_clubb(2)

            rcm_out_zm = zt2zm_api( core_diag%rcm )
            qrl_zm = zt2zm_api(qrl_clubb)
            thlp2_rad_out(:) = 0._r8
            call calculate_thlp2_rad_api(pverp, rcm_out_zm, core_diag%thlprcp, qrl_zm, thlp2_rad_out)
            core_prog%thlp2 = core_prog%thlp2 + thlp2_rad_out * dtime
            core_prog%thlp2 = max(thl_tol**2,core_prog%thlp2)
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
            core_prog%wp2     = zm2zt_api( core_prog%wp2     )
            core_prog%wpthlp  = zm2zt_api( core_prog%wpthlp  )
            core_prog%wprtp   = zm2zt_api( core_prog%wprtp   )
            core_prog%up2     = zm2zt_api( core_prog%up2     )
            core_prog%vp2     = zm2zt_api( core_prog%vp2     )
            core_prog%thlp2   = zm2zt_api( core_prog%thlp2   )
            core_prog%rtp2    = zm2zt_api( core_prog%rtp2    )
            core_prog%rtpthlp = zm2zt_api( core_prog%rtpthlp )

            do k=1,pverp
               core_prog%thlp2(k) = max( thl_tol**2,core_prog%thlp2(k) )
               core_prog%rtp2 (k) = max( rt_tol**2, core_prog%rtp2 (k) )
               core_prog%wp2  (k) = max( w_tol_sqd, core_prog%wp2  (k) )
               core_prog%up2  (k) = max( w_tol_sqd, core_prog%up2  (k) )
               core_prog%vp2  (k) = max( w_tol_sqd, core_prog%vp2  (k) )
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
            sfc_v_diff_tau(i) = um_pert(i,2) - core_prog%um(2)
         else
            sfc_v_diff_tau(i) = ((um_pert(i,2) - core_prog%um(2))*cam_in%wsx(i) &
                 + (vm_pert(i,2) - core_prog%vm(2))*cam_in%wsy(i)) &
                 / hypot(cam_in%wsx(i), cam_in%wsy(i))
         end if
      end if

      !  Arrays need to be "flipped" to CAM grid
      do k=1,pverp

          host_mean% um(i,k) = real(core_prog%um(pverp-k+1), kind = r8)
          host_mean% vm(i,k) = real(core_prog%vm(pverp-k+1), kind = r8)

          host_mean% thlm(i,k) = real( core_prog%thlm(pverp-k+1), kind=r8 )
          host_mean% rtm (i,k) = real( core_prog% rtm(pverp-k+1), kind=r8 )
          host_mean% rcm (i,k) = real( core_diag% rcm(pverp-k+1), kind=r8 )

          host_mnts% upwp(i,k)         = real(core_prog%upwp(pverp-k+1), kind = r8)
          host_mnts% vpwp(i,k)         = real(core_prog%vpwp(pverp-k+1), kind = r8)

          host_mnts% wpthvp(i,k)       = real( core_diag%wpthvp  (pverp-k+1), kind = r8)
          host_mnts% wp2thvp(i,k)      = real( core_diag%wp2thvp (pverp-k+1), kind = r8)
          host_mnts% rtpthvp(i,k)      = real( core_diag%rtpthvp (pverp-k+1), kind = r8)
          host_mnts% thlpthvp(i,k)     = real( core_diag%thlpthvp(pverp-k+1), kind = r8)

          host_mnts% up2(i,k)          = real(core_prog%up2(pverp-k+1), kind = r8)
          host_mnts% vp2(i,k)          = real(core_prog%vp2(pverp-k+1), kind = r8)
          host_mnts% wp2(i,k)          = real(core_prog%wp2(pverp-k+1), kind = r8)
          host_mnts% wp3(i,k)          = real(core_prog%wp3(pverp-k+1), kind = r8)

          host_mnts% rtp2(i,k)         = real( core_prog%rtp2   (pverp-k+1), kind = r8)
          host_mnts% thlp2(i,k)        = real( core_prog%thlp2  (pverp-k+1), kind = r8)
          host_mnts% wprtp(i,k)        = real( core_prog%wprtp  (pverp-k+1), kind = r8)
          host_mnts% wpthlp(i,k)       = real( core_prog%wpthlp (pverp-k+1), kind = r8)
          host_mnts% rtpthlp(i,k)      = real( core_prog%rtpthlp(pverp-k+1), kind = r8)

          wprcp(i,k)        = real( core_diag%wprcp(pverp-k+1), kind = r8)
          cloud_frac(i,k)   = min(real(core_diag%cloud_frac(pverp-k+1), kind = r8),1._r8)
          rcm_in_layer(i,k) = real(rcm_in_layer_out(pverp-k+1), kind = r8)
          cloud_cover(i,k)  = min(real(cloud_cover_out(pverp-k+1), kind = r8),1._r8)

          c2h%khzm(i,k)     = real( core_diag%khzm(pverp-k+1), kind=r8 )
          c2h%khzt(i,k)     = real( core_diag%khzt(pverp-k+1), kind=r8 )

          qclvar(i,k)   = real( min( 1._core_rknd, core_diag%qclvar(pverp-k+1) ), kind=r8 )

          sclrpthvp(i,k,:)  = real(sclrpthvp_inout(pverp-k+1,:), kind = r8)
          pdf_zm_w_1(i,k) = pdf_zm_w_1_inout(pverp-k+1)
          pdf_zm_w_2(i,k) = pdf_zm_w_2_inout(pverp-k+1)
          pdf_zm_varnce_w_1(i,k) = pdf_zm_varnce_w_1_inout(pverp-k+1)
          pdf_zm_varnce_w_2(i,k) = pdf_zm_varnce_w_2_inout(pverp-k+1)
          pdf_zm_mixt_frac(i,k) =  pdf_zm_mixt_frac_inout(pverp-k+1)

      enddo

      edsclr_out(i,:,1:edsclr_dim)     = 0._r8
      do k=1,pverp
         edsclr_out(i,k,1:edsclr_dim) = real( core_prog%edsclr(pverp-k+1,1:edsclr_dim), kind = r8 )
      enddo

      !=================================================
      call clubb_core_fld_dealloc( core_auxil, core_prog, core_diag, core_forcing, core_sfc )
