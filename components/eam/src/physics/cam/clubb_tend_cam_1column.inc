      invrs_gravit = 1._r8 / gravit
 
      !------------------------------------------------------------------------------
      !  Set time_elapsed to host model time step, this is for CLUBB's budget stats
      !------------------------------------------------------------------------------
      time_elapsed = hdtime_core_rknd

      !------------------------------------------------------------------------------
      ! CLUBB's z-grid
      !------------------------------------------------------------------------------
      !  Define the CLUBB momentum grid (in height, units of m)
      do k=1,pverp
         dum1 = state1%zi(i,pverp-k+1)-state1%zi(i,pver+1)
         zi_g(k) = real(dum1, kind = core_rknd)
      enddo

      !  Define the CLUBB thermodynamic grid (in units of m)
      do k=1,pver
         dum1 = state1%zm(i,pver-k+1)-state1%zi(i,pver+1)
         zt_g(k+1) = real(dum1, kind = core_rknd)
         dz_g(k) = state1%zi(i,k)-state1%zi(i,k+1)  ! compute thickness
      enddo

      !  Thermodynamic ghost point is below surface
      zt_g(1) = -1._core_rknd*zt_g(2)

      !  Set the elevation of the surface
      sfc_elevation = real(state1%zi(i,pver+1), kind = core_rknd)

      !------------------------------------------------------------------------------
      !  Compute virtual potential temperature, which is needed for CLUBB

      call cnst_get_ind('Q',ixq)
      call cnst_get_ind('CLDLIQ',ixcldliq)
      call cnst_get_ind('CLDICE',ixcldice)

      do k=1,pver
         thv1d(k) = state1%t(i,k)*exner_clubb(i,k)*(1._r8+zvir*state1%q(i,k,ixq)&
                   -state1%q(i,k,ixcldliq))
      enddo

      !  Compute thermodynamic stuff needed for CLUBB on thermo levels.
      !  Inputs for the momentum levels are set below setup_clubb core
      do k=1,pver
         p_in_Pa(k+1)         = real(state1%pmid(i,pver-k+1), kind = core_rknd)             ! Pressure profile
         exner(k+1)           = 1._core_rknd/real(exner_clubb(i,pver-k+1), kind = core_rknd)
         rho(i,k+1)           = invrs_gravit*state1%pdel(i,pver-k+1)/dz_g(pver-k+1)
         rho_ds_zt(k+1)       = real(rho(i,k+1), kind = core_rknd)
         invrs_rho_ds_zt(k+1) = 1._core_rknd/(rho_ds_zt(k+1))                               ! Inverse ds rho at thermo

         thv_ds_zt(k+1)       = real(thv1d(pver-k+1), kind = core_rknd)                       ! thetav on thermo

         rfrzm(k+1)           = real(state1%q(i,pver-k+1,ixcldice), kind = core_rknd)
         radf(k+1)            = real(radf_clubb(i,pver-k+1), kind = core_rknd)
         dum1                 = qrl(i,pver-k+1)/(cpair*state1%pdel(i,pver-k+1))
         qrl_clubb(k+1)       = real(dum1, kind = core_rknd)
      enddo

      !  Below computes the same stuff for the ghost point.  May or may
      !  not be needed, just to be safe to avoid NaN's
      rho_ds_zt(1)       = rho_ds_zt(2)
      invrs_rho_ds_zt(1) = invrs_rho_ds_zt(2)
      rho(i,1)           = rho(i,2)     !rho_ds_zt(2)
      thv_ds_zt(1)       = thv_ds_zt(2)
      rho_zt(:)          = rho_ds_zt(:) !rho(i,:)
      p_in_Pa(1)         = p_in_Pa(2)
      exner(1)           = exner(2)
      rfrzm(1)           = rfrzm(2)
      radf(1)            = radf(2)
      qrl_clubb(1)       = qrl_clubb(2)

      !  Compute mean w wind on thermo grid, convert from omega to w
      wm_zt(1) = 0._core_rknd
      do k=1,pver
        dum1 = -1._r8*state1%omega(i,pver-k+1)*real(invrs_rho_ds_zt(k+1), kind = r8)*invrs_gravit
        wm_zt(k+1) = real(dum1, kind = core_rknd)
      enddo

      ! ------------------------------------------------- !
      ! Begin case specific code for SCAM cases.          !
      ! This section of code block NOT called in          !
      ! global simulations                                !
      ! ------------------------------------------------- !

      if (single_column) then

        !  Initialize zo if variable ustar is used

        if (cam_in%landfrac(i) .ge. 0.5_r8) then
           zo = 0.035_r8
        else
           zo = 0.0001_r8
        endif

        !  Compute surface wind (ubar)
        ubar = sqrt(host_mean%um(i,pver)**2+host_mean%vm(i,pver)**2)
        if (ubar .lt. 0.25_r8) ubar = 0.25_r8

        !  Below denotes case specifics for surface momentum
        !  and thermodynamic fluxes, depending on the case

        !  Define ustar (based on case, if not variable)
        ustar = 0.25_r8   ! Initialize ustar in case no case

        if(trim(scm_clubb_iop_name) .eq. 'BOMEX_5day') then
           ustar = 0.28_r8
        endif

        if(trim(scm_clubb_iop_name) .eq. 'ATEX_48hr') then
           ustar = 0.30_r8
        endif

        if(trim(scm_clubb_iop_name) .eq. 'RICO_3day') then
           ustar = 0.28_r8
        endif

        if(trim(scm_clubb_iop_name) .eq. 'arm97' .or. trim(scm_clubb_iop_name) .eq. 'gate' .or. &
           trim(scm_clubb_iop_name) .eq. 'toga' .or. trim(scm_clubb_iop_name) .eq. 'mpace' .or. &
           trim(scm_clubb_iop_name) .eq. 'ARM_CC') then

             dum1   = real(zt_g(2), kind = r8)
             bflx22 = (gravit/real(theta0, kind = r8))*real(wpthlp_sfc, kind = r8)
             ustar  = diag_ustar(dum1,bflx22,ubar,zo)
        endif

        !  Compute the surface momentum fluxes, if this is a SCAM simulation
        upwp_sfc = -real((host_mean%um(i,pver)*ustar**2/ubar), kind = core_rknd)
        vpwp_sfc = -real((host_mean%vm(i,pver)*ustar**2/ubar), kind = core_rknd)

      endif

      !  Set stats output and increment equal to CLUBB and host dt
      stats_tsamp = dtime
      stats_tout  = hdtime_core_rknd

      !  Heights need to be set at each timestep.  Therefore, recall
      !  setup_grid and setup_parameters for this.

      !  Read in parameters for CLUBB.  Pack the default and updated (via nml)
      !  tunable parameters into clubb_params
      call read_parameters_api( -99, "", clubb_params )

      !  Set-up CLUBB core at each CLUBB call because heights can change
      call setup_grid_heights_api(l_implemented, grid_type, zi_g(2), &
         zi_g(1), zi_g, zt_g)

      call setup_parameters_api(zi_g(2), clubb_params, pverp, grid_type, &
        zi_g, zt_g, err_code)

      !  Compute some inputs from the thermodynamic grid
      !  to the momentum grid
      rho_ds_zm       = zt2zm_api(rho_ds_zt)
      rho_zm          = zt2zm_api(rho_zt)
      invrs_rho_ds_zm = zt2zm_api(invrs_rho_ds_zt)
      thv_ds_zm       = zt2zm_api(thv_ds_zt)
      wm_zm           = zt2zm_api(wm_zt)

      !  Surface fluxes provided by host model
      wpthlp_sfc = real(cam_in%shf(i), kind = core_rknd)/(real(cpair, kind = core_rknd)*rho_ds_zm(1)) ! Sensible heat flux
      wprtp_sfc  = real(cam_in%cflx(i,1), kind = core_rknd)/rho_ds_zm(1)                              ! Latent heat flux
      upwp_sfc   = real(cam_in%wsx(i), kind = core_rknd)/rho_ds_zm(1)                                 ! Surface meridional momentum flux
      vpwp_sfc   = real(cam_in%wsy(i), kind = core_rknd)/rho_ds_zm(1)                                 ! Surface zonal momentum flux

      ! ------------------------------------------------- !
      ! Apply TMS                                         !
      ! ------------------------------------------------- !
       if ( do_tms ) then
         dum_core_rknd = real((ksrftms(i)*state1%u(i,pver)), kind = core_rknd)
         upwp_sfc      = upwp_sfc-(dum_core_rknd/rho_ds_zm(1))
         dum_core_rknd = real((ksrftms(i)*state1%v(i,pver)), kind = core_rknd)
         vpwp_sfc      = vpwp_sfc-(dum_core_rknd/rho_ds_zm(1))
       endif

      !  Need to flip arrays around for CLUBB core
      do k=1,pverp
         um_in(k)      = real(host_mean%um(i,pverp-k+1), kind = core_rknd)
         vm_in(k)      = real(host_mean%vm(i,pverp-k+1), kind = core_rknd)
         upwp_in(k)    = real(upwp(i,pverp-k+1), kind = core_rknd)
         vpwp_in(k)    = real(vpwp(i,pverp-k+1), kind = core_rknd)
         up2_in(k)     = real(up2(i,pverp-k+1), kind = core_rknd)
         vp2_in(k)     = real(vp2(i,pverp-k+1), kind = core_rknd)
         wp2_in(k)     = real(wp2(i,pverp-k+1), kind = core_rknd)
         wp3_in(k)     = real(wp3(i,pverp-k+1), kind = core_rknd)
         rtp2_in(k)    = real(rtp2(i,pverp-k+1), kind = core_rknd)
         thlp2_in(k)   = real(thlp2(i,pverp-k+1), kind = core_rknd)
         thlm_in(k)    = real(host_mean%thlm(i,pverp-k+1), kind = core_rknd)
         rtm_in(k)     = real(host_mean%rtm(i,pverp-k+1), kind = core_rknd)
         rvm_in(k)     = real(rvm(i,pverp-k+1), kind = core_rknd)
         wprtp_in(k)   = real(wprtp(i,pverp-k+1), kind = core_rknd)
         wpthlp_in(k)  = real(wpthlp(i,pverp-k+1), kind = core_rknd)
         rtpthlp_in(k) = real(rtpthlp(i,pverp-k+1), kind = core_rknd)

         wpthvp_inout(k)     = real(wpthvp(i,pverp-k+1), kind = core_rknd)
         wp2thvp_inout(k)    = real(wp2thvp(i,pverp-k+1), kind = core_rknd)
         rtpthvp_inout(k)    = real(rtpthvp(i,pverp-k+1), kind = core_rknd)
         thlpthvp_inout(k)   = real(thlpthvp(i,pverp-k+1), kind = core_rknd)
         rcm_inout(k)        = real(host_mean%rcm(i,pverp-k+1), kind = core_rknd)
         cloud_frac_inout(k) = real(cloud_frac(i,pverp-k+1), kind = core_rknd)

         ! also flip the arrays for the pdf_params_zm variables to have
         ! consistent vertical orientation for all variables in restart file
         ! also need to flip back after calling advance_clubb_core

         pdf_zm_w_1_inout(k) = pdf_zm_w_1(i,pverp-k+1)
         pdf_zm_w_2_inout(k) = pdf_zm_w_2(i,pverp-k+1)
         pdf_zm_varnce_w_1_inout(k) = pdf_zm_varnce_w_1(i,pverp-k+1)
         pdf_zm_varnce_w_2_inout(k) = pdf_zm_varnce_w_2(i,pverp-k+1)
         pdf_zm_mixt_frac_inout(k) =  pdf_zm_mixt_frac(i,pverp-k+1)

         !  Higher order scalar inouts, set to zero
         sclrpthvp_inout(k,:)= 0._core_rknd
         sclrm(k,:)          = 0._core_rknd
         wpsclrp(k,:)        = 0._core_rknd
         sclrp2(k,:)         = 0._core_rknd
         sclrprtp(k,:)       = 0._core_rknd
         sclrpthlp(k,:)      = 0._core_rknd

         if (k .ne. 1) then
            pre_in(k)    = real(prer_evap(i,pverp-k+1), kind = core_rknd)
         endif

      enddo

      if (linearize_pbl_winds) then
         ! Each host model time step, reset the perturbed variables to be equal to
         ! the unperturbed values.
         if (macmic_it == 1) then
            um_pert(i,:) = um_in
            vm_pert(i,:) = vm_in
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
            upwp_sfc_pert = upwp_sfc + pert_tau / rho_ds_zm(1)
            vpwp_sfc_pert = vpwp_sfc
         else
            upwp_sfc_pert = upwp_sfc + cam_in%wsx(i) * &
                 (pert_tau / (rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
            vpwp_sfc_pert = vpwp_sfc + cam_in%wsy(i) * &
                 (pert_tau / (rho_ds_zm(1) * hypot(cam_in%wsx(i), cam_in%wsy(i))))
         end if
      else
         nullify(upwp_sfc_pert)
         nullify(vpwp_sfc_pert)
         nullify(um_pert_col)
         nullify(vm_pert_col)
         nullify(upwp_pert_col)
         nullify(vpwp_pert_col)
      end if

      pre_in(1) = pre_in(2)

      !  Initialize these to prevent crashing behavior
      edsclr_in(:,:)      = 0._core_rknd
      edsclr_out(:,:)     = 0._r8


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

      !  Do the same for tracers
      icnt=0
      do ixind=1,pcnst
         if (lq(ixind))  then
            icnt=icnt+1
            do k=1,pver
               edsclr_in(k+1,icnt) = real(state1%q(i,pver-k+1,ixind), kind = core_rknd)
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

      rho_in(:) = real(rho(i,:), kind = core_rknd)

      ! --------------------------------------------------------- !
      ! Compute cloud-top radiative cooling contribution to CLUBB !
      ! --------------------------------------------------------- !

      ! Sandbox version of code to take into account meso organization

      if (clubb_do_deep) then
         orgparam = 0._r8
         delpavg = 0._r8

         do k = 1, pver
           if (abs(prer_evap(i,k)) .gt. 0._r8) then
             orgparam = orgparam + (abs(prer_evap(i,k)) * 1000._r8 * 1000._r8 * 2._r8 ) * state1%pdel(i,k)
             delpavg = delpavg + state1%pdel(i,k)
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

      call t_startf('adv_clubb_core_ts_loop')
      do t=1,n_clubb_core_step    ! do needed number of "sub" timesteps for each CAM step

         !  Increment the statistics then being stats timestep
         if (l_stats) then
            time_elapsed = time_elapsed+dtime
            call stats_begin_timestep_api(time_elapsed, 1, 1)
         endif

         !  Advance CLUBB CORE one timestep in the future
         call t_startf('advance_clubb_core')
         !Balli- to do: check whether initent-ins and intent-inouts are actually what they say

         call advance_clubb_core_api &
              ( l_implemented, dtime, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
              thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &         ! intent(in)
              sclrm_forcing, edsclrm_forcing, wprtp_forcing, &             ! intent(in)
              wpthlp_forcing, rtp2_forcing, thlp2_forcing, &               ! intent(in)
              rtpthlp_forcing, wm_zm, wm_zt, &                             ! intent(in)
              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &                 ! intent(in)
              wpsclrp_sfc, wpedsclrp_sfc, &                                ! intent(in)
              p_in_Pa, rho_zm, rho_in, exner, &                            ! intent(in)
              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                     ! intent(in)
              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &           ! intent(in)
              rfrzm, radf, &                                               ! intent(in)
#ifdef CLUBBND_CAM
              varmu2, &                                                    ! intent(in)
#endif
              wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &                ! intent(in)
              host_dx, host_dy, &                                          ! intent(in)
              um_in, vm_in, upwp_in, &                                     ! intent(inout)
              vpwp_in, up2_in, vp2_in, &                                   ! intent(inout)
              thlm_in, rtm_in, wprtp_in, wpthlp_in, &                      ! intent(inout)
              wp2_in, wp3_in, rtp2_in, &                                   ! intent(inout)
              rtp3_in, thlp2_in, thlp3_in, rtpthlp_in, &                   ! intent(inout)
              sclrm,   &                                                   ! intent(inout)
              sclrp2, sclrprtp, sclrpthlp, &                               ! intent(inout)
              wpsclrp, edsclr_in, err_code, &                              ! intent(inout)
              rcm_inout, cloud_frac_inout, &                               ! intent(inout)
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
                 "LAT (radians): ", state1%lat(i), &
                 "LON (radians): ", state1%lon(i), &
                 "LAT (degrees): ", rad_to_deg*state1%lat(i), &
                 "LON (degrees): ", rad_to_deg*state1%lon(i), &
                 "Global Column Number: ", get_gcol_p(lchnk,i)
            call endrun('clubb_tend_cam:  Fatal error in CLUBB library'//errmsg(__FILE__,__LINE__))
         end if

         pdf_zm_w_1_inout = pdf_params_zm%w_1
         pdf_zm_w_2_inout = pdf_params_zm%w_2
         pdf_zm_varnce_w_1_inout = pdf_params_zm%varnce_w_1
         pdf_zm_varnce_w_2_inout = pdf_params_zm%varnce_w_2
         pdf_zm_mixt_frac_inout = pdf_params_zm%mixt_frac

         if (do_rainturb) then
            rvm_in = rtm_in - rcm_inout
            call update_xp2_mc_api(pverp, dtime, cloud_frac_inout, &
            rcm_inout, rvm_in, thlm_in, wm_zt, exner, pre_in, pdf_params, &
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

         if (do_cldcool) then

            rcm_out_zm = zt2zm_api(rcm_inout)
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
            sfc_v_diff_tau(i) = um_pert(i,2) - um_in(2)
         else
            sfc_v_diff_tau(i) = ((um_pert(i,2) - um_in(2))*cam_in%wsx(i) &
                 + (vm_pert(i,2) - vm_in(2))*cam_in%wsy(i)) &
                 / hypot(cam_in%wsx(i), cam_in%wsy(i))
         end if
      end if

      !  Arrays need to be "flipped" to CAM grid
      do k=1,pverp

          host_mean%um(i,k) = real(um_in(pverp-k+1), kind = r8)
          host_mean%vm(i,k) = real(vm_in(pverp-k+1), kind = r8)
          upwp(i,k)         = real(upwp_in(pverp-k+1), kind = r8)
          vpwp(i,k)         = real(vpwp_in(pverp-k+1), kind = r8)
          wpthvp(i,k)       = real(wpthvp_inout(pverp-k+1), kind = r8)
          wp2thvp(i,k)      = real(wp2thvp_inout(pverp-k+1), kind = r8)
          rtpthvp(i,k)      = real(rtpthvp_inout(pverp-k+1), kind = r8)
          thlpthvp(i,k)     = real(thlpthvp_inout(pverp-k+1), kind = r8)
          up2(i,k)          = real(up2_in(pverp-k+1), kind = r8)
          vp2(i,k)          = real(vp2_in(pverp-k+1), kind = r8)
          host_mean%thlm(i,k) = real(thlm_in(pverp-k+1), kind = r8)
          host_mean%rtm(i,k)  = real(rtm_in(pverp-k+1), kind = r8)
          wprtp(i,k)        = real(wprtp_in(pverp-k+1), kind = r8)
          wpthlp(i,k)       = real(wpthlp_in(pverp-k+1), kind = r8)
          wp2(i,k)          = real(wp2_in(pverp-k+1), kind = r8)
          wp3(i,k)          = real(wp3_in(pverp-k+1), kind = r8)
          rtp2(i,k)         = real(rtp2_in(pverp-k+1), kind = r8)
          thlp2(i,k)        = real(thlp2_in(pverp-k+1), kind = r8)
          rtpthlp(i,k)      = real(rtpthlp_in(pverp-k+1), kind = r8)
          host_mean%rcm(i,k) = real(rcm_inout(pverp-k+1), kind = r8)
          wprcp(i,k)        = real(wprcp_out(pverp-k+1), kind = r8)
          cloud_frac(i,k)   = min(real(cloud_frac_inout(pverp-k+1), kind = r8),1._r8)
          rcm_in_layer(i,k) = real(rcm_in_layer_out(pverp-k+1), kind = r8)
          cloud_cover(i,k)  = min(real(cloud_cover_out(pverp-k+1), kind = r8),1._r8)
          zt_out(i,k)       = real(zt_g(pverp-k+1), kind = r8)
          zi_out(i,k)       = real(zi_g(pverp-k+1), kind = r8)
          khzm(i,k)         = real(khzm_out(pverp-k+1), kind = r8)
          khzt(i,k)         = real(khzt_out(pverp-k+1), kind = r8)
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

      !  Fill up arrays needed for McICA.  Note we do not want the ghost point,
      !   thus why the second loop is needed.

      zi_out(i,1) = 0._r8

