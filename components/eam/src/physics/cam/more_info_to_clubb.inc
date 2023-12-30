
      !--------------------------------------------------------------------------

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







      !---------------------------------------------------------
      ! also flip the arrays for the pdf_params_zm variables to have
      ! consistent vertical orientation for all variables in restart file
      ! also need to flip back after calling advance_clubb_core
      do k=1,pverp
         pdf_zm_w_1_inout(k)        = c2h%pdf_zm_w_1       (i,pverp-k+1)
         pdf_zm_w_2_inout(k)        = c2h%pdf_zm_w_2       (i,pverp-k+1)
         pdf_zm_varnce_w_1_inout(k) = c2h%pdf_zm_varnce_w_1(i,pverp-k+1)
         pdf_zm_varnce_w_2_inout(k) = c2h%pdf_zm_varnce_w_2(i,pverp-k+1)
         pdf_zm_mixt_frac_inout(k)  = c2h%pdf_zm_mixt_frac (i,pverp-k+1)
      enddo


      pdf_params    => pdf_params_chnk(i,lchnk)
      pdf_params_zm => pdf_params_zm_chnk(i,lchnk)

      if ( is_first_restart_step() .and. ipdf_call_placement .eq. ipdf_post_advance_fields ) then
         ! assign the values read back from restart file
         ! This is necessary when ipdf_call_placement = 2
         pdf_params_zm%w_1        = pdf_zm_w_1_inout
         pdf_params_zm%w_2        = pdf_zm_w_2_inout
         pdf_params_zm%varnce_w_1 = pdf_zm_varnce_w_1_inout
         pdf_params_zm%varnce_w_2 = pdf_zm_varnce_w_2_inout
         pdf_params_zm%mixt_frac  = pdf_zm_mixt_frac_inout
      end if
      !---------------------------------------------------------




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


