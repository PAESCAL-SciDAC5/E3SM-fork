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

