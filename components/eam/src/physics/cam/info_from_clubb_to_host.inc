
      !=====================================================
      ! Put values back to host model's data structures
      !=====================================================

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

