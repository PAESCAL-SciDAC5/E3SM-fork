
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
      enddo


      edsclr_out(i,:,1:edsclr_dim)     = 0._r8
      do k=1,pverp
         edsclr_out(i,k,1:edsclr_dim) = real( core_prog%edsclr(pverp-k+1,1:edsclr_dim), kind = r8 )
      enddo

      !----------------------------------------------
      ! PDF parameters: save to pbuf for BFB restart
      !----------------------------------------------
      do k=1,pverp
          c2h%pdf_zm_w_1       (i,k) = pdf_params_zm%w_1       (pverp-k+1)
          c2h%pdf_zm_w_2       (i,k) = pdf_params_zm%w_2       (pverp-k+1)
          c2h%pdf_zm_varnce_w_1(i,k) = pdf_params_zm%varnce_w_1(pverp-k+1)
          c2h%pdf_zm_varnce_w_2(i,k) = pdf_params_zm%varnce_w_2(pverp-k+1)
          c2h%pdf_zm_mixt_frac (i,k) = pdf_params_zm%mixt_frac (pverp-k+1)
      enddo
