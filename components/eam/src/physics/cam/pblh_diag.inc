
   ! --------------------------------------------------------------------------------- !
   !  DIAGNOSE THE PBL DEPTH                                                           !
   !  this is needed for aerosol code                                                  !
   ! --------------------------------------------------------------------------------- !

   call t_startf('pbl_depth_diag')
   do i=1,ncol
      do k=1,pver
         th(i,k) = state1%t(i,k)*state1%exner(i,k)
         if (use_sgv) then
           thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq) &
                    - state1%q(i,k,ixcldliq))  !PMA corrects thv formula
         else
           thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq))
         end if
      enddo
   enddo

   ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
   do i=1,ncol
      rrho = invrs_gravit*(state1%pdel(i,pver)/dz_g(pver))
      call calc_ustar( state1%t(i,pver), state1%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                       rrho, ustar2(i) )
      call calc_obklen( th(i,pver), thv(i,pver), cam_in%cflx(i,1), cam_in%shf(i), rrho, ustar2(i), &
                        kinheat(i), kinwat(i), kbfs(i), obklen(i) )
   enddo

   dummy2(:) = 0._r8
   dummy3(:) = 0._r8

   where (kbfs .eq. -0.0_r8) kbfs = 0.0_r8

   !  Compute PBL depth according to Holtslag-Boville Scheme
   call pblintd(ncol, thv, state1%zm, state1%u, state1%v, &
                ustar2, obklen, kbfs, pblh, dummy2, &
                state1%zi, cloud_frac(:,1:pver), 1._r8-cam_in%landfrac, dummy3)
   call t_stopf('pbl_depth_diag')

   !  Output the PBL depth
   call outfld('PBLH', pblh, pcols, lchnk)
