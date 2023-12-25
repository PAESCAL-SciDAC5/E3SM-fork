
   ! ------------------------------------------------- !
   ! Diagnose relative cloud water variance            !
   ! ------------------------------------------------- !
   call pbuf_get_field(pbuf, relvar_idx,  relvar)

   if (deep_scheme .eq. 'CLUBB_SGS') then
      relvarmax = 2.0_r8
   else
      relvarmax = 10.0_r8
   endif

   relvar(:,:) = relvarmax  ! default
!
!PMA c20161114: The lower bound of 0.7 is the mean of scattered Cu in Barker et al (1996).
!     With the new formulation the lower bound and is rarely reached.
!
   relvarmin   = 0.7_r8

!PMA c20161114: Xue Zheng identified the issue with small relvar: the original
!               code uses grid mean variance and water content instead of in-cloud
!               quantities.
!               Following equation A7 in Guo et al (2014), relvar is now  calculated
!               using in-cloud variance and in-cloud total water instead of grid
!               mean. This effectively reduces autoconversion rate especially
!               for thin clouds.
!
!

   relvarc(:ncol,:pver)=fillvalue

   if (deep_scheme .ne. 'CLUBB_SGS') then
      if (relvar_fix) then
         where ( host_mean%rcm(:ncol,:pver) > qsmall .and. qclvar(:ncol,:pver) /= 0._r8)  &
              relvar(:ncol,:pver) = min(relvarmax,max(relvarmin, host_mean%rcm(:ncol,:pver)**2/max(qsmall,  &
              cloud_frac(:ncol,:pver)*qclvar(:ncol,:pver)-  &
              (1._r8-cloud_frac(:ncol,:pver))* host_mean%rcm(:ncol,:pver)**2)))
              relvarc(:ncol,:pver) = min(relvarmax,max(relvarmin, host_mean%rcm(:ncol,:pver)**2/max(qsmall,  &
              cloud_frac(:ncol,:pver)*qclvar(:ncol,:pver)-  &
              (1._r8-cloud_frac(:ncol,:pver))* host_mean%rcm(:ncol,:pver)**2)))
      else

         where ( host_mean%rcm(:ncol,:pver) /= 0 .and. qclvar(:ncol,:pver) /= 0) &
              relvar(:ncol,:pver) = min(relvarmax,max(0.001_r8, host_mean%rcm(:ncol,:pver)**2/qclvar(:ncol,:pver)))
      endif
   endif
   call outfld( 'RELVAR',           relvar,                  pcols, lchnk )
   call outfld( 'RELVARC',          relvarc,                 pcols, lchnk )

