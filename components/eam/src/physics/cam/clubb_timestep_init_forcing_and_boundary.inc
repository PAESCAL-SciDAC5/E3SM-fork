
      !=========================
      ! Non-zero surface fluxes
      !=========================

      !-----------------------------------------
      !  Surface fluxes provided by host model
      !-----------------------------------------
      wpthlp_sfc = real(cam_in%shf(i), kind = core_rknd)/(real(cpair, kind = core_rknd)* core_state%rho_ds_zm(1)) ! Sensible heat flux
      wprtp_sfc  = real(cam_in%cflx(i,1),    kind = core_rknd)/ core_state%rho_ds_zm(1) ! Latent heat flux
      upwp_sfc   = real(cam_in%wsx(i),       kind = core_rknd)/ core_state%rho_ds_zm(1) ! Surface meridional momentum flux
      vpwp_sfc   = real(cam_in%wsy(i),       kind = core_rknd)/ core_state%rho_ds_zm(1) ! Surface zonal momentum flux

      !-----------
      ! Apply TMS
      !-----------
       if ( do_tms ) then
         dum_core_rknd = real((ksrftms(i)*host_mean%um(i,pver)), kind = core_rknd)
         upwp_sfc      = upwp_sfc-(dum_core_rknd/core_state%rho_ds_zm(1))
         dum_core_rknd = real((ksrftms(i)*host_mean%vm(i,pver)), kind = core_rknd)
         vpwp_sfc      = vpwp_sfc-(dum_core_rknd/core_state%rho_ds_zm(1))
       endif

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

             dum1   = real(zt_bot, kind = r8)
             bflx22 = (gravit/real(theta0, kind = r8))*real(wpthlp_sfc, kind = r8)
             ustar  = diag_ustar(dum1,bflx22,ubar,zo)
        endif

        !  Compute the surface momentum fluxes, if this is a SCAM simulation
        upwp_sfc = -real((host_mean%um(i,pver)*ustar**2/ubar), kind = core_rknd)
        vpwp_sfc = -real((host_mean%vm(i,pver)*ustar**2/ubar), kind = core_rknd)

      endif


   !========================================
   ! Fluxes and forcing terms that are zero
   !========================================

   !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
   !  forcings already applied through CAM
   thlm_forcing(1:pverp) = 0._core_rknd
   rtm_forcing(1:pverp)  = 0._core_rknd
   um_forcing(1:pverp)   = 0._core_rknd
   vm_forcing(1:pverp)   = 0._core_rknd

   wprtp_forcing(1:pverp)   = 0._core_rknd
   wpthlp_forcing(1:pverp)  = 0._core_rknd
   rtp2_forcing(1:pverp)    = 0._core_rknd
   thlp2_forcing(1:pverp)   = 0._core_rknd
   rtpthlp_forcing(1:pverp) = 0._core_rknd

   ! rtp3_in and thlp3_in are not currently used in CLUBB's default code.
   rtp3_in(:)  = 0.0_core_rknd
   thlp3_in(:) = 0.0_core_rknd

   !  Higher order scalar inouts, set to zero
   sclrpthvp_inout(1:pverp,:) = 0._core_rknd
             sclrm(1:pverp,:) = 0._core_rknd
           wpsclrp(1:pverp,:) = 0._core_rknd
            sclrp2(1:pverp,:) = 0._core_rknd
          sclrprtp(1:pverp,:) = 0._core_rknd
         sclrpthlp(1:pverp,:) = 0._core_rknd

   !  Define surface sources for transported variables for diffusion, will
   !  be zero as these tendencies are done in clubb_surface
   do ixind=1,edsclr_dim
      wpedsclrp_sfc(ixind) = 0._core_rknd
   enddo

   ice_supersat_frac(1:pverp) = 0._core_rknd

   !  Higher order scalar inputs, set to zero
   wpsclrp_sfc(:)      = 0._core_rknd
   hydromet(:,:)       = 0._core_rknd
   wphydrometp(:,:)    = 0._core_rknd
   wp2hmp(:,:)         = 0._core_rknd
   rtphmp_zt(:,:)      = 0._core_rknd
   thlphmp_zt(:,:)     = 0._core_rknd

   !  Initialize forcings for transported scalars to zero
   sclrm_forcing(:,:)   = 0._core_rknd
   edsclrm_forcing(:,:) = 0._core_rknd

   !  Determine Coriolis force at given latitude.  This is never used
   !  when CLUBB is implemented in a host model, therefore just set
   !  to zero.
   fcoriolis = 0._core_rknd
