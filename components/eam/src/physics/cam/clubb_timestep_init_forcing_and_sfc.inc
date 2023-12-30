
      !=========================
      ! Non-zero surface fluxes
      !=========================

      !-----------------------------------------
      !  Surface fluxes provided by host model
      !-----------------------------------------
      core_sfc% wpthlp = real(cam_in%shf(i), kind = core_rknd)/(real(cpair, kind = core_rknd)* core_state%rho_ds_zm(1)) ! Sensible heat flux
      core_sfc% wprtp  = real(cam_in%cflx(i,1),    kind = core_rknd)/ core_state%rho_ds_zm(1) ! Latent heat flux
      core_sfc% upwp   = real(cam_in%wsx(i),       kind = core_rknd)/ core_state%rho_ds_zm(1) ! Surface meridional momentum flux
      core_sfc% vpwp   = real(cam_in%wsy(i),       kind = core_rknd)/ core_state%rho_ds_zm(1) ! Surface zonal momentum flux

      !-----------
      ! Apply TMS
      !-----------
       if ( do_tms ) then
         dum_core_rknd = real((ksrftms(i)*host_mean%um(i,pver)), kind = core_rknd)
         core_sfc%upwp = core_sfc%upwp-(dum_core_rknd/core_state%rho_ds_zm(1))
         dum_core_rknd = real((ksrftms(i)*host_mean%vm(i,pver)), kind = core_rknd)
         core_sfc%vpwp = core_sfc%vpwp-(dum_core_rknd/core_state%rho_ds_zm(1))
       endif

      !  Define surface sources for transported variables for diffusion, will
      !  be zero as these tendencies are done in clubb_surface
      !  ---------------
      !  Hui Wan's question, 2023-12:
      !   Is the statement above correct?
      !    - If edsclr contains "eddy passive scalars" such as chemical gases,
      !      then in principle we could have eddy fluxes of such gases.
      !    - If do_expldiff = .t., then thlm and rtm will be transported as "eddy passive tracers",
      !      then shouldn't they have surface fluxes?
      !  ---------------
      core_sfc% wpedsclrp(1:edsclr_dim) = 0._core_rknd

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
             bflx22 = (gravit/real(theta0, kind = r8))*real(core_sfc%wpthlp, kind = r8)
             ustar  = diag_ustar(dum1,bflx22,ubar,zo)
        endif

        !  Compute the surface momentum fluxes, if this is a SCAM simulation
        core_sfc% upwp = -real((host_mean%um(i,pver)*ustar**2/ubar), kind = core_rknd)
        core_sfc% vpwp = -real((host_mean%vm(i,pver)*ustar**2/ubar), kind = core_rknd)

      endif


   !======================================================================
   ! Forcing terms that are zero in the current implementation
   !======================================================================
   !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
   !  forcings already applied through CAM

   core_forcing% um  (1:pverp) = 0._core_rknd
   core_forcing% vm  (1:pverp) = 0._core_rknd

   core_forcing% thlm(1:pverp) = 0._core_rknd
   core_forcing% rtm (1:pverp) = 0._core_rknd

   core_forcing% rtp2   (1:pverp) = 0._core_rknd
   core_forcing% thlp2  (1:pverp) = 0._core_rknd
   core_forcing% rtpthlp(1:pverp) = 0._core_rknd
   core_forcing% wprtp  (1:pverp) = 0._core_rknd
   core_forcing% wpthlp (1:pverp) = 0._core_rknd

   core_forcing% edsclr(:,:) = 0._core_rknd


   !-------------------------------------------------------------------
   !  Determine Coriolis force at given latitude.  This is never used
   !  when CLUBB is implemented in a host model, therefore just set to zero.
   !-------------------------------------------------------------------
   fcoriolis = 0._core_rknd

