
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
   fcor = 0._core_rknd
