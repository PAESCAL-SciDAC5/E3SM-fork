
   !==================================================================
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

   !  Higher order scalar inputs, set to zero
   wpsclrp_sfc(:)      = 0._core_rknd
   hydromet(:,:)       = 0._core_rknd
   wphydrometp(:,:)    = 0._core_rknd
   wp2hmp(:,:)         = 0._core_rknd
   rtphmp_zt(:,:)      = 0._core_rknd
   thlphmp_zt(:,:)     = 0._core_rknd

   !  Initialize forcings for transported scalars to zero
   sclrm_forcing(:,:)   = 0._core_rknd

