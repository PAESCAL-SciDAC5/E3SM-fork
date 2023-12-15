   !  Output calls of variables goes here
   call outfld( 'WP2_CLUBB',        wp2,                     pcols, lchnk )
   call outfld( 'UP2_CLUBB',        up2,                     pcols, lchnk )
   call outfld( 'VP2_CLUBB',        vp2,                     pcols, lchnk )

   call outfld( 'UPWP_CLUBB',       upwp,                    pcols, lchnk )
   call outfld( 'VPWP_CLUBB',       vpwp,                    pcols, lchnk )

   tmp_array = rtp2(:ncol,:)*1000._r8
   call outfld( 'RTP2_CLUBB',       tmp_array,               ncol,  lchnk )
   call outfld( 'THLP2_CLUBB',      thlp2,                   pcols, lchnk )

   tmp_array = rcm(:ncol,:)*1000._r8
   call outfld( 'RCM_CLUBB',        tmp_array,               ncol,  lchnk )
   tmp_array = wprcp(:ncol,:)*latvap
   call outfld( 'WPRCP_CLUBB',      tmp_array,               ncol,  lchnk )
   call outfld( 'CLOUDFRAC_CLUBB',  alst,                    pcols, lchnk )
   tmp_array = rcm_in_layer(:ncol,:)*1000._r8
   call outfld( 'RCMINLAYER_CLUBB', tmp_array,               ncol,  lchnk )
   call outfld( 'CLOUDCOVER_CLUBB', cloud_frac,              pcols, lchnk )

   tmp_array = 1._r8*zt_out(:ncol,:)
   call outfld( 'ZT_CLUBB',         tmp_array,               ncol,  lchnk )
   tmp_array = 1._r8*zi_out(:ncol,:)
   call outfld( 'ZM_CLUBB',         tmp_array,               ncol,  lchnk )
   call outfld( 'UM_CLUBB',         um,                      pcols, lchnk )
   call outfld( 'VM_CLUBB',         vm,                      pcols, lchnk )

