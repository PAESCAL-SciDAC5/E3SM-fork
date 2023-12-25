   !------
   ! TKE
   !------
   call pbuf_get_field(pbuf, tke_idx,     tke)
   tke(1:ncol,1:pverp) = 0.5_r8*(up2(:ncol,:pverp)+vp2(:ncol,:pverp)+wp2(:ncol,:pverp))  !  turbulent kinetic energy

   !--------------------------------
   ! relative liquid water variance
   !--------------------------------
#include "relvar.inc"
   ! input: rcm, qcvar, cloud_frac
   ! output:  relvar (in pbuf), relvarc (outfld)
   !----------------------------------------
   ! Optional Accretion enhancement factor
   !----------------------------------------
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)
    accre_enhan(:ncol,:pver) = micro_mg_accre_enhan_fac !default is 1._r8




   call outfld( 'CLOUDFRAC_CLUBB',  cloud_frac,              pcols, lchnk )
   call outfld( 'WP2_CLUBB',        wp2,                     pcols, lchnk )
   call outfld( 'UP2_CLUBB',        up2,                     pcols, lchnk )
   call outfld( 'VP2_CLUBB',        vp2,                     pcols, lchnk )

   call outfld( 'UPWP_CLUBB',       upwp,                    pcols, lchnk )
   call outfld( 'VPWP_CLUBB',       vpwp,                    pcols, lchnk )

   tmp_array = rtp2(:ncol,:)*1000._r8
   call outfld( 'RTP2_CLUBB',       tmp_array,               ncol,  lchnk )
   call outfld( 'THLP2_CLUBB',      thlp2,                   pcols, lchnk )

   tmp_array = host_mean%rcm(:ncol,:)*1000._r8
   call outfld( 'RCM_CLUBB',        tmp_array,               ncol,  lchnk )

   tmp_array = wprcp(:ncol,:)*latvap
   call outfld( 'WPRCP_CLUBB',      tmp_array,               ncol,  lchnk )

   tmp_array = rcm_in_layer(:ncol,:)*1000._r8
   call outfld( 'RCMINLAYER_CLUBB', tmp_array,               ncol,  lchnk )

   tmp_array = 1._r8*zt_out(:ncol,:)
   call outfld( 'ZT_CLUBB',         tmp_array,               ncol,  lchnk )

   tmp_array = 1._r8*zi_out(:ncol,:)
   call outfld( 'ZM_CLUBB',         tmp_array,               ncol,  lchnk )

   call outfld( 'UM_CLUBB',         host_mean%um,            pcols, lchnk )
   call outfld( 'VM_CLUBB',         host_mean%vm,            pcols, lchnk )

   ! ------------------------------------------------- !
   ! Diagnose some output variables                    !
   ! ------------------------------------------------- !
   !  density
   rho(:ncol,1:pver) = state1%pmid(:ncol,1:pver)/(rair*state1%t(:ncol,1:pver))
   rho(:ncol,pverp)  = state1%ps(:ncol)/(rair*state1%t(:ncol,pver))
   call outfld( 'RHO_CLUBB',        rho,                     pcols, lchnk )

   !-------------
   ! k = 1,pver
   !-------------
   !  buoyancy flux
   eps = rair/rh2o

   wpthvp_diag(:,:) = 0.0_r8
   wpthvp_diag(:ncol,:pver) = (wpthlp(:ncol,:pver)-(apply_const*wpthlp_const))+((1._r8-eps)/eps)*real(theta0, kind = r8)* &
                               (wprtp(:ncol,:pver)-(apply_const*wprtp_const))+((latvap/cpair)* &
                         state1%exner(:ncol,:pver)-(1._r8/eps)*real(theta0, kind = r8))*wprcp(:ncol,:pver)
   tmp_array = wpthvp_diag(:ncol,:)*cpair
   call outfld( 'WPTHVP_CLUBB',     tmp_array,               ncol,  lchnk )

   ! total water mixing ratio
   qt_output(:ncol,:pver) = state1%q(:ncol,:pver,ixq)+state1%q(:ncol,:pver,ixcldliq)+state1%q(:ncol,:pver,ixcldice)
   call outfld( 'QT',               qt_output,               pcols, lchnk )

   ! liquid water potential temperature
   thetal_output(:ncol,:pver) = (state1%t(:ncol,:pver)*state1%exner(:ncol,:pver))-(latvap/cpair)*state1%q(:ncol,:pver,ixcldliq)
   call outfld( 'THETAL',           thetal_output,           pcols, lchnk )

   ! liquid water static energy
   sl_output(:ncol,:pver) = cpair*state1%t(:ncol,:pver)+gravit*state1%zm(:ncol,:pver)-latvap*state1%q(:ncol,:pver,ixcldliq)
   call outfld( 'SL',               sl_output,               pcols, lchnk )

   !-------------
   ! k = 1,pverp
   !-------------

   ! liquid water potential temperature flux
   wpthlp_output(:ncol,:pverp)  = (wpthlp(:ncol,:pverp)-(apply_const*wpthlp_const))*rho(:ncol,:pverp)*cpair
   call outfld( 'WPTHLP_CLUBB',     wpthlp_output,           pcols, lchnk )
 
   ! wp3 output
   wp3_output(:ncol,:pverp) = wp3(:ncol,:pverp) - (apply_const*wp3_const)
   call outfld( 'WP3_CLUBB',        wp3_output,              pcols, lchnk )

   ! rtpthlp output
   rtpthlp_output(:ncol,:pverp) = rtpthlp(:ncol,:pverp)-(apply_const*rtpthlp_const)
   tmp_array = rtpthlp_output(:ncol,:)*1000._r8
   call outfld( 'RTPTHLP_CLUBB',    tmp_array,               ncol,  lchnk )

   ! total water mixig ratio flux
   wprtp_output(:ncol,:pverp) = (wprtp(:ncol,:pverp)-(apply_const*wprtp_const))*rho(:ncol,:pverp)*latvap
   call outfld( 'WPRTP_CLUBB',      wprtp_output,            pcols, lchnk )


