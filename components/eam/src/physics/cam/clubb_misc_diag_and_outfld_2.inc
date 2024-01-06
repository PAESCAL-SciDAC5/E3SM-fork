   !-------------------------------
   ! Diagnose some output variables
   !-------------------------------
   !-------------
   ! k = 1,pverp     
   !-------------
   !  density
   rho(:ncol,1:pver) = state1%pmid(:ncol,1:pver)/(rair*state1%t(:ncol,1:pver))
   rho(:ncol,pverp)  = state1%ps(:ncol)/(rair*state1%t(:ncol,pver))
   call outfld( 'RHO_CLUBB',        rho,                     pcols, lchnk )

   ! liquid water potential temperature flux
   wpthlp_output(:ncol,:pverp)  = (host_mnts%wpthlp(:ncol,:pverp)-(apply_const*wpthlp_const))*rho(:ncol,:pverp)*cpair
   call outfld( 'WPTHLP_CLUBB',     wpthlp_output,           pcols, lchnk )
 
   ! total water mixig ratio flux
   wprtp_output(:ncol,:pverp) = (host_mnts%wprtp(:ncol,:pverp)-(apply_const*wprtp_const))*rho(:ncol,:pverp)*latvap
   call outfld( 'WPRTP_CLUBB',      wprtp_output,            pcols, lchnk )


   !-------------
   ! k = 1,pver     
   !-------------
   !  buoyancy flux 
   eps = rair/rh2o
   
   wpthvp_diag(:,:) = 0.0_r8
   wpthvp_diag(:ncol,:pver) = ( host_mnts% wpthlp(:ncol,:pver)-(apply_const*wpthlp_const) ) &
                             +  ((1._r8-eps)/eps)*real(theta0, kind = r8)                   &
                               *(host_mnts% wprtp(:ncol,:pver)-(apply_const*wprtp_const))   &
                             +  ((latvap/cpair)*state1%exner(:ncol,:pver)-(1._r8/eps)*real(theta0, kind = r8)) &
                               *wprcp(:ncol,:pver)
   tmp_array = wpthvp_diag(:ncol,:)*cpair
   call outfld( 'WPTHVP_CLUBB',     tmp_array,               pcols, lchnk )
   
   ! total water mixing ratio
   qt_output(:ncol,:pver) = state1%q(:ncol,:pver,ixq)+state1%q(:ncol,:pver,ixcldliq)+state1%q(:ncol,:pver,ixcldice)
   call outfld( 'QT',               qt_output,               pcols, lchnk )
   
   ! liquid water potential temperature
   thetal_output(:ncol,:pver) = (state1%t(:ncol,:pver)*state1%exner(:ncol,:pver))-(latvap/cpair)*state1%q(:ncol,:pver,ixcldliq)
   call outfld( 'THETAL',           thetal_output,           pcols, lchnk )
   
   ! liquid water static energy
   sl_output(:ncol,:pver) = cpair*state1%t(:ncol,:pver)+gravit*state1%zm(:ncol,:pver)-latvap*state1%q(:ncol,:pver,ixcldliq)
   call outfld( 'SL',               sl_output,               pcols, lchnk )

