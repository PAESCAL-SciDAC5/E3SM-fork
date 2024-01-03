   !-------------------------------
   ! Diagnose some output variables
   !-------------------------------
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

