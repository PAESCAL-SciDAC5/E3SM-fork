 subroutine clubb_related_diags_on_host_grid( pbuf, ptend_loc, host_mean, host_mnts, c2h,    &! in
                                              qclvar, deep_scheme, clubb_do_deep, varmu,     &! in
                                              l_pbl_winds_diag, sfc_v_diff_tau, pert_tau,    &! in
                                              apply_const,                                   &! in
                                              zi_out, zt_out, rcm_in_layer, wprcp            )! in

   use physics_buffer,        only: pbuf_get_index, pbuf_get_field, physics_buffer_desc 
   use physics_types,         only: physics_ptend
   use clubb_intr_host_types, only: clubb_mean_2d_t, clubb_mnts_2d_t, clubb_to_host_t
   use constituents,          only: cnst_get_ind

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_ptend),  intent(in) :: ptend_loc
   type(clubb_mean_2d_t),intent(in) :: host_mean
   type(clubb_mnts_2d_t),intent(in) :: host_mnts
   type(clubb_to_host_t),intent(in) :: c2h

   real(r8) :: qclvar(pcols,pverp)              ! cloud water variance                          [kg^2/kg^2]
   character(len=*), intent(in) :: deep_scheme

   logical, intent(in) :: clubb_do_deep
   real(r8),intent(in) :: varmu(pcols)

   logical, intent(in) :: l_pbl_winds_diag
   real(r8),intent(in) :: sfc_v_diff_tau(pcols)
   real(r8),intent(in) :: pert_tau 

   real(r8),intent(in) :: apply_const 

   real(r8),intent(in) :: zt_out(pcols,pverp)              ! output for the thermo CLUBB grid              [m]
   real(r8),intent(in) :: zi_out(pcols,pverp)              ! output for momentum CLUBB grid                [m]
   real(r8),intent(in) :: rcm_in_layer(pcols,pverp)        ! CLUBB in-cloud liquid water mixing ratio      [kg/kg]
   real(r8),intent(in) :: wprcp(pcols,pverp)               ! CLUBB liquid water flux                       [m/s kg/kg]

   
   real(r8), pointer, dimension(:,:) :: cmeliq 
   real(r8), pointer, dimension(:,:) :: accre_enhan ! accretion enhancement factor              [-]
   real(r8), pointer, dimension(:,:) :: tke      ! turbulent kinetic energy                     [m^2/s^2]

   real(r8)  relvarc(pcols,pver)
   real(r8), pointer, dimension(:,:) :: relvar   ! relative cloud water variance                [-]
   real(r8) :: relvarmax,relvarmin

   real(r8) ::      tmp_array(pcols,pverp)
   real(r8) ::    wpthvp_diag(pcols,pverp)       ! CLUBB buoyancy flux                           [W/m^2]
   real(r8) ::     wp3_output(pcols,pverp)       ! wp3 output                                    [m^3/s^3]
   real(r8) :: rtpthlp_output(pcols,pverp)       ! rtpthlp ouptut                                [K kg/kg]

   integer :: ncol, lchnk
   integer :: ixq, ixcldliq, ixcldice
   integer :: i

   !---------------------------
   ncol  = host_mean%ncol
   lchnk = host_mean%lchnk

   call cnst_get_ind('Q',ixq) 
   call cnst_get_ind('CLDLIQ',ixcldliq)
   call cnst_get_ind('CLDICE',ixcldice)
   !---------------------------
   ! Output various tendencies
   !---------------------------
   call outfld( 'RVMTEND_CLUBB', ptend_loc%q(:,:,ixq), pcols, lchnk)
   call outfld( 'RCMTEND_CLUBB', ptend_loc%q(:,:,ixcldliq), pcols, lchnk)
   call outfld( 'RIMTEND_CLUBB', ptend_loc%q(:,:,ixcldice), pcols, lchnk)
   call outfld( 'TTEND_CLUBB',   ptend_loc%s/cpair,pcols, lchnk)
   call outfld( 'UTEND_CLUBB',   ptend_loc%u,pcols, lchnk)
   call outfld( 'VTEND_CLUBB',   ptend_loc%v,pcols, lchnk)

   call pbuf_get_field(pbuf, pbuf_get_index('CMELIQ'), cmeliq)
   cmeliq(:,:) = ptend_loc%q(:,:,ixcldliq)
   call outfld( 'CMELIQ',        cmeliq, pcols, lchnk)

   !------------------------------------------------------------
   ! Diagnose TKE. Save in pbuf so that it can be used by other
   ! parameterizations (e.g., aerosol activation)
   !------------------------------------------------------------
   call pbuf_get_field(pbuf, tke_idx,     tke)
   tke(1:ncol,1:pverp) = 0.5_r8*( host_mnts% up2(:ncol,:pverp) &
                                 +host_mnts% vp2(:ncol,:pverp) &
                                 +host_mnts% wp2(:ncol,:pverp) )  !  turbulent kinetic energy

   !--------------------------------------------------------------------------------
   ! Diagnose relative liquid water variance; save to pbuf for cloud microphysics
   !--------------------------------------------------------------------------------
   ! input: rcm, qcvar, cloud_frac
   ! output:  relvar (in pbuf), relvarc (outfld)

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
              c2h%cloud_frac(:ncol,:pver)*qclvar(:ncol,:pver)-  &
              (1._r8-c2h%cloud_frac(:ncol,:pver))* host_mean%rcm(:ncol,:pver)**2)))
              relvarc(:ncol,:pver) = min(relvarmax,max(relvarmin, host_mean%rcm(:ncol,:pver)**2/max(qsmall,  &
              c2h%cloud_frac(:ncol,:pver)*qclvar(:ncol,:pver)-  &
              (1._r8-c2h%cloud_frac(:ncol,:pver))* host_mean%rcm(:ncol,:pver)**2)))
      else

         where ( host_mean%rcm(:ncol,:pver) /= 0 .and. qclvar(:ncol,:pver) /= 0) &
              relvar(:ncol,:pver) = min(relvarmax,max(0.001_r8, host_mean%rcm(:ncol,:pver)**2/qclvar(:ncol,:pver)))
      endif
   endif
   call outfld( 'RELVAR',           relvar,                  pcols, lchnk )
   call outfld( 'RELVARC',          relvarc,                 pcols, lchnk )

   !----------------------------------------
   ! Optional Accretion enhancement factor
   !----------------------------------------
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)
   accre_enhan(:ncol,:pver) = micro_mg_accre_enhan_fac !default is 1._r8

   !----------------------------------------
   !----------------------------------------
   if (clubb_do_deep) call outfld( 'MU_CLUBB',      varmu      ,pcols, lchnk)

   ! Linearization of PBL winds: diagnose wsresp and tau_est for history output.
   ! The outfld calls are in cam_diagnostics.F90

   if (l_pbl_winds_diag) then
      do i = 1, ncol
         c2h% wsresp(i) = sfc_v_diff_tau(i) / pert_tau
         ! Estimated tau in balance with wind is the tau we just used.
         if (host_mean%wsx(i) == 0._r8 .or. host_mean%wsy(i) == 0._r8) then
            ! Work around an odd FPE issue with intel compiler.
            c2h% tau_est(i) = abs(host_mean%wsx(i)) + abs(host_mean%wsy(i))
         else
            c2h% tau_est(i) = hypot(host_mean%wsx(i), host_mean%wsy(i))
         end if
      end do
   end if

   
   !----------------------------------------
   ! Miscellaneous diagnostics from CLUBB
   !----------------------------------------
   call outfld( 'UM_CLUBB',         host_mean%um,            pcols, lchnk )
   call outfld( 'VM_CLUBB',         host_mean%vm,            pcols, lchnk )

   call outfld( 'THLP2_CLUBB',      host_mnts% thlp2,        pcols, lchnk )
   call outfld( 'WP2_CLUBB',        host_mnts% wp2,          pcols, lchnk )
   call outfld( 'UP2_CLUBB',        host_mnts% up2,          pcols, lchnk )
   call outfld( 'VP2_CLUBB',        host_mnts% vp2,          pcols, lchnk )

   call outfld( 'UPWP_CLUBB',       host_mnts% upwp,         pcols, lchnk )
   call outfld( 'VPWP_CLUBB',       host_mnts% vpwp,         pcols, lchnk )

   call outfld( 'CLOUDFRAC_CLUBB',  c2h% cloud_frac,         pcols, lchnk )

   tmp_array = host_mnts%rtp2(:ncol,:)*1000._r8
   call outfld( 'RTP2_CLUBB',       tmp_array,               pcols, lchnk )

   tmp_array = host_mean%rcm(:ncol,:)*1000._r8
   call outfld( 'RCM_CLUBB',        tmp_array,               pcols, lchnk )

   tmp_array = wprcp(:ncol,:)*latvap
   call outfld( 'WPRCP_CLUBB',      tmp_array,               pcols, lchnk )

   tmp_array = rcm_in_layer(:ncol,:)*1000._r8
   call outfld( 'RCMINLAYER_CLUBB', tmp_array,               pcols, lchnk )

   tmp_array = 1._r8*zt_out(:ncol,:)
   call outfld( 'ZT_CLUBB',         tmp_array,               pcols, lchnk )

   tmp_array = 1._r8*zi_out(:ncol,:)
   call outfld( 'ZM_CLUBB',         tmp_array,               pcols, lchnk )

   ! ------------------------------------------------- !
   ! Diagnose some output variables                    !
   ! ------------------------------------------------- !

   !-------------
   ! k = 1,pverp
   !-------------
   ! wp3 output
   wp3_output(:ncol,:pverp) = host_mnts%wp3(:ncol,:pverp) - (apply_const*wp3_const)
   call outfld( 'WP3_CLUBB',        wp3_output,              pcols, lchnk )

   ! rtpthlp output
   rtpthlp_output(:ncol,:pverp) = host_mnts%rtpthlp(:ncol,:pverp)-(apply_const*rtpthlp_const)
   tmp_array = rtpthlp_output(:ncol,:)*1000._r8
   call outfld( 'RTPTHLP_CLUBB',    tmp_array,               pcols, lchnk )

 end subroutine clubb_related_diags_on_host_grid
