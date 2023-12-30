

#include "clubb_terms_not_in_eqns.inc"

         !Balli- to do: check whether initent-ins and intent-inouts are actually what they say

         call advance_clubb_core_api &
              ( l_implemented, dtime, fcoriolis, sfc_elevation, hydromet_dim, & ! intent(in)
              !
              core_forcing%thlm, &
              core_forcing%rtm,  &
              core_forcing%um,   &
              core_forcing%vm,   &! intent(in)
              sclrm_forcing,     &
              core_forcing%edsclr,  &
              core_forcing%wprtp,   &! intent(in)
              core_forcing%wpthlp,  &
              core_forcing%rtp2,    &
              core_forcing%thlp2,   &! intent(in)
              core_forcing%rtpthlp, &
              !
              core_auxil%wm_zm, &
              core_auxil%wm_zt, &! intent(in)
              !
              core_sfc%wpthlp,    &
              core_sfc%wprtp,     &
              core_sfc%upwp,      &
              core_sfc%vpwp,      &! intent(in)
              wpsclrp_sfc,        &
              core_sfc%wpedsclrp, &! intent(in)
              !
              core_auxil%p_in_Pa,         &
              core_auxil%rho_zm,          &
              core_auxil%rho_zt,          &
              core_auxil%exner,           &! intent(in)
              core_auxil%rho_ds_zm,       &
              core_auxil%rho_ds_zt,       &
              core_auxil%invrs_rho_ds_zm, &! intent(in)
              core_auxil%invrs_rho_ds_zt, &
              core_auxil%thv_ds_zm,       &
              core_auxil%thv_ds_zt,       &
              !
              hydromet,    &! intent(in)
              rfrzm, radf, &! intent(in)
#ifdef CLUBBND_CAM
              varmu2,      &! intent(in)
#endif
              wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &! intent(in)
              host_dx, host_dy, &                          ! intent(in)
              !
              core_prog%um,     &
              core_prog%vm,     &
              core_prog%upwp,   &! intent(inout)
              core_prog%vpwp,   &
              core_prog%up2,    &
              core_prog%vp2,    &! intent(inout)
              core_prog%thlm,   &
              core_prog%rtm,    &
              core_prog%wprtp,  &
              core_prog%wpthlp, &! intent(inout)
              core_prog%wp2,    &
              core_prog%wp3,    &
              core_prog%rtp2,   &! intent(inout)
              rtp3_in,          &
              core_prog%thlp2,  &
              thlp3_in,         &
              core_prog%rtpthlp,&! intent(inout)
              sclrm,                       &! intent(inout)
              sclrp2, sclrprtp, sclrpthlp, &! intent(inout)
              wpsclrp, &
              core_prog%edsclr, &
              err_code,         &! intent(inout)
              !
              core_diag%rcm,        &
              core_diag%cloud_frac, &! intent(inout)
              !
              core_diag%wpthvp,  &
              core_diag%wp2thvp, &
              core_diag%rtpthvp, &
              core_diag%thlpthvp,&! intent(inout)
              !
              sclrpthvp_inout, &! intent(inout)
              pdf_params, pdf_params_zm, &! intent(inout)
              !
              core_diag%khzm,    &
              core_diag%khzt,    &
              core_diag%qclvar,  &
              core_diag%thlprcp, &! intent(out)
              core_diag%wprcp,   &
              !
              ice_supersat_frac,                 &! intent(out)
              rcm_in_layer_out, cloud_cover_out, &! intent(out)
              upwp_sfc_pert, vpwp_sfc_pert,      &! intent(in)
              um_pert_col, vm_pert_col,          &! 
              upwp_pert_col, vpwp_pert_col)       ! intent(inout)
         call t_stopf('advance_clubb_core')

         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Fatal error in CLUBB: at timestep ", get_nstep(), &
                 "LAT (radians): ", host_mean%lat(i), &
                 "LON (radians): ", host_mean%lon(i), &
                 "LAT (degrees): ", rad_to_deg*host_mean%lat(i), &
                 "LON (degrees): ", rad_to_deg*host_mean%lon(i), &
                 "Global Column Number: ", get_gcol_p(lchnk,i)
            call endrun('clubb_tend_cam:  Fatal error in CLUBB library'//errmsg(__FILE__,__LINE__))
         end if


