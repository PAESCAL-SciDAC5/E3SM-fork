
   !-----------------------------
   ! Initialize ptend for CLUBB
   !-----------------------------
   call physics_ptend_init(ptend_loc, state1%psetcols, 'clubb', ls=.true., lu=.true., lv=.true., lq=lq)

   !-------------------------------------------------------------
   ! Load input (initial and boundary conditions) from host model's data structures
   !-------------------------------------------------------------
   call  gather_mean_fields( state1, itim_old, pbuf, host_mean ) !in, in, inout, out

#include "clubb_timestep_init_moments.inc"
! input: state or pbuf
! output: 9 moments
!-------------------------------------------------------------
#include "clubb_timestep_init_forcing_and_boundary.inc"
#include "clubb_misc_inout_fields.inc"


   call t_startf('adv_clubb_core_col_loop')
   !----------------------------------------------------------
   !  Loop over all columns in lchnk to advance CLUBB core
   !----------------------------------------------------------
   do i=1,ncol   ! loop over columns

      ! Calculate column total energy before CLUBB. The value is saved
      ! for the energy fixer at the end of the column loop.
 
      call column_total_energy( state1, cam_in, i, hdtime, &! in
                                te_b                       )! out, scalar

      ! CLUBB can be substepped with respect to the host model's timestep hdtime.
      ! Determine the substepsize and number of substeps for CLUBB.
      ! Note that in the standard EAM which uses the same vertical grid for all processes,
      ! this call only needs to be done once at the beginning of an integration.
      ! But we place the call here (meaning repeating the calculation every timestep and
      ! in each grid column), as we are preparing for the implementation of AMR.

      call determine_clubb_dtime( hdtime, clubb_timestep, hdtime_core_rknd, dtime, n_clubb_core_step ) ! 2xin, 3xout

#include "clubb_tend_cam_1column.inc"

      ! After advancing CLUBB in time, apply a total energy fixer by adjusting the
      ! ending s (dry static energy).

      call column_total_energy_fixer( te_b, host_mean%thlm, host_mean%rtm, host_mean%rcm, &
                                      host_mean%exner_clubb,          &! in
                                      host_mean%um, host_mean%vm, wp2, state1, i, &! in
                                      clubb_s                                     )! out
#include "clubb_ptend_cal.inc"

   enddo  ! end column loop
   call t_stopf('adv_clubb_core_col_loop')

   deallocate( host_mean%exner_clubb )
   ! ------------------------------------------------------------ !
   ! Add constant to ghost point so that output is not corrupted
   ! (Note that this is only done to k = pverp.)
   ! ------------------------------------------------------------ !
   if (clubb_do_adv) then
      if (macmic_it .eq. cld_macmic_num_steps) then
         wp3(:,pverp) = wp3(:,pverp) + wp3_const
         rtpthlp(:,pverp) = rtpthlp(:,pverp) + rtpthlp_const
         wpthlp(:,pverp) = wpthlp(:,pverp) + wpthlp_const
         wprtp(:,pverp) = wprtp(:,pverp) + wprtp_const
      endif
   endif

   !============================
   ! Linearization of PBL winds
   !============================
   if (linearize_pbl_winds .and. macmic_it == cld_macmic_num_steps) then
      do i = 1, ncol
         wsresp(i) = sfc_v_diff_tau(i) / pert_tau
         ! Estimated tau in balance with wind is the tau we just used.
         if (cam_in%wsx(i) == 0._r8 .or. cam_in%wsy(i) == 0._r8) then
            ! Work around an odd FPE issue with intel compiler.
            tau_est(i) = abs(cam_in%wsx(i)) + abs(cam_in%wsy(i))
         else
            tau_est(i) = hypot(cam_in%wsx(i), cam_in%wsy(i))
         end if
      end do
   end if

   !=========
   ! Output
   !=========
#include "clubb_stats_output.inc"

   ! ------------------------------------------------- !
   ! End column computation of CLUBB, begin to apply   !
   ! and compute output, etc                           !
   ! ------------------------------------------------- !
   !  Output CLUBB tendencies
   call outfld( 'RVMTEND_CLUBB', ptend_loc%q(:,:,ixq), pcols, lchnk)
   call outfld( 'RCMTEND_CLUBB', ptend_loc%q(:,:,ixcldliq), pcols, lchnk)
   call outfld( 'RIMTEND_CLUBB', ptend_loc%q(:,:,ixcldice), pcols, lchnk)
   call outfld( 'TTEND_CLUBB',   ptend_loc%s/cpair,pcols, lchnk)
   call outfld( 'UTEND_CLUBB',   ptend_loc%u,pcols, lchnk)
   call outfld( 'VTEND_CLUBB',   ptend_loc%v,pcols, lchnk)

   call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)
   cmeliq(:,:) = ptend_loc%q(:,:,ixcldliq)
   call outfld( 'CMELIQ',        cmeliq, pcols, lchnk)

   if (clubb_do_deep) call outfld( 'MU_CLUBB',      varmu      ,pcols, lchnk)
