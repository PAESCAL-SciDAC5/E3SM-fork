
  !===============================================================================
  !
  !===============================================================================
  subroutine gather_host_fields( cld_macmic_num_steps, macmic_it, itim_old,      &
                                 state1, pbuf, host_mean, host_mnts, apply_const )

   use constituents,    only: cnst_get_ind
   use physics_types,   only: physics_state
   use physics_buffer,  only: physics_buffer_desc, pbuf_get_field
   use clubb_intr_types,only: clubb_mean_2d_t, clubb_mnts_2d_t

   use constituents,    only: pcnst, cnst_type
   use co2_cycle,       only: co2_cycle_set_cnst_type
   use physics_types,   only: set_dry_to_wet

   use time_manager,    only: is_first_step

   integer, intent(in) :: cld_macmic_num_steps     ! number of mac-mic iterations
   integer, intent(in) :: macmic_it                ! number of mac-mic iterations
   
   type(physics_state),target,intent(inout) :: state1
   type(physics_buffer_desc), pointer       :: pbuf(:)
   type(clubb_mean_2d_t)                 :: host_mean
   type(clubb_mnts_2d_t)                 :: host_mnts
   real(r8),intent(out)                  :: apply_const

   integer,intent(in) :: itim_old

   integer :: ncol
   integer :: ixq, ixcldliq
   character(len=3), dimension(pcnst) :: cnst_type_loc ! local override option for constituents cnst_type

   ncol            = state1%ncol
   host_mean%ncol  = state1%ncol
   host_mean%lchnk = state1%lchnk

   host_mean%lat   => state1%lat
   host_mean%lon   => state1%lon

   !===================================================================================================
   !  At each CLUBB call, initialize mean momentum and thermo CLUBB state from the host model's state
   !===================================================================================================
   ! From host model: u,  v;  qv, ql, t, pmid
   ! To CLUBB:        um, vm; rtm, thlm, rcm 
   !--------------------------------------------------------

   call pbuf_get_field(pbuf, um_idx, host_mean%um, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vm_idx, host_mean%vm, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   host_mean%um(:ncol,:pver) = state1%u(:ncol,:pver)
   host_mean%vm(:ncol,:pver) = state1%v(:ncol,:pver)

   !---------------------------
   call pbuf_get_field(pbuf, thlm_idx, host_mean%thlm, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf,  rtm_idx, host_mean%rtm,  start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf,  rcm_idx, host_mean%rcm,  start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)

   host_mean%rtm(:ncol,:pver) = state1%q(:ncol,:pver,ixq)+state1%q(:ncol,:pver,ixcldliq)

   !  Compute exner function consistent with CLUBB's definition, which uses a constant
   !  surface pressure.  CAM's exner (in state does not).  Therefore, for consistent
   !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables
   !  (such as thlm), use host_mean%exner_clubb; otherwise use the exner in state

   allocate( host_mean%exner_clubb(pcols,pverp) )
   host_mean%exner_clubb(:ncol,:pver) = (real(p0_clubb, kind = r8 )/state1%pmid(:ncol,:pver))**(rair/cpair)
   host_mean%thlm(:ncol,:pver) = state1%t(:ncol,:pver)*host_mean%exner_clubb(:ncol,:pver) &
                                 -(latvap/cpair)*state1%q(:ncol,:pver,ixcldliq)

   !---------------------------------
   ! Assign value to ghost level - this can be done later?
   !---------------------------------

   host_mean%um(1:ncol,pverp) = host_mean%um(1:ncol,pver)
   host_mean%vm(1:ncol,pverp) = host_mean%vm(1:ncol,pver)

   host_mean%rtm(1:ncol,pverp)  = host_mean%rtm(1:ncol,pver)
   host_mean%thlm(1:ncol,pverp) = host_mean%thlm(1:ncol,pver)

   !===================================================================================================
   ! Associate pointers to additional variables in host model's state vector that are needed below
   !===================================================================================================
   host_mean%zi    => state1%zi
   host_mean%zm    => state1%zm
   host_mean%pmid  => state1%pmid
   host_mean%pdel  => state1%pdel
   host_mean%t     => state1%t
   host_mean%omega => state1%omega
   host_mean%q     => state1%q

   
   !=====================================================
   ! Passive tracers:
   ! Constituents are all treated as wet mmr by clubb
   ! don't convert co2 tracers to wet mixing ratios
   !=====================================================
   cnst_type_loc(:) = cnst_type(:)
   call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
   call set_dry_to_wet(state1, cnst_type_loc)

   !===========================================================================================
   !  At each CLUBB call, initialize CLUBB's PROGNOSTIC moments
   !===========================================================================================
   ! Associate pointers to arrays in pbuf
   !-------------------------------------------------------------------------------------------

   call pbuf_get_field(pbuf, wp2_idx,     host_mnts%wp2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wp3_idx,     host_mnts%wp3,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wpthlp_idx,  host_mnts%wpthlp,  start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wprtp_idx,   host_mnts%wprtp,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtpthlp_idx, host_mnts%rtpthlp, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtp2_idx,    host_mnts%rtp2,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, thlp2_idx,   host_mnts%thlp2,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, up2_idx,     host_mnts%up2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vp2_idx,     host_mnts%vp2,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   !----------------------------------------------------------------------------------------------------
   ! If CLUBB's moments are transported by the host model, copy the moments from the host's state vector
   !----------------------------------------------------------------------------------------------------
   if (.not.clubb_do_adv) then ! no need to apply constant if CLUBB's moments are not advected

     apply_const = 0._r8

   else !-----------------------------------------------

     ! Intialize the apply_const variable (note special logic is due to eularian backstepping)
     if (is_first_step() .or. all(host_mnts%wpthlp(1:ncol,1:pver) .eq. 0._r8)) then
        apply_const = 0._r8  ! On first time through do not remove constant
                             !  from moments since it has not been added yet
     else
        apply_const = 1._r8  ! Initialize to one, only if CLUBB's moments are advected
     endif

     if (macmic_it==1) then
        !  This is the first macmic substep. Get the moments from the host model.
        !  Note that some of the moments below can be positive or negative.
        !    Remove a constant that was added to prevent dynamics from clipping
        !    them to prevent dynamics from making them positive.
        host_mnts%   thlp2(:ncol,:pver) = state1%q(:ncol,:pver,ixthlp2)
        host_mnts%    rtp2(:ncol,:pver) = state1%q(:ncol,:pver,ixrtp2)
        host_mnts% rtpthlp(:ncol,:pver) = state1%q(:ncol,:pver,ixrtpthlp) - (rtpthlp_const*apply_const)
        host_mnts%  wpthlp(:ncol,:pver) = state1%q(:ncol,:pver,ixwpthlp)  - (wpthlp_const*apply_const)
        host_mnts%   wprtp(:ncol,:pver) = state1%q(:ncol,:pver,ixwprtp)   - (wprtp_const*apply_const)
        host_mnts%     wp2(:ncol,:pver) = state1%q(:ncol,:pver,ixwp2)
        host_mnts%     wp3(:ncol,:pver) = state1%q(:ncol,:pver,ixwp3)     - (wp3_const*apply_const)
        host_mnts%     up2(:ncol,:pver) = state1%q(:ncol,:pver,ixup2)
        host_mnts%     vp2(:ncol,:pver) = state1%q(:ncol,:pver,ixvp2)
     endif

     ! Assign value to ghost level
     host_mnts%   thlp2(1:ncol,pverp) = host_mnts%   thlp2(1:ncol,pver)
     host_mnts%    rtp2(1:ncol,pverp) = host_mnts%    rtp2(1:ncol,pver)
     host_mnts% rtpthlp(1:ncol,pverp) = host_mnts% rtpthlp(1:ncol,pver)
     host_mnts%  wpthlp(1:ncol,pverp) = host_mnts%  wpthlp(1:ncol,pver)
     host_mnts%   wprtp(1:ncol,pverp) = host_mnts%   wprtp(1:ncol,pver)
     host_mnts%     wp2(1:ncol,pverp) = host_mnts%     wp2(1:ncol,pver)
     host_mnts%     wp3(1:ncol,pverp) = host_mnts%     wp3(1:ncol,pver)
     host_mnts%     up2(1:ncol,pverp) = host_mnts%     up2(1:ncol,pver)
     host_mnts%     vp2(1:ncol,pverp) = host_mnts%     vp2(1:ncol,pver)

     ! If not last step of macmic loop then set apply_const back to
     !   zero to prevent output from being corrupted.
     if (macmic_it .eq. cld_macmic_num_steps) then
       apply_const = 1._r8
     else
       apply_const = 0._r8
     endif

   endif  ! if (not.clubb_do_adv) !-----------------------------------------------

   !=====================
   ! upwp and vpwp
   !=====================
   call pbuf_get_field(pbuf, upwp_idx, host_mnts%upwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vpwp_idx, host_mnts%vpwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   !=====================
   ! Diagnostic moments
   !=====================
   call pbuf_get_field(pbuf, wpthvp_idx,   host_mnts%wpthvp,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, wp2thvp_idx,  host_mnts%wp2thvp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, rtpthvp_idx,  host_mnts%rtpthvp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, thlpthvp_idx, host_mnts%thlpthvp,   start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   !=====================
   ! Miscellaneous fields
   !=====================

   call pbuf_get_field(pbuf, prer_evap_idx, host_mean%prer_evap)

   if (do_cldcool) then
      call pbuf_get_field(pbuf, qrl_idx, host_mean%qrl)
   end if

  end subroutine gather_host_fields

  !===============================================================================
  !===============================================================================
  subroutine setup_clubb_to_host_fields( pbuf, c2h, itim_old )

   use physics_buffer,  only: physics_buffer_desc, pbuf_get_field
   use clubb_intr_types,only: clubb_to_host_t

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(clubb_to_host_t)              :: c2h
   integer,intent(in)                 :: itim_old

   call pbuf_get_field( pbuf, kvm_idx, c2h%khzt )
   call pbuf_get_field( pbuf, kvh_idx, c2h%khzm )

   call pbuf_get_field( pbuf, cloud_frac_idx,        c2h% cloud_frac       , start=(/1,1,itim_old/), kount=(/pcols,pverp,1/) )

   call pbuf_get_field( pbuf, pdf_zm_w_1_idx,        c2h% pdf_zm_w_1       , start=(/1,1,itim_old/), kount=(/pcols,pverp,1/) )
   call pbuf_get_field( pbuf, pdf_zm_w_2_idx,        c2h% pdf_zm_w_2       , start=(/1,1,itim_old/), kount=(/pcols,pverp,1/) )
   call pbuf_get_field( pbuf, pdf_zm_varnce_w_1_idx, c2h% pdf_zm_varnce_w_1, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/) )
   call pbuf_get_field( pbuf, pdf_zm_varnce_w_2_idx, c2h% pdf_zm_varnce_w_2, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/) )
   call pbuf_get_field( pbuf, pdf_zm_mixt_frac_idx,  c2h% pdf_zm_mixt_frac , start=(/1,1,itim_old/), kount=(/pcols,pverp,1/) )

  end subroutine setup_clubb_to_host_fields


  !===============================================================================
  !===============================================================================
  subroutine map_1column_of_info_from_host_to_clubb( i, host_mean, host_mnts, c2h, &
                                                     pverp, nz, edsclr_dim, lq,ixq, ixcldliq, ixcldice, &
                                                     macmic_it, clubb_do_adv, do_expldiff,      &
                                                     core_auxil, core_prog, core_diag,          &
                                                     core_forcing, core_sfc, clubb_misc,        &
                                                     pdf_params_zm, dz_g, varmu  )

    use shr_kind_mod,     only: r8=>shr_kind_r8
    use clubb_precision,  only: core_rknd

    use clubb_intr_types,      only: clubb_mean_2d_t, clubb_mnts_2d_t, clubb_to_host_t
    use clubb_intr_core_types, only: core_auxil_t, core_prog_t, core_diag_t
    use clubb_intr_core_types, only: core_forcing_t, core_sfc_t, clubb_misc_t
    use clubb_api_module,      only: zt2zm_api

    use physconst,        only: gravit

    use time_manager,     only: is_first_restart_step
    use model_flags,      only: ipdf_call_placement
    use advance_clubb_core_module, only: ipdf_post_advance_fields
      
    type(clubb_mean_2d_t),intent(in) :: host_mean
    type(clubb_mnts_2d_t),intent(in) :: host_mnts
    type(clubb_to_host_t),intent(in) :: c2h

    integer, intent(in) :: i       ! column index 
    integer, intent(in) :: pverp   ! number of vertical levels in column in the host model
    integer, intent(in) :: nz      ! number of vertical levels in column in CLUBB 
    integer, intent(in) :: edsclr_dim 
    logical, intent(in) :: lq(:)
    integer, intent(in) :: ixq, ixcldliq, ixcldice
    integer, intent(in) :: macmic_it
    logical, intent(in) :: clubb_do_adv
    logical, intent(in) :: do_expldiff

    type(core_auxil_t  ),intent(inout) :: core_auxil
    type(core_prog_t   ),intent(inout) :: core_prog
    type(core_diag_t   ),intent(inout) :: core_diag
    type(core_forcing_t),intent(inout) :: core_forcing
    type(core_sfc_t    ),intent(inout) :: core_sfc
    type(clubb_misc_t  ),intent(inout) :: clubb_misc
    type(pdf_parameter), pointer :: pdf_params_zm ! PDF parameters on momentum levs. [units vary]

    real(r8),intent(out) :: dz_g(:)
    real(r8),intent(out) :: varmu

    integer :: pver 
    integer :: icnt, ixind
    integer :: k,kkhost
    real(r8) :: invrs_gravit 
    real(r8) :: thvtmp, dum1

      if (pverp/=nz) call endrun('pverp/=nz not supported by current code')

      pver = pverp-1

      !---------------------------------------------------------------
      ! Grid column info that is used only for debugging currently
      !---------------------------------------------------------------
      core_auxil% icol  = i 
      core_auxil% lat   = host_mean% lat(i)
      core_auxil% lon   = host_mean% lon(i)
      core_auxil% lchnk = host_mean% lchnk

      !------------------------------------------------------------------
      invrs_gravit = 1._r8 / gravit

      do k=1,pver
         dz_g(k) = host_mean%zi(i,k)-host_mean%zi(i,k+1)  ! compute thickness
      enddo

      ! k = pver
      ! dz_g_bot = host_mean%zi(i,k)-host_mean%zi(i,k+1)  ! compute thickness

      !------------------------------------------------------------------
      !  Compute thermodynamic stuff needed for CLUBB on thermo levels.
      !  Inputs for the momentum levels are set below setup_clubb core
      !------------------------------------------------------------------
      do k=1,pver
         kkhost = pver-k+1

         core_auxil%p_in_Pa(k+1)   = real( host_mean%pmid(i,kkhost), kind=core_rknd )
         core_auxil%exner(k+1)     = real( 1._r8/host_mean%exner_clubb(i,kkhost), kind=core_rknd )

         thvtmp =  host_mean%t(i,kkhost) &
                  *host_mean%exner_clubb(i,kkhost) &
                  *( 1._r8+zvir*host_mean%q(i,kkhost,ixq)-host_mean%q(i,kkhost,ixcldliq) )
         core_auxil%thv_ds_zt(k+1) = real(thvtmp, kind = core_rknd)                       ! thetav on thermo
      enddo

      core_auxil%p_in_Pa  (1) = core_auxil%p_in_Pa(2)
      core_auxil%exner    (1) = core_auxil%exner  (2)
      core_auxil%thv_ds_zt(1) = core_auxil%thv_ds_zt(2)

      !----------
      do k=1,pver
         kkhost = pver-k+1
         core_auxil%rho_ds_zt(k+1) = real( invrs_gravit*host_mean%pdel(i,kkhost)/dz_g(kkhost), kind=core_rknd )
      enddo
      core_auxil%rho_ds_zt(1) = core_auxil%rho_ds_zt(2)


      !----------
      core_auxil%rho_zt(:)  = core_auxil%rho_ds_zt(:)

      core_auxil%invrs_rho_ds_zt(2:pverp) = 1._core_rknd/core_auxil%rho_ds_zt(2:pverp) ! Inverse ds rho at thermo
      core_auxil%invrs_rho_ds_zt(1)       = core_auxil%invrs_rho_ds_zt(2)

      !  Compute some inputs from the thermodynamic grid to the momentum grid
      !  Note: this can only be done after "call setup_grid_heights_api"

      core_auxil%rho_ds_zm       = zt2zm_api( core_auxil%rho_ds_zt )
      core_auxil%rho_zm          = zt2zm_api( core_auxil%rho_zt    )
      core_auxil%invrs_rho_ds_zm = zt2zm_api( core_auxil%invrs_rho_ds_zt)
      core_auxil%thv_ds_zm       = zt2zm_api( core_auxil%thv_ds_zt )

      !----------
      !  Compute mean w wind on thermo grid, convert from omega to w

      core_auxil%wm_zt(1) = 0._core_rknd
      do k=1,pver
        kkhost = pver-k+1
        dum1 = -1._r8*host_mean%omega(i,kkhost) * real(core_auxil%invrs_rho_ds_zt(k+1),kind=r8) *invrs_gravit
        core_auxil%wm_zt(k+1) = real(dum1,kind=core_rknd)
      enddo

      core_auxil%wm_zm = zt2zm_api( core_auxil%wm_zt )


      !---------------------------------------------------------------
      !  Need to flip arrays around for CLUBB core
      do k=1,pverp

         core_prog%thlm(k) = real(host_mean%thlm(i,pverp-k+1), kind = core_rknd)
         core_prog%rtm (k) = real(host_mean%rtm(i,pverp-k+1), kind = core_rknd)

         core_prog%um(k)      = real(host_mean%um(i,pverp-k+1), kind = core_rknd)
         core_prog%vm(k)      = real(host_mean%vm(i,pverp-k+1), kind = core_rknd)

         core_prog%up2(k)     = real(host_mnts%up2(i,pverp-k+1), kind = core_rknd)
         core_prog%vp2(k)     = real(host_mnts%vp2(i,pverp-k+1), kind = core_rknd)
         core_prog%wp2(k)     = real(host_mnts%wp2(i,pverp-k+1), kind = core_rknd)
         core_prog%wp3(k)     = real(host_mnts%wp3(i,pverp-k+1), kind = core_rknd)

         core_prog%upwp(k)    = real(host_mnts%upwp(i,pverp-k+1), kind = core_rknd)
         core_prog%vpwp(k)    = real(host_mnts%vpwp(i,pverp-k+1), kind = core_rknd)

         core_prog%rtp2   (k) = real(host_mnts%rtp2(i,pverp-k+1), kind = core_rknd)
         core_prog%thlp2  (k) = real(host_mnts%thlp2(i,pverp-k+1), kind = core_rknd)
         core_prog%wprtp  (k) = real(host_mnts%wprtp(i,pverp-k+1), kind = core_rknd)
         core_prog%wpthlp (k) = real(host_mnts%wpthlp(i,pverp-k+1), kind = core_rknd)
         core_prog%rtpthlp(k) = real(host_mnts%rtpthlp(i,pverp-k+1), kind = core_rknd)

         core_diag%wpthvp  (k) = real(host_mnts%wpthvp  (i,pverp-k+1), kind = core_rknd)
         core_diag%wp2thvp (k) = real(host_mnts%wp2thvp (i,pverp-k+1), kind = core_rknd)
         core_diag%rtpthvp (k) = real(host_mnts%rtpthvp (i,pverp-k+1), kind = core_rknd)
         core_diag%thlpthvp(k) = real(host_mnts%thlpthvp(i,pverp-k+1), kind = core_rknd)

         core_diag%rcm       (k) = real(host_mean%rcm(i,pverp-k+1), kind = core_rknd)
         core_diag%cloud_frac(k) = real( c2h% cloud_frac(i,pverp-k+1), kind = core_rknd)
      enddo

      !-----------------------------------------------
      ! Get input values for passive tracers
      !-----------------------------------------------
      core_prog%edsclr(:,:) = 0._core_rknd
 
      icnt=0
      do ixind=1,pcnst
         if (lq(ixind))  then
            icnt=icnt+1
            do k=1,pver
               core_prog%edsclr(k+1,icnt) = real(host_mean%q(i,pver-k+1,ixind), kind = core_rknd)
            enddo
            core_prog%edsclr(1,icnt) = core_prog%edsclr(2,icnt)
         end if
      enddo

      if (do_expldiff) then
        do k=1,pver
          core_prog%edsclr(k+1,icnt+1) = real(host_mean%thlm(i,pver-k+1), kind = core_rknd)
          core_prog%edsclr(k+1,icnt+2) = real(host_mean%rtm (i,pver-k+1), kind = core_rknd)
        enddo

        core_prog%edsclr(1,icnt+1) = core_prog%edsclr(2,icnt+1)
        core_prog%edsclr(1,icnt+2) = core_prog%edsclr(2,icnt+2)
      endif

      !-----------------------------------------------
      if (clubb_do_adv) then
        if (macmic_it .eq. 1) then

          core_prog%up2 = zt2zm_api( core_prog%up2 )
          core_prog%vp2 = zt2zm_api( core_prog%vp2 )
          core_prog%wp2 = zt2zm_api( core_prog%wp2 )

          core_prog%wpthlp  = zt2zm_api( core_prog%wpthlp  )
          core_prog%wprtp   = zt2zm_api( core_prog%wprtp   )
          core_prog%thlp2   = zt2zm_api( core_prog%thlp2   )
          core_prog%rtp2    = zt2zm_api( core_prog%rtp2    )
          core_prog%rtpthlp = zt2zm_api( core_prog%rtpthlp )

          do k=1,pverp
            core_prog%thlp2(k)=max(thl_tol**2,core_prog%thlp2(k))
            core_prog%rtp2 (k)=max( rt_tol**2,core_prog%rtp2 (k))

            core_prog%up2(k)=max(w_tol_sqd,core_prog%up2(k))
            core_prog%vp2(k)=max(w_tol_sqd,core_prog%vp2(k))
            core_prog%wp2(k)=max(w_tol_sqd,core_prog%wp2(k))
          enddo
        endif
      endif

      !---------------------------------------------------------------------------------------------
      ! In a restart run, assign values read back from restart file to PDF parameters on zm levels.
      ! This is necessary when ipdf_call_placement = 2.
      ! Note the vertical indices need to be flipped.
      !---------------------------------------------------------------------------------------------
      if ( is_first_restart_step() .and. ipdf_call_placement .eq. ipdf_post_advance_fields ) then

         do k=1,pverp
            pdf_params_zm% w_1       (k) = c2h% pdf_zm_w_1       (i,pverp-k+1)
            pdf_params_zm% w_2       (k) = c2h% pdf_zm_w_2       (i,pverp-k+1)
            pdf_params_zm% varnce_w_1(k) = c2h% pdf_zm_varnce_w_1(i,pverp-k+1)
            pdf_params_zm% varnce_w_2(k) = c2h% pdf_zm_varnce_w_2(i,pverp-k+1)
            pdf_params_zm% mixt_frac (k) = c2h% pdf_zm_mixt_frac (i,pverp-k+1)
         enddo

      end if

      !------------------------
      ! Cloud ice mixing ratio
      !------------------------
      do k=1,pver
         clubb_misc%rfrzm(k+1) = real( host_mean%q(i,pver-k+1,ixcldice), kind=core_rknd )
      enddo
      clubb_misc%rfrzm(1) = clubb_misc%rfrzm(2)

      !--------------------------------------------------------------------
      ! Get QRL for use in subroutine enhance_thlp2_for_cloud_top_cooling
      !--------------------------------------------------------------------
      if (do_rainturb) then
         do k=2,pverp
            clubb_misc% prer_evap(k) = real( host_mean%prer_evap(i,pverp-k+1), kind=core_rknd )
         enddo
         clubb_misc% prer_evap(1) = clubb_misc% prer_evap(2)
      end if

      if (do_cldcool) then
         do k=1,pver
            kkhost = pver-k+1
            clubb_misc% qrl_zt(k+1) = real( host_mean%qrl(i,kkhost)/(cpair*host_mean%pdel(i,kkhost)), kind=core_rknd )
         enddo
         clubb_misc% qrl_zt(1) = clubb_misc% qrl_zt(2)
      end if

      !--------------------------------------------------------------------
      ! Sandbox version of code to take into account meso organization
      !--------------------------------------------------------------------
      if (clubb_do_deep) then
         call varmu_cal( pver, host_mean%prer_evap(i,:), host_mean%pdel(i,:), & 
                         varmu, clubb_misc%varmu )
      endif

  end subroutine map_1column_of_info_from_host_to_clubb

  !=========================================================
  ! Code separated from clubb_tend_cam by Hui Wan, 2023-12
  !=========================================================
  subroutine varmu_cal( pver, prer_evap, pdel, varmu, varmu2 )

    use parameters_tunable, only: mu

    integer, intent(in) :: pver
    real(r8),intent(in) :: prer_evap(:)
    real(r8),intent(in) :: pdel(:)

    real(r8),       intent(out) :: varmu 
    real(core_rknd),intent(out) :: varmu2

    real(r8) :: orgparam, delpavg
    integer  :: k

    !------------------
    orgparam = 0._r8
    delpavg = 0._r8

    do k = 1, pver
       if (abs(prer_evap(k)) .gt. 0._r8) then
          orgparam = orgparam + (abs(prer_evap(k)) * 1000._r8 * 1000._r8 * 2._r8 ) * pdel(k)
          delpavg = delpavg + pdel(k)
       endif
    enddo

    if (delpavg .gt. 0._r8) then
      orgparam = orgparam/delpavg
    endif

    ! Now compute new entrainment rate based on organization

    varmu   = mu / (1._r8 + orgparam * 100._r8)
    varmu2  = real( varmu, kind=core_rknd )

  end subroutine varmu_cal
