
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

   ncol           = state1%ncol
   host_mean%ncol = state1%ncol

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

   host_mean%lat   => state1%lat
   host_mean%lon   => state1%lon
   
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
  end subroutine gather_host_fields

  !===============================================================================
  !===============================================================================
  subroutine setup_clubb_to_host_fields( pbuf, c2h )

   use physics_buffer,  only: physics_buffer_desc, pbuf_get_field
   use clubb_intr_types,only: clubb_to_host_t

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(clubb_to_host_t)              :: c2h

   call pbuf_get_field( pbuf, kvm_idx, c2h%khzt )
   call pbuf_get_field( pbuf, kvh_idx, c2h%khzm )

  end subroutine setup_clubb_to_host_fields
