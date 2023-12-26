
 !subroutine gather_moment_fields( pbuf, itim_old, clubb_do_adv, apply_const, host_mnts )
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
