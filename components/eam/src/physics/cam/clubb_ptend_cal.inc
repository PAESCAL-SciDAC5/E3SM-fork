  subroutine derive_ptend( host_mean, host_mnts, state1, clubb_s, edsclr_out, hdtime,    &! in
                           ixq, ixcldliq, clubb_do_adv, macmic_it, cld_macmic_num_steps, &! in
                           ptend_loc )! inout
  
    use physics_types, only: physics_state, physics_ptend
    use clubb_intr_types, only: clubb_mean_2d_t, clubb_mnts_2d_t

    type(clubb_mean_2d_t),intent(in) :: host_mean
    type(clubb_mnts_2d_t),intent(in) :: host_mnts
    type(physics_state),intent(in) :: state1

    real(r8),intent(in) :: clubb_s   (:,:)
    real(r8),intent(in) :: edsclr_out(:,:,:)

    real(r8),intent(in) :: hdtime
    integer, intent(in) :: ixq, ixcldliq
    logical, intent(in) :: clubb_do_adv
    integer, intent(in) :: macmic_it, cld_macmic_num_steps

    type(physics_ptend),intent(inout) :: ptend_loc

    real(r8) :: invrs_hdtime
    integer  :: i, k
    integer  :: icnt, ixind

    !---------------------------------------------------------------------------------
    !  Now compute the tendencies of CLUBB to CAM, note that pverp is the ghost point
    !  for all variables and therefore is never called in this loop
    !---------------------------------------------------------------------------------
    invrs_hdtime = 1._r8 / hdtime

    do i=1,state1%ncol
      do k=1,pver

         ptend_loc%u(i,k)   = ( host_mean%um(i,k) - state1%u(i,k) )*invrs_hdtime    ! east-west wind
         ptend_loc%v(i,k)   = ( host_mean%vm(i,k) - state1%v(i,k) )*invrs_hdtime    ! north-south wind
         ptend_loc%q(i,k,ixq) = (host_mean%rtm(i,k)-host_mean%rcm(i,k)-state1%q(i,k,ixq))*invrs_hdtime  ! water vapor
         ptend_loc%q(i,k,ixcldliq) = (host_mean%rcm(i,k)-state1%q(i,k,ixcldliq))*invrs_hdtime ! Tendency of liquid water
         ptend_loc%s(i,k)   = (clubb_s(i,k)-state1%s(i,k))*invrs_hdtime             ! Tendency of static energy

         if (clubb_do_adv) then
            if (macmic_it .eq. cld_macmic_num_steps) then

               ! Here add a constant to moments which can be either positive or
               !  negative.  This is to prevent clipping when dynamics tries to
               !  make all constituents positive
               host_mnts%     wp3(i,k) = host_mnts%     wp3(i,k) +     wp3_const
               host_mnts% rtpthlp(i,k) = host_mnts% rtpthlp(i,k) + rtpthlp_const
               host_mnts%  wpthlp(i,k) = host_mnts%  wpthlp(i,k) +  wpthlp_const
               host_mnts%   wprtp(i,k) = host_mnts%   wprtp(i,k) +   wprtp_const

               ptend_loc%q(i,k,ixthlp2)   = (host_mnts%  thlp2(i,k)-state1%q(i,k,ixthlp2  ))*invrs_hdtime   ! THLP Variance
               ptend_loc%q(i,k,ixrtp2)    = (host_mnts%   rtp2(i,k)-state1%q(i,k,ixrtp2   ))*invrs_hdtime   ! RTP Variance
               ptend_loc%q(i,k,ixrtpthlp) = (host_mnts%rtpthlp(i,k)-state1%q(i,k,ixrtpthlp))*invrs_hdtime   ! RTP THLP covariance
               ptend_loc%q(i,k,ixwpthlp)  = (host_mnts% wpthlp(i,k)-state1%q(i,k,ixwpthlp ))*invrs_hdtime   ! WPTHLP
               ptend_loc%q(i,k,ixwprtp)   = (host_mnts%  wprtp(i,k)-state1%q(i,k,ixwprtp  ))*invrs_hdtime   ! WPRTP
               ptend_loc%q(i,k,ixwp2)     = (host_mnts%    wp2(i,k)-state1%q(i,k,ixwp2    ))*invrs_hdtime   ! WP2
               ptend_loc%q(i,k,ixwp3)     = (host_mnts%    wp3(i,k)-state1%q(i,k,ixwp3    ))*invrs_hdtime   ! WP3
               ptend_loc%q(i,k,ixup2)     = (host_mnts%    up2(i,k)-state1%q(i,k,ixup2    ))*invrs_hdtime   ! UP2
               ptend_loc%q(i,k,ixvp2)     = (host_mnts%    vp2(i,k)-state1%q(i,k,ixvp2    ))*invrs_hdtime   ! VP2
            else
               ptend_loc%q(i,k,ixthlp2)=0._r8
               ptend_loc%q(i,k,ixrtp2)=0._r8
               ptend_loc%q(i,k,ixrtpthlp)=0._r8
               ptend_loc%q(i,k,ixwpthlp)=0._r8
               ptend_loc%q(i,k,ixwprtp)=0._r8
               ptend_loc%q(i,k,ixwp2)=0._r8
               ptend_loc%q(i,k,ixwp3)=0._r8
               ptend_loc%q(i,k,ixup2)=0._r8
               ptend_loc%q(i,k,ixvp2)=0._r8
            endif

         endif

         !  Apply tendencies to ice mixing ratio, liquid and ice number, and aerosol constituents.
         !  Loading up this array doesn't mean the tendencies are applied.
         ! edsclr_out is compressed with just the constituents being used, ptend and state are not compressed

         icnt=0
         do ixind=1,pcnst
            if (lq(ixind)) then
               icnt=icnt+1
               if ((ixind /= ixq)       .and. (ixind /= ixcldliq) .and.&
                   (ixind /= ixthlp2)   .and. (ixind /= ixrtp2)   .and.&
                   (ixind /= ixrtpthlp) .and. (ixind /= ixwpthlp) .and.&
                   (ixind /= ixwprtp)   .and. (ixind /= ixwp2)    .and.&
                   (ixind /= ixwp3)     .and. (ixind /= ixup2)    .and. (ixind /= ixvp2) ) then
                       ptend_loc%q(i,k,ixind) = (edsclr_out(i,k,icnt)-state1%q(i,k,ixind))*invrs_hdtime ! transported constituents
               end if
            end if
         enddo

      enddo
    end do

    !--------------------------------------------------------------
    ! Only for k = pverp:
    ! Add constant to ghost point so that output is not corrupted.
    !--------------------------------------------------------------
    if (clubb_do_adv) then
       if (macmic_it .eq. cld_macmic_num_steps) then
          host_mnts%     wp3(:,pverp) = host_mnts%     wp3(:,pverp) +     wp3_const
          host_mnts% rtpthlp(:,pverp) = host_mnts% rtpthlp(:,pverp) + rtpthlp_const
          host_mnts%  wpthlp(:,pverp) = host_mnts%  wpthlp(:,pverp) +  wpthlp_const
          host_mnts%   wprtp(:,pverp) = host_mnts%   wprtp(:,pverp) +   wprtp_const
       endif
    endif

  end subroutine derive_ptend
