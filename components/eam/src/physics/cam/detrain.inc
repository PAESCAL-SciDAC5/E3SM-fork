
   ! --------------------------------------------------------------------------------- !
   !  COMPUTE THE ICE CLOUD DETRAINMENT                                                !
   !  Detrainment of convective condensate into the environment or stratiform cloud    !
   ! --------------------------------------------------------------------------------- !
   invrs_gravit = 1._r8 / gravit

   det_s(:)   = 0.0_r8
   det_ice(:) = 0.0_r8

   !  Initialize the shallow convective detrainment rate, will always be zero
   dlf2(:,:) = 0.0_r8

   call cnst_get_ind('NUMLIQ',ixnumliq)
   call cnst_get_ind('NUMICE',ixnumice)

   lqice(:)        = .false.
   lqice(ixcldliq) = .true.
   lqice(ixcldice) = .true.
   lqice(ixnumliq) = .true.
   lqice(ixnumice) = .true.

   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb_det', ls=.true., lq=lqice)

   call t_startf('ice_cloud_detrain_diag')
   do k=1,pver
      do i=1,ncol
         if( state1%t(i,k) > clubb_tk1 ) then
            dum1 = 0.0_r8
         elseif ( state1%t(i,k) < clubb_tk2 ) then
            dum1 = 1.0_r8
         else
            !Note: Denominator is changed from 30.0_r8 to (clubb_tk1 - clubb_tk2),
            !(clubb_tk1 - clubb_tk2) is also 30.0 but it introduced a non-bfb change
            dum1 = ( clubb_tk1 - state1%t(i,k) ) /(clubb_tk1 - clubb_tk2)
         endif

         ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
         ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
         ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) &
                                     / (4._r8*3.14_r8* clubb_liq_deep**3*997._r8) + & ! Deep    Convection
                                     3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) &
                                     / (4._r8*3.14_r8*clubb_liq_sh**3*997._r8)     ! Shallow Convection
         ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) &
                                     / (4._r8*3.14_r8*clubb_ice_deep**3*500._r8) + & ! Deep    Convection
                                     3._r8 * (                         dlf2(i,k)    *  dum1 ) &
                                     / (4._r8*3.14_r8*clubb_ice_sh**3*500._r8)     ! Shallow Convection
         ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice

         ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
         !   track of the integrals of ice and static energy that is effected from conversion to ice
         !   so that the energy checker doesn't complain.
         det_s(i)                  = det_s(i) + ptend_loc%s(i,k)*state1%pdel(i,k)*invrs_gravit
         det_ice(i)                = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state1%pdel(i,k)*invrs_gravit

      enddo
   enddo

   det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water
   call t_stopf('ice_cloud_detrain_diag')

   call outfld( 'DPDLFLIQ', ptend_loc%q(:,:,ixcldliq), pcols, lchnk)
   call outfld( 'DPDLFICE', ptend_loc%q(:,:,ixcldice), pcols, lchnk)
   call outfld( 'DPDLFT',   ptend_loc%s(:,:)/cpair, pcols, lchnk)
