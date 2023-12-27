module macrop_driver_with_clubb

  use perf_mod,      only: t_startf, t_stopf
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use cam_history,   only: outfld
  use ppgrid,        only: pver, pcols
  use constituents,  only: pcnst

  implicit none

  private
  public :: ice_supersat_adj_tend
  public :: deepcu_detrainment_tend
  public :: pblh_diag
  public :: gustiness

contains

  !---------------------------------------
  ! Ice Saturation Adjustment Computation 
  !---------------------------------------
  subroutine ice_supersat_adj_tend( state1, pbuf, hdtime, ptend_loc)

    use physics_types,  only: physics_state, physics_ptend, &
                              physics_state_copy, physics_ptend_init
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
    use constituents,   only: cnst_get_ind
    use physconst,      only: latvap, latice, cpair
    use ref_pres,       only: top_lev => trop_cloud_top_lev
    use macrop_driver,  only: ice_macro_tend

    type(physics_state),intent(in)    :: state1
    type(physics_buffer_desc),pointer :: pbuf(:)
    real(r8),intent(in)               :: hdtime

    type(physics_ptend),intent(out)   :: ptend_loc   ! Local tendency from process

    integer :: ncol, lchnk
    logical :: lq2(pcnst)
    integer :: ixcldice, ixnumice
    real(r8) :: latsub

    real(r8)  stend(pcols,pver)
    real(r8)  qvtend(pcols,pver)
    real(r8)  qitend(pcols,pver)
    real(r8)  initend(pcols,pver)

    real(r8), pointer, dimension(:,:) :: naai
    !---------------------------------------
     ncol  = state1%ncol
     lchnk = state1%lchnk

     call cnst_get_ind('CLDICE',ixcldice)
     call cnst_get_ind('NUMICE',ixnumice)

     lq2(:)  = .FALSE.
     lq2(1)  = .TRUE.
     lq2(ixcldice) = .TRUE.
     lq2(ixnumice) = .TRUE.

     latsub = latvap + latice

     call physics_ptend_init(ptend_loc, state1%psetcols, 'iceadj', ls=.true., lq=lq2 )

       stend(:ncol,:)=0._r8
      qvtend(:ncol,:)=0._r8
      qitend(:ncol,:)=0._r8
     initend(:ncol,:)=0._r8

     call pbuf_get_field(pbuf, pbuf_get_index('NAAI'), naai)

     call t_startf('ice_macro_tend')
     call ice_macro_tend(naai(:ncol,top_lev:pver),state1%t(:ncol,top_lev:pver), &
        state1%pmid(:ncol,top_lev:pver),state1%q(:ncol,top_lev:pver,1),state1%q(:ncol,top_lev:pver,ixcldice),&
        state1%q(:ncol,top_lev:pver,ixnumice),latsub,hdtime,&
        stend(:ncol,top_lev:pver),qvtend(:ncol,top_lev:pver),qitend(:ncol,top_lev:pver),&
        initend(:ncol,top_lev:pver))
     call t_stopf('ice_macro_tend')

     ! Update local copy of state with the tendencies

     ptend_loc%q(:ncol,top_lev:pver,1)        =  qvtend(:ncol,top_lev:pver)
     ptend_loc%q(:ncol,top_lev:pver,ixcldice) =  qitend(:ncol,top_lev:pver)
     ptend_loc%q(:ncol,top_lev:pver,ixnumice) = initend(:ncol,top_lev:pver)
     ptend_loc%s(:ncol,top_lev:pver)          =   stend(:ncol,top_lev:pver)

     ! Write output for tendencies
     call outfld( 'TTENDICE',  stend/cpair, pcols, lchnk )
     call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
     call outfld( 'QITENDICE', qitend, pcols, lchnk )
     call outfld( 'NITENDICE', initend, pcols, lchnk )

  end subroutine ice_supersat_adj_tend

  !---------------------------------------------------------------------------------
  ! Detrainment of convective condensate into the environment or stratiform cloud   
  !---------------------------------------------------------------------------------
  subroutine deepcu_detrainment_tend( state1, ixcldliq, ixcldice, ixnumliq, ixnumice, &
                                      dlf, clubb_tk1, clubb_tk2,    &
                                      clubb_liq_deep, clubb_liq_sh, &
                                      clubb_ice_deep, clubb_ice_sh, &
                                      ptend_loc, det_s, det_ice )

   use physconst,      only: cpair, gravit, latice
   use physics_types,  only: physics_state, physics_ptend, physics_ptend_init

   type(physics_state),intent(in) :: state1
   integer,  intent(in)  :: ixcldliq,ixcldice,ixnumliq,ixnumice
   real(r8), intent(in)  :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection [kg/ks/s]
   real(r8), intent(in)  :: clubb_tk1, clubb_tk2
   real(r8), intent(in)  :: clubb_liq_deep
   real(r8), intent(in)  :: clubb_liq_sh
   real(r8), intent(in)  :: clubb_ice_deep
   real(r8), intent(in)  :: clubb_ice_sh

   type(physics_ptend),intent(out) :: ptend_loc
   real(r8), intent(out) :: det_s  (pcols)    ! Integral of detrained static energy from ice
   real(r8), intent(out) :: det_ice(pcols)    ! Integral of detrained ice for energy check   

   integer  :: ncol, lchnk
   integer  :: i,k
   real(r8) :: invrs_gravit
   real(r8) :: dlf2(pcols,pver)  ! Detraining cld H20 from shallow convection    [kg/kg/day]
   logical  :: lqice(pcnst)
   real(r8) :: dum1

   ncol  = state1%ncol
   lchnk = state1%lchnk

   invrs_gravit = 1._r8 / gravit

   det_s(:)   = 0.0_r8
   det_ice(:) = 0.0_r8

   dlf2(:,:) = 0.0_r8  ! shallow convective detrainment rate is zero when CLUBB is used.

   lqice(:)        = .false.
   lqice(ixcldliq) = .true.
   lqice(ixcldice) = .true.
   lqice(ixnumliq) = .true.
   lqice(ixnumice) = .true.

   call physics_ptend_init(ptend_loc,state1%psetcols, 'clubb_det', ls=.true., lq=lqice)

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

  end subroutine deepcu_detrainment_tend

  !--------------------------------------------------------------------------
  ! Diagnose PBL height
  !--------------------------------------------------------------------------
  ! Hui Wan's note from 2023-12:
  !  This subroutine contains code lines separated from clubb_tend_cam.
  !  There seems to be a bug in the original code: the variable dz_g_bot
  !  here, which correspond to dz_g(pver) in the original code,
  !  should have a column dimension.
  !--------------------------------------------------------------------------
  subroutine pblh_diag( state1, cam_in, cloud_frac, dz_g_bot, use_sgv, pblh )

   use physconst,      only: gravit,zvir
   use physics_types,  only: physics_state
   use camsrfexch,     only: cam_in_t
   use constituents,   only: cnst_get_ind

   use pbl_utils,      only: calc_ustar, calc_obklen
   use hb_diff,        only: pblintd 

   type(physics_state),intent(in)  :: state1
   type(cam_in_t),     intent(in)  :: cam_in
   real(r8),           intent(in)  :: cloud_frac(:,:)
   real(r8),           intent(in)  :: dz_g_bot
   logical,            intent(in)  :: use_sgv
   real(r8),           intent(out) :: pblh(:)

   real(r8) :: th(pcols,pver)   ! potential temperature                         [K]
   real(r8) :: thv(pcols,pver)  ! virtual potential temperature                 [K]
   real(r8) :: ustar2(pcols)    ! Surface stress for PBL height                 [m2/s2]

   real(r8) :: invrs_gravit
   real(r8) :: rrho             ! Inverse of air density                        [1/kg/m^3]

   real(r8) :: kinheat(pcols)   ! Kinematic Surface heat flux                   [K m/s]
   real(r8) ::  kinwat(pcols)   ! Kinematic water vapor flux                    [m/s]
   real(r8) ::    kbfs(pcols)   ! Kinematic Surface heat flux                   [K m/s]
   real(r8) ::  obklen(pcols)   ! Obukov length                                 [m]
   real(r8) ::  dummy2(pcols)   ! dummy variable
   real(r8) ::  dummy3(pcols)   ! dummy variable

   integer :: ixcldliq, ixq
   integer :: i,k, ncol
   
   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)

   ncol = state1%ncol
 
   ! --------------------------------------------------------------------------------- !
   !  DIAGNOSE THE PBL DEPTH                                                           !
   !  this is needed for aerosol code                                                  !
   ! --------------------------------------------------------------------------------- !

   do i=1,ncol
      do k=1,pver
         th(i,k) = state1%t(i,k)*state1%exner(i,k)
         if (use_sgv) then
           thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq) &
                    - state1%q(i,k,ixcldliq))  !PMA corrects thv formula
         else
           thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq))
         end if
      enddo
   enddo

   ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)

   invrs_gravit = 1._r8 / gravit

   do i=1,ncol
      rrho = invrs_gravit*(state1%pdel(i,pver)/dz_g_bot)
      call calc_ustar( state1%t(i,pver), state1%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                       rrho, ustar2(i) )
      call calc_obklen( th(i,pver), thv(i,pver), cam_in%cflx(i,1), cam_in%shf(i), rrho, ustar2(i), &
                        kinheat(i), kinwat(i), kbfs(i), obklen(i) )
   enddo

   dummy2(:) = 0._r8
   dummy3(:) = 0._r8

   where (kbfs .eq. -0.0_r8) kbfs = 0.0_r8

   !  Compute PBL depth according to Holtslag-Boville Scheme
   call pblintd(ncol, thv, state1%zm, state1%u, state1%v, &
                ustar2, obklen, kbfs, pblh, dummy2, &
                state1%zi, cloud_frac(:,1:pver), 1._r8-cam_in%landfrac, dummy3)

  end subroutine pblh_diag
  !------------------------

  !----------------------------------------------------------------------------------
  ! The gustiness subroutine contains code lines separated from clubb_tend_cam.
  ! Note in the original subroutine: !PMA adds gustiness and tpert
  !----------------------------------------------------------------------------------
  subroutine gustiness( use_sgv, state1, cam_in, up2b, vp2b, prec_dp, snow_dp, pblh, thlp2, &! in
                        vmag_gust, tpert )

    use physics_types,  only: physics_state
    use camsrfexch,     only: cam_in_t
    use constituents,   only: cnst_get_ind
    use physconst,      only: latvap, cpair

    logical,            intent(in)  :: use_sgv
    type(physics_state),intent(in)  :: state1
    type(cam_in_t),     intent(in)  :: cam_in
    real(r8),           intent(in)  :: up2b(:), vp2b(:)
    real(r8),           intent(in)  :: thlp2(:,:)
    real(r8),           intent(in)  :: prec_dp(:), snow_dp(:)
    real(r8),           intent(in)  :: pblh(:)
    real(r8),           intent(out) :: vmag_gust(:)
    real(r8),           intent(out) :: tpert(:)

    integer  :: ktopi(pcols)
    real(r8) :: umb(pcols), vmb(pcols)
    real(r8) :: prec_gust(pcols)
    real(r8) :: gust_fac(pcols)
    real(r8) :: vmag_gust_dp(pcols),vmag_gust_cl(pcols)
    real(r8) :: vmag(pcols)
    integer :: ixcldliq
    integer :: ncol, lchnk, i,k

    real(r8),parameter :: gust_facl = 1.2_r8 !gust fac for land
    real(r8),parameter :: gust_faco = 0.9_r8 !gust fac for ocean
    real(r8),parameter :: gust_facc = 1.5_r8 !gust fac for clubb

    ncol  = state1%ncol
    lchnk = state1%lchnk

    call cnst_get_ind('CLDLIQ',ixcldliq)

    vmag_gust(:)    = 0._r8

    vmag_gust_dp(:) = 0._r8
    vmag_gust_cl(:) = 0._r8
    ktopi(:)        = pver

    if (use_sgv) then
       do i=1,ncol
          !up2b(i)          = up2(i,pver)
          !vp2b(i)          = vp2(i,pver)
           umb(i)           = state1%u(i,pver)
           vmb(i)           = state1%v(i,pver)
           prec_gust(i)     = max(0._r8,prec_dp(i)-snow_dp(i))*1.e3_r8
           if (cam_in%landfrac(i).gt.0.95_r8) then
             gust_fac(i)   = gust_facl
           else
             gust_fac(i)   = gust_faco
           endif
           vmag(i)         = max(1.e-5_r8,sqrt( umb(i)**2._r8 + vmb(i)**2._r8))
           vmag_gust_dp(i) = ugust(min(prec_gust(i),6.94444e-4_r8),gust_fac(i)) ! Limit for the ZM gustiness equation set in Redelsperger et al. (2000)
           vmag_gust_dp(i) = max(0._r8, vmag_gust_dp(i) )!/ vmag(i))
           vmag_gust_cl(i) = gust_facc*(sqrt(max(0._r8,up2b(i)+vp2b(i))+vmag(i)**2._r8)-vmag(i))
           vmag_gust_cl(i) = max(0._r8, vmag_gust_cl(i) )!/ vmag(i))
           vmag_gust(i)    = vmag_gust_cl(i) + vmag_gust_dp(i)

          do k=1,pver
             if (state1%zi(i,k)>pblh(i).and.state1%zi(i,k+1)<=pblh(i)) then
                ktopi(i) = k
                exit
             end if
          end do

          tpert(i) = min(2._r8,(sqrt(thlp2(i,ktopi(i)))+(latvap/cpair)*state1%q(i,ktopi(i),ixcldliq)) &
                    /max(state1%exner(i,ktopi(i)),1.e-3_r8)) !proxy for tpert
       end do
    end if

   call outfld('VMAGDP', vmag_gust_dp, pcols, lchnk)
   call outfld('VMAGCL', vmag_gust_cl, pcols, lchnk)
   call outfld('VMAGGUST', vmag_gust, pcols, lchnk)
   call outfld('TPERTBLT', tpert, pcols, lchnk)

   contains

     ! ZM gustiness equation below from Redelsperger et al. (2000)
     ! numbers are coefficients of the empirical equation.
     function ugust(gprec,gfac)

       real(r8) :: ugust  ! function: gustiness as a function of convective rainfall
       real(r8) :: gfac
       real(r8) :: gprec 
           
       ugust = gfac*log(1._R8+57801.6_R8*gprec-3.55332096e7_R8*(gprec**2.0_R8))

     end function

  end subroutine gustiness
  !------------------------

end module macrop_driver_with_clubb
