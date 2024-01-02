
  subroutine set_clubb_forcing_and_sfc( icol, pverp, nz, host_mean, core_auxil, &! in
                                        do_tms, ksrftms,                        &! in
                                        core_forcing, core_sfc, fcoriolis       )! inout or out

    use physconst,       only: cpair
    use clubb_intr_types,only: clubb_mean_2d_t
    use clubb_intr_core_types, only: core_auxil_t, core_forcing_t, core_sfc_t

    integer,intent(in) :: icol
    integer,intent(in) :: pverp
    integer,intent(in) :: nz

    type(clubb_mean_2d_t), intent(in) :: host_mean
    type(core_auxil_t),    intent(in) :: core_auxil
    logical,               intent(in) :: do_tms
    real(r8),              intent(in) :: ksrftms(:)
    type(core_forcing_t),  intent(inout) :: core_forcing
    type(core_sfc_t),      intent(inout) :: core_sfc
    real(core_rknd),       intent(out)   :: fcoriolis

    real(core_rknd) :: dum_core_rknd

      !-----------

      if (pverp/=nz) call endrun('pverp/=nz not supported by current code')

      !==================================================================================================
      ! Non-zero surface fluxes
      !  Hui Wan's note, 2023-12: in an AMR setting, we need to be careful to use the correct rho values
      !==================================================================================================
      !  Surface fluxes provided by host model
      !-----------------------------------------
      core_sfc% wpthlp = real(host_mean%shf (icol),   kind = core_rknd)/ &
                        (real(cpair, kind = core_rknd)* core_auxil%rho_ds_zm(1))                 ! Sensible heat flux
      core_sfc% wprtp  = real(host_mean%cflx(icol,1), kind = core_rknd)/ core_auxil%rho_ds_zm(1) ! Latent heat flux
      core_sfc% upwp   = real(host_mean%wsx (icol),   kind = core_rknd)/ core_auxil%rho_ds_zm(1) ! Surface meridional momentum flux
      core_sfc% vpwp   = real(host_mean%wsy (icol),   kind = core_rknd)/ core_auxil%rho_ds_zm(1) ! Surface zonal momentum flux

      !-----------
      ! Apply TMS
      !-----------
       if ( do_tms ) then
         dum_core_rknd = real((ksrftms(icol)*host_mean%um(icol,pver)), kind = core_rknd)
         core_sfc%upwp = core_sfc%upwp-(dum_core_rknd/core_auxil%rho_ds_zm(1))
         dum_core_rknd = real((ksrftms(icol)*host_mean%vm(icol,pver)), kind = core_rknd)
         core_sfc%vpwp = core_sfc%vpwp-(dum_core_rknd/core_auxil%rho_ds_zm(1))
       endif

      !  Define surface sources for transported variables for diffusion, will
      !  be zero as these tendencies are done in clubb_surface
      !  ---------------
      !  Hui Wan's question, 2023-12:
      !   Is the statement above correct?
      !    - If edsclr contains "eddy passive scalars" such as chemical gases,
      !      then in principle we could have eddy fluxes of such gases.
      !    - If do_expldiff = .t., then thlm and rtm will be transported as "eddy passive tracers",
      !      then shouldn't they have surface fluxes?
      !  ---------------
      core_sfc% wpedsclrp(1:edsclr_dim) = 0._core_rknd

   !======================================================================
   ! Forcing terms that are zero in the current implementation
   !======================================================================
   !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
   !  forcings already applied through CAM

   core_forcing% um  (1:nz) = 0._core_rknd
   core_forcing% vm  (1:nz) = 0._core_rknd

   core_forcing% thlm(1:nz) = 0._core_rknd
   core_forcing% rtm (1:nz) = 0._core_rknd

   core_forcing% rtp2   (1:nz) = 0._core_rknd
   core_forcing% thlp2  (1:nz) = 0._core_rknd
   core_forcing% rtpthlp(1:nz) = 0._core_rknd
   core_forcing% wprtp  (1:nz) = 0._core_rknd
   core_forcing% wpthlp (1:nz) = 0._core_rknd

   core_forcing% edsclr(:,:) = 0._core_rknd

   !-------------------------------------------------------------------
   !  Determine Coriolis force at given latitude.  This is never used
   !  when CLUBB is implemented in a host model, therefore just set to zero.
   !-------------------------------------------------------------------
   fcoriolis = 0._core_rknd

  end subroutine set_clubb_forcing_and_sfc


  !============================================================
  subroutine set_clubb_sfc_for_scm( scm_clubb_iop_name, icol, pver, host_mean, zt_bot, theta0, core_sfc )
 
    use physconst,             only: gravit 
    use clubb_intr_types,      only: clubb_mean_2d_t
    use clubb_intr_core_types, only: core_sfc_t

    character(len=*),     intent(in)    :: scm_clubb_iop_name 
    integer,              intent(in)    :: icol
    integer,              intent(in)    :: pver
    type(clubb_mean_2d_t),intent(in)    :: host_mean
    real(core_rknd),      intent(in)    :: zt_bot
    real(core_rknd),      intent(in)    :: theta0
    type(core_sfc_t),     intent(inout) :: core_sfc

    real(r8) :: zo
    real(r8) :: ubar, ustar
    real(r8) :: dum1, bflx22

    !----------------------------------------------
    !  Initialize zo if variable ustar is used

    if (host_mean%landfrac(icol) .ge. 0.5_r8) then
       zo = 0.035_r8
    else
       zo = 0.0001_r8
    endif

    !  Compute surface wind (ubar)
    ubar = sqrt(host_mean%um(icol,pver)**2+host_mean%vm(icol,pver)**2)
    if (ubar .lt. 0.25_r8) ubar = 0.25_r8

    !  Below denotes case specifics for surface momentum
    !  and thermodynamic fluxes, depending on the case

    !  Define ustar (based on case, if not variable)
    ustar = 0.25_r8   ! Initialize ustar in case no case

    if(trim(scm_clubb_iop_name) .eq. 'BOMEX_5day') then
       ustar = 0.28_r8
    endif

    if(trim(scm_clubb_iop_name) .eq. 'ATEX_48hr') then
       ustar = 0.30_r8
    endif

    if(trim(scm_clubb_iop_name) .eq. 'RICO_3day') then
       ustar = 0.28_r8
    endif

    if(trim(scm_clubb_iop_name) .eq. 'arm97' .or. trim(scm_clubb_iop_name) .eq. 'gate' .or. &
       trim(scm_clubb_iop_name) .eq. 'toga' .or. trim(scm_clubb_iop_name) .eq. 'mpace' .or. &
       trim(scm_clubb_iop_name) .eq. 'ARM_CC') then

         dum1   = real(zt_bot, kind = r8)
         bflx22 = (gravit/real(theta0, kind = r8))*real(core_sfc%wpthlp, kind = r8)
         ustar  = diag_ustar(dum1,bflx22,ubar,zo)
    endif

    !  Compute the surface momentum fluxes, if this is a SCAM simulation
    core_sfc% upwp = -real( (host_mean%um(icol,pver)*ustar**2/ubar), kind=core_rknd )
    core_sfc% vpwp = -real( (host_mean%vm(icol,pver)*ustar**2/ubar), kind=core_rknd )

  end subroutine set_clubb_sfc_for_scm
