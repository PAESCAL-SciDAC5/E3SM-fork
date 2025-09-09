module misc_diagnostics

use shr_kind_mod,   only: r8 => shr_kind_r8

implicit none
public

contains


!------------------------------------------------
! saturation specific humidity wrt ice.
! The calculation in this subroutine follows 
! what is used in the model for 
!  - history output
!  - ice nucleation parameterization
! 
subroutine qsat_ice( ncol, pver, tair, pair, qsati )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)

  real(r8),intent(out) :: qsati(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  qsati = qsatw*esi/esl

end subroutine qsat_ice
!----------------------------------

!-----------------------------------------------------------
! supersaturation wrt water (liquid) given as mixing ratio
!
subroutine supersat_q_water( ncol, pver, tair, pair, qv, qssatw )

  use wv_saturation, only: qsat_water

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: qssatw(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver) !not used

  call qsat_water( tair, pair, esl, qsatw )
  qssatw = qv-qsatw

end subroutine supersat_q_water
!----------------------------------

!-----------------------------------------------------------
! supersaturation wrt ice given as mixing ratio
!
subroutine supersat_q_ice( ncol, pver, tair, pair, qv, qssati )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: qssati(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  qssati = qv - qsatw*esi/esl


end subroutine supersat_q_ice
!----------------------------------

subroutine relhum_water_percent( ncol, pver, tair, pair, qv,  rhw_percent )

  use wv_saturation, only: qsat_water

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: rhw_percent(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver) !not used

  call qsat_water( tair, pair, esl, qsatw )
  rhw_percent = qv/qsatw*100._r8

end subroutine relhum_water_percent


!----------------------------------
subroutine relhum_ice_percent( ncol, pver, tair, pair, qv,  rhi_percent )

  use wv_saturation, only: qsat_water, svp_ice

  integer, intent(in)  :: ncol, pver

  real(r8),intent(in)  :: tair(ncol,pver)
  real(r8),intent(in)  :: pair(ncol,pver)
  real(r8),intent(in)  ::   qv(ncol,pver)

  real(r8),intent(out) :: rhi_percent(ncol,pver)

  real(r8) :: qsatw(ncol,pver)
  real(r8) ::   esl(ncol,pver)
  real(r8) ::   esi(ncol,pver)

  integer :: i,k

  call qsat_water( tair, pair, esl, qsatw )

  do i=1,ncol
  do k=1,pver
     esi(i,k)=svp_ice(tair(i,k))
  end do
  end do

  rhi_percent = 100._r8* qv/qsatw* esl/esi

end subroutine relhum_ice_percent

subroutine compute_cape( state, pbuf, pcols, pver, cape )
!----------------------------------------------------------------------
! Purpose: compute CAPE (convecitve available potential energy)
!          using subroutine buoyan_dilute from the ZM deep convection
!          parameterization (module file zm_conv.F90)
! History: first version by Hui Wan, 2021-05
!----------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: cpair, gravit, rair, latvap
  use zm_conv,        only: buoyan_dilute, limcnv

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pver
  integer,                  intent(in) :: pcols
  real(r8),                intent(out) :: cape(pcols)

  ! local variables used for providing input to subroutine buoyan_dilute

  real(r8) :: pmid_in_hPa(pcols,pver)
  real(r8) :: pint_in_hPa(pcols,pver+1)

  real(r8),pointer ::    qv(:,:)
  real(r8),pointer ::  temp(:,:)

  real(r8) ::   zs(pcols)
  real(r8) :: pblt(pcols)
  real(r8) :: zmid_above_sealevel(pcols,pver)

  real(r8),pointer :: tpert(:), pblh(:)

  integer :: idx, kk, lchnk, ncol, msg

  ! variables returned by buoyan_dilute but not needed here

  real(r8) ::   ztp(pcols,pver) ! parcel temperatures.
  real(r8) :: zqstp(pcols,pver) ! grid slice of parcel temp. saturation mixing ratio.
  real(r8) ::   ztl(pcols)      ! parcel temperature at lcl.
  integer  ::  zlcl(pcols)      ! base level index of deep cumulus convection.
  integer  ::  zlel(pcols)      ! index of highest theoretical convective plume.
  integer  ::  zlon(pcols)      ! index of onset level for deep convection.
  integer  :: zmaxi(pcols)      ! index of level with largest moist static energy.

  !----------------------------------------------------------------------- 
  ncol  = state%ncol
  lchnk = state%lchnk

  msg = limcnv - 1  ! limcnv is the top interface level limit for convection

  idx = pbuf_get_index('tpert') ; call pbuf_get_field( pbuf, idx, tpert )

  pmid_in_hPa(1:ncol,:) = state%pmid(1:ncol,:) * 0.01_r8
  pint_in_hPa(1:ncol,:) = state%pint(1:ncol,:) * 0.01_r8

  qv   => state%q(:,:,1)
  temp => state%t

  ! Surface elevation (m) is needed to calculate height above sea level (m) 
  ! Note that zm (and zi) stored in state are height above surface. 
  ! The layer midpoint height provided to buoyan_dilute is height above sea level. 

  zs(1:ncol) = state%phis(1:ncol)/gravit

  !------------------------------------------
  ! Height above sea level at layer midpoints
  !------------------------------------------
  do kk = 1,pver
     zmid_above_sealevel(1:ncol,kk) = state%zm(1:ncol,kk)+zs(1:ncol)
  end do

  !--------------------------
  ! layer index for PBL top
  !--------------------------
  idx = pbuf_get_index('pblh')  ; call pbuf_get_field( pbuf, idx, pblh )

  pblt(:) = pver
  do kk = pver-1, msg+1, -1
     where( abs(zmid_above_sealevel(:ncol,kk)-zs(:ncol)-pblh(:ncol))   &
            < (state%zi(:ncol,kk)-state%zi(:ncol,kk+1))*0.5_r8       ) & 
     pblt(:ncol) = kk
  end do

  !----------------
  ! Calculate CAPE
  !----------------
  call buoyan_dilute(lchnk ,ncol, qv, temp,            &! in
                     pmid_in_hPa, zmid_above_sealevel, &! in
                     pint_in_hPa,                      &! in
                     ztp, zqstp, ztl,                  &! out
                     latvap, cape, pblt,               &! in, out, in
                     zlcl, zlel, zlon, zmaxi,          &! out
                     rair, gravit, cpair, msg, tpert,  &! in
                     .true.                            )! in (.t. = std CAPE calculation)

 end subroutine compute_cape
!---------------------------

subroutine ncic_diag( state, pbuf, pcols, pver, ncic )
!----------------------------------------------------------------------
! Purpose: diagnose in-cloud droplet number concentration (unit: 1/cc)
! History: first version by Hui Wan, 2022-06
!----------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: rair
  use micro_mg_utils, only: mincld
  use constituents,   only: cnst_get_ind

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pcols
  integer,                  intent(in) :: pver
  real(r8),                intent(out) :: ncic(pcols,pver)

  integer :: ncol, ixnumliq

  integer :: idx
  real(r8), pointer :: liqcldf(:,:)  ! liquid cloud fraction 

  real(r8) :: lcldm(pcols,pver)      ! liquid cloud fraction, clipped to mincld
  real(r8) ::   rho(pcols,pver)      ! air density 

  !-----------------------
  ncol = state%ncol

  ! Assume the liquid cloud fraction (liqcldf) equals the total stratiform cloud fraction (ast).
  ! This assumption follows the subroutine micro_mg_cam:micro_mg_cam_tend.

  idx = pbuf_get_index('AST') ; call pbuf_get_field(pbuf, idx, liqcldf)

  ! Clip the liquid cloud fraction to avoid division by zero.
  ! This treatment follows the subroutine micro_mg2_0:micro_mg_tend.

  lcldm(:ncol,:) = max( liqcldf(:ncol,:), mincld )

  ! Calculate air density, to be used in unit conversion from 1/kg to 1/m3
  ! following the subroutine micro_mg2_0:micro_mg_tend.

  rho(:ncol,:) = state%pmid(:ncol,:) / ( rair * state%t(:ncol,:) )

  ! Now, calculate the in-cloud droplet number concentration

  call cnst_get_ind( 'NUMLIQ', ixnumliq )
  ncic(:ncol,:) = state%q(:ncol,:,ixnumliq) * rho(:ncol,:) / lcldm(:ncol,:)

end subroutine ncic_diag


subroutine qcic_diag( state, pbuf, pcols, pver, qcic )
!----------------------------------------------------------------------
! Purpose: diagnose in-cloud droplet mass concentration (unit: g/cc)
! History: first version by Hui Wan, 2022-06
!----------------------------------------------------------------------

  use physics_types,  only: physics_state
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use physconst,      only: rair
  use micro_mg_utils, only: mincld
  use constituents,   only: cnst_get_ind

  type(physics_state),intent(in),target:: state
  type(physics_buffer_desc),pointer    :: pbuf(:)
  integer,                  intent(in) :: pcols
  integer,                  intent(in) :: pver
  real(r8),                intent(out) :: qcic(pcols,pver)

  integer :: ncol, ixcldliq

  integer :: idx
  real(r8), pointer :: liqcldf(:,:)  ! liquid cloud fraction 

  real(r8) :: lcldm(pcols,pver)      ! liquid cloud fraction, clipped to mincld
  real(r8) ::   rho(pcols,pver)      ! air density 

  !-----------------------
  ncol = state%ncol

  ! Assume the liquid cloud fraction (liqcldf) equals the total stratiform cloud fraction (ast).
  ! This assumption follows the subroutine micro_mg_cam:micro_mg_cam_tend.

  idx = pbuf_get_index('AST') ; call pbuf_get_field(pbuf, idx, liqcldf)

  ! Clip the liquid cloud fraction to avoid division by zero.
  ! This treatment follows the subroutine micro_mg2_0:micro_mg_tend.

  lcldm(:ncol,:) = max( liqcldf(:ncol,:), mincld )

  ! Calculate air density, to be used in unit conversion from 1/kg to 1/m3
  ! following the subroutine micro_mg2_0:micro_mg_tend.

  rho(:ncol,:) = state%pmid(:ncol,:) / ( rair * state%t(:ncol,:) )

  ! Now, calculate the in-cloud droplet mass concentration

  call cnst_get_ind( 'CLDLIQ', ixcldliq )
  qcic(:ncol,:) = 1000._r8* state%q(:ncol,:,ixcldliq) * rho(:ncol,:) / lcldm(:ncol,:)

end subroutine qcic_diag

subroutine tmp_numliq_update_after_activation( state_in, pbuf, dt, state_out )

  use physics_types,  only: physics_state, physics_state_copy, physics_ptend, physics_ptend_init
  use physics_update_mod, only: physics_update
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use constituents,   only: pcnst, cnst_get_ind

  type(physics_state),intent(in)  :: state_in
  type(physics_buffer_desc),pointer :: pbuf(:)
  real(r8) :: dt
  type(physics_state),intent(out) :: state_out

  ! local variables

  integer :: ncol, idx
  real(r8), pointer :: npccn(:,:)
  logical :: lq(pcnst) 
  integer :: ixnumliq
  type(physics_ptend) :: ptend

  !--------------------------
  ncol = state_in%ncol

  ! Retrieve droplet number tendency predicted by the aerosol activation parameterization

  idx = pbuf_get_index('NPCCN') ; call pbuf_get_field(pbuf, idx, npccn)

  ! Transfer info to ptend

  call cnst_get_ind('NUMLIQ', ixnumliq )
  lq(:)        = .false.
  lq(ixnumliq) = .true.
  call physics_ptend_init(ptend, state_in%psetcols, 'actdiag', lq=lq)

  ptend%q(:ncol,:,ixnumliq) = npccn(:ncol,:)

  ! Use npccn to advance NUMLIQ in time; save results to the output state

  call physics_state_copy( state_in, state_out )
  call physics_update( state_out, ptend, dt )

end subroutine tmp_numliq_update_after_activation

end module misc_diagnostics