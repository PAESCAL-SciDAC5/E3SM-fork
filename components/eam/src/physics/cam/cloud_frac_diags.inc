
 subroutine cloudfrac_diags(state1, cam_in, cmfmc, itim_old, liqcf_fix, pbuf, cloud_frac, alst_o)

   use physics_types,  only: physics_state
   use camsrfexch,     only: cam_in_t
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field
   use time_manager,   only: is_first_step
   use cldfrc2m,       only: aist_vector
   use ppgrid,         only: pver, pverp, pcols
   use constituents,   only: cnst_get_ind
   use scamMOD,        only: single_column,scm_clubb_iop_name

   type(physics_state),intent(in) :: state1
   type(cam_in_t),     intent(in) :: cam_in
   real(r8),           intent(in) :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c           [kg/m2/s]
   integer,            intent(in) :: itim_old
   logical,            intent(in) :: liqcf_fix

   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8),intent(inout) :: cloud_frac(:,:)

   real(r8), intent(out) :: alst_o(pcols,pver)  ! H. Wang: for old liquid status fraction

   real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction                               [fraction]
   real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: ast      ! stratiform cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: alst     ! liquid stratiform cloud fraction             [fraction]
   real(r8), pointer, dimension(:,:) :: aist     ! ice stratiform cloud fraction                [fraction]
   real(r8), pointer, dimension(:,:) :: qlst     ! Physical in-stratus LWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: qist     ! Physical in-stratus IWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: deepcu   ! deep convection cloud fraction               [fraction]
   real(r8), pointer, dimension(:,:) :: shalcu   ! shallow convection cloud fraction            [fraction]

   real(r8), pointer, dimension(:,:) :: dp_icwmr ! deep convection in cloud mixing ratio        [kg/kg]
   real(r8), pointer, dimension(:,:) :: cmfmc_sh ! Shallow convective mass flux--m subc (pcols,pverp) [kg/m2/s/]


   real(r8), pointer, dimension(:,:) :: rcm        ! CLUBB cloud water mixing ratio               [kg/kg]


   real(r8) :: frac_limit, ic_limit
   real(r8) :: cldthresh                        ! threshold to determin cloud fraction          [kg/kg]

   integer :: ncol, lchnk
   integer :: i,k
   integer :: ixq, ixcldice, ixnumice

   ncol  = state1%ncol
   lchnk = state1%lchnk

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDICE',ixcldice)
   call cnst_get_ind('NUMICE',ixnumice)

   ! --------------------------------------------------------------------------------- !
   !  Diagnose some quantities that are computed in macrop_tend here.                  !
   !  These are inputs required for the microphysics calculation.                      !
   !                                                                                   !
   !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM CLUBB CLOUD FRACTION      !
   ! --------------------------------------------------------------------------------- !
   call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qlst_idx,    qlst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, rcm_idx,     rcm,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))


  ! HW: set alst to alst_o before getting updated
  if(liqcf_fix) then
      if(.not.is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
   endif

   !  initialize variables
   alst(:,:) = 0.0_r8
   qlst(:,:) = 0.0_r8

   do k=1,pver
      do i=1,ncol
         alst(i,k) = cloud_frac(i,k)
         qlst(i,k) = rcm(i,k)/max(0.01_r8,alst(i,k))  ! Incloud stratus condensate mixing ratio
      enddo
   enddo

   ! HW
   if(liqcf_fix) then
      if(is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
   endif
   !HW

   ! --------------------------------------------------------------------------------- !
   !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION                 !
   ! --------------------------------------------------------------------------------- !
   call pbuf_get_field(pbuf, dp_frac_idx, deepcu)
   call pbuf_get_field(pbuf, sh_frac_idx, shalcu)
   deepcu(:,pver) = 0.0_r8
   shalcu(:,pver) = 0.0_r8

   frac_limit = 0.01_r8
     ic_limit = 1.e-12_r8

   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr)
   call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)

   call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   do k=1,pver-1
      do i=1,ncol
         !  diagnose the deep convective cloud fraction, as done in macrophysics based on the
         !  deep convective mass flux, read in from pbuf.  Since shallow convection is never
         !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
         !  fraction is purely from deep convection scheme.
         deepcu(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+500.0_r8*(cmfmc(i,k+1)-cmfmc_sh(i,k+1))),0.6_r8))
         shalcu(i,k) = 0._r8

         if (deepcu(i,k) <= frac_limit .or. dp_icwmr(i,k) < ic_limit) then
            deepcu(i,k) = 0._r8
         endif

         !  using the deep convective cloud fraction, and CLUBB cloud fraction (variable
         !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
         !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud
         !  from CLUBB plus the deep convective cloud fraction
         concld(i,k) = min(cloud_frac(i,k)-alst(i,k)+deepcu(i,k),0.80_r8)
      enddo
   enddo

   if (single_column) then
      if (trim(scm_clubb_iop_name) .eq. 'ATEX_48hr'       .or. &
          trim(scm_clubb_iop_name) .eq. 'BOMEX_5day'      .or. &
          trim(scm_clubb_iop_name) .eq. 'DYCOMSrf01_4day' .or. &
          trim(scm_clubb_iop_name) .eq. 'DYCOMSrf02_06hr' .or. &
          trim(scm_clubb_iop_name) .eq. 'RICO_3day'       .or. &
          trim(scm_clubb_iop_name) .eq. 'ARM_CC') then

             deepcu(:,:) = 0.0_r8
             concld(:,:) = 0.0_r8

      endif
   endif
   call outfld( 'CONCLD',           concld,                  pcols, lchnk )

   ! --------------------------------------------------------------------------------- !
   !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
   !  use the aist_vector function to compute the ice cloud fraction                   !
   ! --------------------------------------------------------------------------------- !
   call pbuf_get_field(pbuf, aist_idx,    aist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call t_startf('ice_cloud_frac_diag')
   do k=1,pver
      call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
           state1%q(:,k,ixnumice),cam_in%landfrac(:),cam_in%snowhland(:),aist(:,k),ncol)
   enddo
   call t_stopf('ice_cloud_frac_diag')

   ! --------------------------------------------------------------------------------- !
   !  THIS PART COMPUTES THE LIQUID STRATUS FRACTION                                   !
   !                                                                                   !
   !  For now leave the computation of ice stratus fraction from macrop_driver intact  !
   !  because CLUBB does nothing with ice.  Here I simply overwrite the liquid stratus !
   !  fraction that was coded in macrop_driver                                         !
   ! --------------------------------------------------------------------------------- !
   call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qist_idx,    qist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   !  Recompute net stratus fraction using maximum over-lapping assumption, as done
   !  in macrophysics code, using alst computed above and aist read in from physics buffer

   cldthresh=1.e-18_r8

   do k=1,pver
      do i=1,ncol

         ast(i,k) = max(alst(i,k),aist(i,k))

         qist(i,k) = state1%q(i,k,ixcldice)/max(0.01_r8,aist(i,k))
      enddo
   enddo

   !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just
   !  be outputting the shallow convective cloud fraction

   do k=1,pver
      do i=1,ncol
         cloud_frac(i,k) = min(ast(i,k)+deepcu(i,k),1.0_r8)
      enddo
   enddo
   call outfld( 'CLOUDCOVER_CLUBB', cloud_frac, pcols, lchnk )

   ! Assign the first pver levels of cloud_frac back to cld

   call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   cld(:,1:pver) = cloud_frac(:,1:pver)

 end subroutine cloudfrac_diags
