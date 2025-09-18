module dycoms_rad

   implicit none
   private
   public :: dycoms_radiation_tend

   contains

subroutine dycoms_radiation_tend(state, ptend, pbuf, ztodt)

   use shr_kind_mod,        only: r8 => shr_kind_r8
   use ppgrid,              only: pver, pcols, pverp
   use physics_buffer,      only: physics_buffer_desc
   use physics_types,       only: physics_state, &
                                  physics_ptend, &
                                  physics_ptend_init
   use physconst,           only: gravit, cpair
   use physics_buffer,      only: pbuf_get_index, pbuf_get_field
   use cam_history,         only: outfld

   implicit none
   
   type(physics_state), intent(in), target :: state
   type(physics_ptend), intent(out) :: ptend
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in) :: ztodt
   
   !-----------------------------------------------------------------------------------

   integer :: ncol
   integer :: lchnk

   real(r8) :: pint(1:pcols,1:pver+1)

   ! array shape correct?
   real(r8) :: lwp(1:pcols,1:pver+1)
   real(r8) :: radflux(1:pcols,1:pver+1)
   real(r8) :: qrl(1:pcols,1:pver)

   ! Diagnostic outputs - vertical profiles
   real(r8) :: qrl_diag(1:pcols,1:pver) ! LW heating rate for diagnostics (K/s)

   real(r8), parameter :: F0 = 70.0_r8, &
                          F1 = 22.0_r8, &
                          kappa = 85.0_r8

   real(r8), pointer, dimension(:,:) :: iclwp, cld

   integer :: i, k

   ncol = state%ncol
   lchnk = state%lchnk
   pint(1:ncol,1:pver+1) = state%pint(1:ncol,1:pver+1)

   call pbuf_get_field(pbuf, pbuf_get_index('ICLWP'), iclwp)
   call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld)

   lwp = 0.0_r8
   do k=2,pver+1
      do i = 1,ncol
         lwp(i,k) = lwp(i,k-1) + iclwp(i,k-1)*cld(i,k-1)
      end do
   end do

   do k = 1, pver+1
      do i = 1,ncol
         radflux(i, k) = F0*exp(-kappa*lwp(i,k)) + F1*exp(-kappa*(lwp(i,pver+1)-lwp(i,k)))
      end do
   end do

   do k = 1,pver
      do i = 1,ncol
         qrl(i,k) = (radflux(i,k+1)-radflux(i,k))*gravit/(pint(i,k+1)-pint(i,k))
      end do
   end do

   call physics_ptend_init(ptend, state%psetcols, 'dycoms_radheat', ls=.true.) 
   ptend%s(:ncol,:) = qrl(:ncol,:)

   ! Prepare diagnostic outputs following existing radiation scheme patterns
   
   ! Heating rate for diagnostics (convert to K/s)
   do k = 1,pver
      do i = 1,ncol
         qrl_diag(i,k) = qrl(i,k) / cpair
      end do
   end do

   ! Output diagnostic fields to history buffer following existing patterns
   ! LWP vertical profile (kg/m²) - cumulative from TOA to each level
   call outfld('DYCOMS_LWP',   lwp(1:ncol,1:pver+1),     ncol, lchnk)
   
   ! Longwave flux vertical profile (W/m²) at interfaces
   call outfld('DYCOMS_LWFLX', radflux(1:ncol,1:pver+1), ncol, lchnk)  
   
   ! Longwave heating rate vertical profile (K/s)
   call outfld('DYCOMS_QRL',   qrl_diag(1:ncol,1:pver),  ncol, lchnk)

   return

end subroutine dycoms_radiation_tend

end module dycoms_rad