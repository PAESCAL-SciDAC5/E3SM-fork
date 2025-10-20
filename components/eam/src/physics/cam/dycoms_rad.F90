module dycoms_rad

   implicit none
   private
   public :: dycoms_radiation_init
   public :: dycoms_radiation_tend

   contains

subroutine dycoms_radiation_init()

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   ! Initialization specific to the idealized LWP-based radiation calculation
   ! for the DYCOMS-II RF01 case (Stevens et al., 2005).
   !
   !-----------------------------------------------------------------------
   use cam_history,      only: addfld

   ! Add output variables

   call addfld ('DYCOMS_LWPABV',(/ 'ilev' /), 'A', 'kg/m2', 'DYCOMS liquid water path integrated from TOA to the current layer interface')
   call addfld ('DYCOMS_LWFLX', (/ 'ilev' /), 'A', 'W/m2',  'DYCOMS longwave flux profile')
   call addfld ('DYCOMS_QRL',   (/ 'lev' /),  'A', 'K/s',   'DYCOMS longwave heating rate')

end subroutine dycoms_radiation_init


subroutine dycoms_radiation_tend(state, ptend, net_flx)

   use shr_kind_mod,        only: r8 => shr_kind_r8
   use ppgrid,              only: pver, pcols, pverp
   use physics_types,       only: physics_state, &
                                  physics_ptend, &
                                  physics_ptend_init
   use physconst,           only: gravit, cpair
   use constituents,        only: cnst_get_ind
   use cam_history,         only: outfld

   implicit none
   
   type(physics_state), intent(in), target :: state
   type(physics_ptend), intent(out) :: ptend

   real(r8), intent(out) :: net_flx(pcols)
   
   integer :: ncol
   integer :: lchnk
   integer :: ixcldliq

   real(r8) :: lwp_in_layer(pcols)
   real(r8) :: lwp_abv(1:pcols,1:pverp)
   real(r8) :: radflux(1:pcols,1:pverp)
   real(r8) :: qrl(1:pcols,1:pver)

   ! Diagnostic outputs - vertical profiles
   real(r8) :: qrl_diag(1:pcols,1:pver) ! LW heating rate for diagnostics (K/s)

   real(r8), parameter :: F0 = 70.0_r8, &
                          F1 = 22.0_r8, &
                          kappa = 85.0_r8

   integer :: i, k

   !-----------------------------------------------------------------------------------
   ncol = state%ncol
   lchnk = state%lchnk

   ! Calculate liquid water path integrated from model top to a layer interface of interest.
   ! Note that in this highly idealized case, we are not making any assumptions about cloud overlap,
   ! hence grid-box mean cloud liquid mixing ratio is used here.

   call cnst_get_ind('CLDLIQ',ixcldliq)

   lwp_abv(:,:) = 0.0_r8
   do k=2,pverp
      lwp_in_layer(:ncol) = state%q(:ncol,k-1,ixcldliq) * state%pdel(:ncol,k-1) /gravit
      lwp_abv(:ncol,k) = lwp_abv(:ncol,k-1) + lwp_in_layer(:ncol) 
   end do

   ! Calculate upwelling LW flux at layer interfaces
   ! using the first 2 RHS terms in Eq. (3) of Stevens et al. (2005), "Evaluation of
   ! Large-Eddy Simulations via Observations of Nocturnal Marine Stratocumulus"

   do k = 1, pverp
      do i = 1,ncol
         radflux(i, k) = F0*exp(-kappa*lwp_abv(i,k)) &
                       + F1*exp(-kappa*(lwp_abv(i,pverp)-lwp_abv(i,k)))
      end do
   end do

   ! Calculate cp*dT/dt

   do k = 1,pver
      qrl(:ncol,k) = (radflux(:ncol,k+1)-radflux(:ncol,k)) * gravit / state%pdel(:ncol,k)
   end do

   !--------------------------------------------------------------------------
   ! Assign values to arrays needed by the host model
   !--------------------------------------------------------------------------
   call physics_ptend_init(ptend, state%psetcols, 'dycoms_radheat', ls=.true.) 
   ptend%s(:ncol,:) = qrl(:ncol,:)

   ! Net flux into the atm column = net incoming LW at sfc - net outgoing LW at TOP.
   ! SW fluxes are zero.

   net_flx(1:ncol) = radflux(:ncol,pverp) - radflux(:ncol,1)

   !--------------------------------------------------------------------------
   ! Prepare diagnostic outputs following existing radiation scheme patterns
   !--------------------------------------------------------------------------
   ! Heating rate for diagnostics (convert to K/s)
   qrl_diag(:ncol,:pver) = qrl(:ncol,:pver) / cpair

   ! Output diagnostic fields to history buffer following existing patterns
   ! LWP vertical profile (kg/m²) - cumulative from TOA to each level
   call outfld('DYCOMS_LWPABV', lwp_abv(1:ncol,1:pverp), ncol, lchnk)
   
   ! Longwave flux vertical profile (W/m²) at interfaces
   call outfld('DYCOMS_LWFLX', radflux(1:ncol,1:pverp), ncol, lchnk)  
   
   ! Longwave heating rate vertical profile (K/s)
   call outfld('DYCOMS_QRL',  qrl_diag(1:ncol,1:pver), ncol, lchnk)

   return

end subroutine dycoms_radiation_tend

end module dycoms_rad
