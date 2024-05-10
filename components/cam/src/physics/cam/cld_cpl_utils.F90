module cld_cpl_utils

  use shr_kind_mod,     only: r8 => shr_kind_r8

  implicit none
  public

contains

   subroutine cld_cpl_register( cld_cpl_opt )

     use physics_buffer, only: pbuf_add_field, dtype_r8
     use ppgrid,         only: pcols, pver

     integer, intent(in) :: cld_cpl_opt

     integer :: idxtmp   ! pbuf component index

     if (cld_cpl_opt > 1) then
       call pbuf_add_field('S_After_MACMIC' ,  'global', dtype_r8, (/pcols,pver/), idxtmp)
       call pbuf_add_field('T_After_MACMIC' ,  'global', dtype_r8, (/pcols,pver/), idxtmp)
       call pbuf_add_field('Q_After_MACMIC' ,  'global', dtype_r8, (/pcols,pver/), idxtmp)
       call pbuf_add_field('QL_After_MACMIC',  'global', dtype_r8, (/pcols,pver/), idxtmp)
       call pbuf_add_field('QI_After_MACMIC',  'global', dtype_r8, (/pcols,pver/), idxtmp)
       call pbuf_add_field('NL_After_MACMIC',  'global', dtype_r8, (/pcols,pver/), idxtmp)
       call pbuf_add_field('NI_After_MACMIC',  'global', dtype_r8, (/pcols,pver/), idxtmp)
     end if

   end subroutine cld_cpl_register

   !---------------------------------------------------------------------------------------
   subroutine set_state_and_tendencies( state, pbuf, ztodt, ptend_dribble )

   use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
   use physics_buffer,   only: pbuf_get_index
   use constituents,     only: pcnst, cnst_get_ind
   use ppgrid,           only: pcols, pver

   type(physics_state), intent(inout) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8),intent(in) :: ztodt

   type(physics_ptend),intent(out) :: ptend_dribble

   ! local variables

   integer :: ifld
   integer :: ncol
   integer :: ixcldliq, ixcldice, ixnumliq, ixnumice, ixq
   logical :: lq(pcnst)

    real(r8), pointer, dimension(:,:) :: s_after_macmic
    real(r8), pointer, dimension(:,:) :: t_after_macmic
    real(r8), pointer, dimension(:,:) :: q_after_macmic
    real(r8), pointer, dimension(:,:) :: ql_after_macmic
    real(r8), pointer, dimension(:,:) :: qi_after_macmic
    real(r8), pointer, dimension(:,:) :: nl_after_macmic
    real(r8), pointer, dimension(:,:) :: ni_after_macmic

   !-----------------------------------------------
   ncol = state%ncol

   ! Retrieve an old state saved in pbuf

   ifld = pbuf_get_index( 'S_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  s_after_macmic )
   ifld = pbuf_get_index( 'T_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  t_after_macmic )
   ifld = pbuf_get_index( 'Q_After_MACMIC'); call pbuf_get_field(pbuf, ifld,  q_after_macmic )
   ifld = pbuf_get_index('QL_After_MACMIC'); call pbuf_get_field(pbuf, ifld, ql_after_macmic )
   ifld = pbuf_get_index('QI_After_MACMIC'); call pbuf_get_field(pbuf, ifld, qi_after_macmic )
   ifld = pbuf_get_index('NL_After_MACMIC'); call pbuf_get_field(pbuf, ifld, nl_after_macmic )
   ifld = pbuf_get_index('NI_After_MACMIC'); call pbuf_get_field(pbuf, ifld, ni_after_macmic )

   ! Calculate tendencies to be dribbled, and save in ptend_dribble

   call cnst_get_ind('Q',      ixq)
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('NUMLIQ', ixnumliq)
   call cnst_get_ind('NUMICE', ixnumice)

   lq(:)        = .FALSE.
   lq(ixq)      = .TRUE.
   lq(ixcldliq) = .TRUE.
   lq(ixcldice) = .TRUE.
   lq(ixnumliq) = .TRUE.
   lq(ixnumice) = .TRUE.

   call physics_ptend_init(ptend_dribble, state%psetcols, 'macmic_dribble_tend', ls= .true., lq=lq)

   ptend_dribble%s(:ncol,:pver)          = (state%s(:ncol,:pver)          -   s_after_macmic(:ncol,:pver))  / ztodt
   ptend_dribble%q(:ncol,:pver,ixq)      = (state%q(:ncol,:pver,ixq)      -   q_after_macmic(:ncol,:pver))  / ztodt
   ptend_dribble%q(:ncol,:pver,ixcldliq) = (state%q(:ncol,:pver,ixcldliq) -  ql_after_macmic(:ncol,:pver))  / ztodt
   ptend_dribble%q(:ncol,:pver,ixcldice) = (state%q(:ncol,:pver,ixcldice) -  qi_after_macmic(:ncol,:pver))  / ztodt
   ptend_dribble%q(:ncol,:pver,ixnumliq) = (state%q(:ncol,:pver,ixnumliq) -  nl_after_macmic(:ncol,:pver))  / ztodt
   ptend_dribble%q(:ncol,:pver,ixnumice) = (state%q(:ncol,:pver,ixnumice) -  ni_after_macmic(:ncol,:pver))  / ztodt

   ! Reset state back to the snapshot at old time step  

   state%s(:ncol,:pver)          =  s_after_macmic(:ncol,:pver)
   state%t(:ncol,:pver)          =  t_after_macmic(:ncol,:pver)
   state%q(:ncol,:pver,1)        =  q_after_macmic(:ncol,:pver)
   state%q(:ncol,:pver,ixcldliq) = ql_after_macmic(:ncol,:pver)
   state%q(:ncol,:pver,ixcldice) = qi_after_macmic(:ncol,:pver)
   state%q(:ncol,:pver,ixnumliq) = nl_after_macmic(:ncol,:pver)
   state%q(:ncol,:pver,ixnumice) = ni_after_macmic(:ncol,:pver)
 
   end subroutine set_state_and_tendencies

   !-------------------------------------------------------
   subroutine save_state_snapshot_to_pbuf( state, pbuf )

   use physics_types,  only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
   use ppgrid,         only: pver
   use constituents,   only: cnst_get_ind

   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), pointer, dimension(:,:) :: ptr2d
   integer :: ifld
   integer :: ncol
   integer :: ixcldliq, ixcldice, ixnumliq, ixnumice, ixq

   ncol = state%ncol

   call cnst_get_ind('Q',      ixq)
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('NUMLIQ', ixnumliq)
   call cnst_get_ind('NUMICE', ixnumice)

   ifld = pbuf_get_index('S_After_MACMIC')
   call pbuf_get_field(pbuf, ifld, ptr2d)
   ptr2d(:ncol,:pver) = state%s(:ncol,:pver)

   ifld = pbuf_get_index('T_After_MACMIC')
   call pbuf_get_field(pbuf, ifld, ptr2d)
   ptr2d(:ncol,:pver) = state%t(:ncol,:pver)

   ifld = pbuf_get_index('Q_After_MACMIC')
   call pbuf_get_field(pbuf, ifld, ptr2d)
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixq)

   ifld = pbuf_get_index('QL_After_MACMIC')
   call pbuf_get_field(pbuf, ifld, ptr2d )
   ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)

   ifld = pbuf_get_index('QI_After_MACMIC')
   call pbuf_get_field(pbuf, ifld, ptr2d )
   ptr2d(:ncol,:pver)  = state%q(:ncol,:pver,ixcldice)

   ifld = pbuf_get_index('NL_After_MACMIC')
   call pbuf_get_field(pbuf, ifld, ptr2d )
   ptr2d(:ncol,:pver)  = state%q(:ncol,:pver,ixnumliq)

   ifld = pbuf_get_index('NI_After_MACMIC')
   call pbuf_get_field(pbuf, ifld, ptr2d )
   ptr2d(:ncol,:pver)  = state%q(:ncol,:pver,ixnumice)

   end subroutine save_state_snapshot_to_pbuf


end module cld_cpl_utils
