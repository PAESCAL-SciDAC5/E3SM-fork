module macrop_driver_with_clubb

  use perf_mod,      only: t_startf, t_stopf
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use cam_history,   only: outfld
  use ppgrid,        only: pver, pcols
  use constituents,  only: pcnst

  private
  public :: ice_supersat_adj_tend

contains

  !---------------------------------------
  ! Ice Saturation Adjustment Computation 
  !---------------------------------------
  subroutine ice_supersat_adj_tend( state1, pbuf, hdtime, ptend_loc)

    use physics_types,  only: physics_state, physics_ptend, &
                              physics_state_copy, physics_ptend_init
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
    use constituents,   only: cnst_get_ind
    use physconst,      only: latvap, latice
    use ref_pres,       only: top_lev => trop_cloud_top_lev
    use macrop_driver,  only: ice_macro_tend

    type(physics_state),intent(in)    :: state1
    type(physics_buffer_desc),pointer :: pbuf(:)
    real(r8),intent(in)               :: hdtime

    type(physics_ptend),intent(out)   :: ptend_loc   ! Local tendency from process

    integer :: ncol, lchnk
    logical :: lq2(pcnst)
    integer :: ixcldice, ixnumliq
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
     call cnst_get_ind('NUMLIQ',ixnumliq)

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

end module macrop_driver_with_clubb
