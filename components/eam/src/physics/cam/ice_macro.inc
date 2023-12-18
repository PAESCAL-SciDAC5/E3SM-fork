
     ! -------------------------------------- !
     ! Ice Saturation Adjustment Computation  !
     ! -------------------------------------- !
     call cnst_get_ind('CLDICE',ixcldice)
     call cnst_get_ind('NUMLIQ',ixnumliq)

     lq2(:)  = .FALSE.
     lq2(1)  = .TRUE.
     lq2(ixcldice) = .TRUE.
     lq2(ixnumice) = .TRUE.

     latsub = latvap + latice

     call physics_ptend_init(ptend_loc, state%psetcols, 'iceadj', ls=.true., lq=lq2 )

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

     ! update local copy of state with the tendencies
     ptend_loc%q(:ncol,top_lev:pver,1)=qvtend(:ncol,top_lev:pver)
     ptend_loc%q(:ncol,top_lev:pver,ixcldice)=qitend(:ncol,top_lev:pver)
     ptend_loc%q(:ncol,top_lev:pver,ixnumice)=initend(:ncol,top_lev:pver)
     ptend_loc%s(:ncol,top_lev:pver)=stend(:ncol,top_lev:pver)

    ! Add the ice tendency to the output tendency
     call physics_ptend_sum(ptend_loc, ptend_all, ncol)

    ! ptend_loc is reset to zero by this call
     call physics_update(state1, ptend_loc, hdtime)

    !Write output for tendencies:
    !        oufld: QVTENDICE,QITENDICE,NITENDICE
     call outfld( 'TTENDICE',  stend/cpair, pcols, lchnk )
     call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
     call outfld( 'QITENDICE', qitend, pcols, lchnk )
     call outfld( 'NITENDICE', initend, pcols, lchnk )
