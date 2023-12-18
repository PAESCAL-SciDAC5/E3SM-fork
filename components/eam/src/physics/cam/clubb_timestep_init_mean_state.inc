
   !===========================================================================================
   !  At each CLUBB call, initialize mean momentum  and thermo CLUBB state from the CAM state
   !===========================================================================================
   ! From host model: u,  v;  qv, ql, t, pmid
   ! To CLUBB:        um, vm; rtm, thlm (why rvm also?)
   !--------------------------------------------------------

     um(:ncol,:pver) = state1%u(:ncol,:pver)
     vm(:ncol,:pver) = state1%v(:ncol,:pver)

   !---------------------------
   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)

    rtm(:ncol,:pver) = state1%q(:ncol,:pver,ixq)+state1%q(:ncol,:pver,ixcldliq)
    rvm(:ncol,:pver) = state1%q(:ncol,:pver,ixq)

   !  Compute exner function consistent with CLUBB's definition, which uses a constant
   !  surface pressure.  CAM's exner (in state does not).  Therefore, for consistent
   !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables
   !  (such as thlm), use "exner_clubb" other wise use the exner in state

    exner_clubb(:ncol,:pver) = (real(p0_clubb, kind = r8 )/state1%pmid(:ncol,:pver))**(rair/cpair)
           thlm(:ncol,:pver) = state1%t(:ncol,:pver)*exner_clubb(:ncol,:pver)-(latvap/cpair)*state1%q(:ncol,:pver,ixcldliq)

   !---------------------------------
   ! Assign value to ghost level
   !---------------------------------

   um(1:ncol,pverp)   = state1%u(1:ncol,pver)
   vm(1:ncol,pverp)   = state1%v(1:ncol,pver)
   rtm(1:ncol,pverp)  = rtm(1:ncol,pver)
   thlm(1:ncol,pverp) = thlm(1:ncol,pver)

   !-------------------------------------------------
   ! constituents are all treated as wet mmr by clubb
   ! don't convert co2 tracers to wet mixing ratios
   !-------------------------------------------------
   cnst_type_loc(:) = cnst_type(:)
   call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
   call set_dry_to_wet(state1, cnst_type_loc)

