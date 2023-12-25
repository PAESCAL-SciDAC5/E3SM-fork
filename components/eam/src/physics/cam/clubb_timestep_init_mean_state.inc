


   !===========================================================================================
   !  At each CLUBB call, initialize mean momentum  and thermo CLUBB state from the CAM state
   !===========================================================================================
   ! From host model: u,  v;  qv, ql, t, pmid
   ! To CLUBB:        um, vm; rtm, thlm (why rvm also?)
   !--------------------------------------------------------

   call pbuf_get_field(pbuf, um_idx, host_mean%um, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, vm_idx, host_mean%vm, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   host_mean%um(:ncol,:pver) = state1%u(:ncol,:pver)
   host_mean%vm(:ncol,:pver) = state1%v(:ncol,:pver)

   !---------------------------
   call pbuf_get_field(pbuf, thlm_idx, host_mean%thlm, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf,  rtm_idx, host_mean%rtm,  start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf,  rcm_idx, host_mean%rcm,  start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)

   host_mean%rtm(:ncol,:pver) = state1%q(:ncol,:pver,ixq)+state1%q(:ncol,:pver,ixcldliq)

   !  Compute exner function consistent with CLUBB's definition, which uses a constant
   !  surface pressure.  CAM's exner (in state does not).  Therefore, for consistent
   !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables
   !  (such as thlm), use "exner_clubb" other wise use the exner in state

    exner_clubb(:ncol,:pver) = (real(p0_clubb, kind = r8 )/state1%pmid(:ncol,:pver))**(rair/cpair)
    host_mean%thlm(:ncol,:pver) = state1%t(:ncol,:pver)*exner_clubb(:ncol,:pver)-(latvap/cpair)*state1%q(:ncol,:pver,ixcldliq)

   !---------------------------------
   ! Assign value to ghost level
   !---------------------------------

   host_mean%um(1:ncol,pverp) = host_mean%um(1:ncol,pver)
   host_mean%vm(1:ncol,pverp) = host_mean%vm(1:ncol,pver)

   host_mean%rtm(1:ncol,pverp)  = host_mean%rtm(1:ncol,pver)
   host_mean%thlm(1:ncol,pverp) = host_mean%thlm(1:ncol,pver)

   !=====================================================
   ! Passive tracers:
   ! Constituents are all treated as wet mmr by clubb
   ! don't convert co2 tracers to wet mixing ratios
   !=====================================================
   cnst_type_loc(:) = cnst_type(:)
   call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
   call set_dry_to_wet(state1, cnst_type_loc)

