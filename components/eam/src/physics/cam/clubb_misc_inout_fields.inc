
   call pbuf_get_field(pbuf, pdf_zm_w_1_idx, pdf_zm_w_1, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_w_2_idx, pdf_zm_w_2, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_varnce_w_1_idx, pdf_zm_varnce_w_1, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_varnce_w_2_idx, pdf_zm_varnce_w_2, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
   call pbuf_get_field(pbuf, pdf_zm_mixt_frac_idx, pdf_zm_mixt_frac, start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

   call pbuf_get_field(pbuf, qrl_idx,      qrl)
   call pbuf_get_field(pbuf, radf_idx,     radf_clubb)
   call pbuf_get_field(pbuf, prer_evap_idx,prer_evap)

   call pbuf_get_field(pbuf, kvm_idx,     khzt)
   call pbuf_get_field(pbuf, kvh_idx,     khzm) 

   !===========================
   if (linearize_pbl_winds) then
      call pbuf_get_field(pbuf, um_pert_idx, um_pert, start=(/1,1/),          kount=(/pcols,pverp/))
      call pbuf_get_field(pbuf, vm_pert_idx, vm_pert, start=(/1,1/),          kount=(/pcols,pverp/))
      call pbuf_get_field(pbuf, upwp_pert_idx,upwp_pert,start=(/1,1/),        kount=(/pcols,pverp/))
      call pbuf_get_field(pbuf, vpwp_pert_idx,vpwp_pert,start=(/1,1/),        kount=(/pcols,pverp/))
      call pbuf_get_field(pbuf, wsresp_idx, wsresp)
      call pbuf_get_field(pbuf, tau_est_idx, tau_est)
   else
      nullify(um_pert)
      nullify(vm_pert)
      nullify(upwp_pert)
      nullify(vpwp_pert)
   end if
   !===========================
