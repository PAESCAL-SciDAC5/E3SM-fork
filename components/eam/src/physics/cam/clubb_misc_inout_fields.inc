

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
