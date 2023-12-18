
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

