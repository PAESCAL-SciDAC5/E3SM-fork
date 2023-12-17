
!PMA adds gustiness and tpert

    call pbuf_get_field(pbuf, tpert_idx,   tpert)
    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
    call pbuf_get_field(pbuf, vmag_gust_idx, vmag_gust)

    vmag_gust(:)    = 0._r8
    vmag_gust_dp(:) = 0._r8
    vmag_gust_cl(:) = 0._r8
    ktopi(:)        = pver

    if (use_sgv) then
       do i=1,ncol
           up2b(i)          = up2(i,pver)
           vp2b(i)          = vp2(i,pver)
           umb(i)           = state1%u(i,pver)
           vmb(i)           = state1%v(i,pver)
           prec_gust(i)     = max(0._r8,prec_dp(i)-snow_dp(i))*1.e3_r8
           if (cam_in%landfrac(i).gt.0.95_r8) then
             gust_fac(i)   = gust_facl
           else
             gust_fac(i)   = gust_faco
           endif
           vmag(i)         = max(1.e-5_r8,sqrt( umb(i)**2._r8 + vmb(i)**2._r8))
           vmag_gust_dp(i) = ugust(min(prec_gust(i),6.94444e-4_r8),gust_fac(i)) ! Limit for the ZM gustiness equation set in Redelsperger et al. (2000)
           vmag_gust_dp(i) = max(0._r8, vmag_gust_dp(i) )!/ vmag(i))
           vmag_gust_cl(i) = gust_facc*(sqrt(max(0._r8,up2b(i)+vp2b(i))+vmag(i)**2._r8)-vmag(i))
           vmag_gust_cl(i) = max(0._r8, vmag_gust_cl(i) )!/ vmag(i))
           vmag_gust(i)    = vmag_gust_cl(i) + vmag_gust_dp(i)
          do k=1,pver
             if (state1%zi(i,k)>pblh(i).and.state1%zi(i,k+1)<=pblh(i)) then
                ktopi(i) = k
                exit
             end if
          end do
          tpert(i) = min(2._r8,(sqrt(thlp2(i,ktopi(i)))+(latvap/cpair)*state1%q(i,ktopi(i),ixcldliq)) &
                    /max(state1%exner(i,ktopi(i)),1.e-3_r8)) !proxy for tpert
       end do
    end if

   call outfld('VMAGGUST', vmag_gust, pcols, lchnk)
   call outfld('VMAGDP', vmag_gust_dp, pcols, lchnk)
   call outfld('VMAGCL', vmag_gust_cl, pcols, lchnk)
   call outfld('TPERTBLT', tpert, pcols, lchnk)


