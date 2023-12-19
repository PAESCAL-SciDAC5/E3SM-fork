
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Compute integrals of static energy, kinetic energy, water vapor, and liquid water
      ! for the computation of total energy before CLUBB is called.  This is for an
      ! effort to conserve energy since liquid water potential temperature (which CLUBB
      ! conserves) and static energy (which CAM conserves) are not exactly equal.

      invrs_gravit = 1._r8 / gravit

      se_b = 0._r8  ! initialize vertical integrals
      ke_b = 0._r8
      wv_b = 0._r8
      wl_b = 0._r8

      do k=1,pver ! vertical integral
         ! use s=c_pT+g*z, total energy needs term c_pT but not gz
         se_b = se_b + (state1%s(i,k) - gravit*state1%zm(i,k) - state1%phis(i)) &
                      *  state1%pdel(i,k)*invrs_gravit
         ke_b = ke_b + 0.5_r8*(state1%u(i,k)**2+state1%v(i,k)**2)*state1%pdel(i,k)*invrs_gravit
         wv_b = wv_b + state1%q(i,k,ixq)*state1%pdel(i,k)*invrs_gravit
         wl_b = wl_b + state1%q(i,k,ixcldliq)*state1%pdel(i,k)*invrs_gravit
      enddo

      ! Total energy: sum up the components and also take into account the surface fluxes of heat and moisture
      te_b = se_b + ke_b + (latvap+latice)*wv_b + latice*wl_b
      te_b = te_b +(cam_in%shf(i)+(cam_in%cflx(i,1))*(latvap+latice))*hdtime
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
