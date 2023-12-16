
      ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
      ! after CLUBB is called.  This is for energy conservation purposes.
      se_a = 0._r8
      ke_a = 0._r8
      wv_a = 0._r8
      wl_a = 0._r8
      do k=1,pver
         enthalpy = cpair*((thlm(i,k)+(latvap/cpair)*rcm(i,k))/exner_clubb(i,k))
         clubb_s(k) = enthalpy + gravit*state1%zm(i,k)+state1%phis(i)
!         se_a(i) = se_a(i) + clubb_s(k)*state1%pdel(i,k)*invrs_gravit
         se_a(i) = se_a(i) + enthalpy * state1%pdel(i,k)*invrs_gravit
         ke_a(i) = ke_a(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)*invrs_gravit
         wv_a(i) = wv_a(i) + (rtm(i,k)-rcm(i,k))*state1%pdel(i,k)*invrs_gravit
         wl_a(i) = wl_a(i) + (rcm(i,k))*state1%pdel(i,k)*invrs_gravit
      enddo

      ! Based on these integrals, compute the total energy before and after CLUBB call
      ! TE as in Williamson2015, E= \int_{whole domain} (K+c_p*T) +
      ! \int_{surface} p_s\phi_s (up to water forms), but we ignore surface term
      ! under assumption that CLUBB does not change surface pressure
  !   do k=1,pver
         te_a(i) = se_a(i) + ke_a(i) + (latvap+latice)*wv_a(i)+latice*wl_a(i)
  !      te_b(i) = se_b(i) + ke_b(i) + (latvap+latice)*wv_b(i)+latice*wl_b(i)
  !   enddo

  !   ! Take into account the surface fluxes of heat and moisture
  !   te_b(i) = te_b(i)+(cam_in%shf(i)+(cam_in%cflx(i,1))*(latvap+latice))*hdtime

      ! Limit the energy fixer to find highest layer where CLUBB is active
      ! Find first level where wp2 is higher than lowest threshold
      clubbtop = 1
      do while (wp2(i,clubbtop) .eq. w_tol_sqd .and. clubbtop .lt. pver-1)
         clubbtop = clubbtop + 1
      enddo

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis = (te_a(i) - te_b(i))/(state1%pint(i,pverp)-state1%pint(i,clubbtop))

      ! Apply this fixer throughout the column evenly, but only at layers where
      ! CLUBB is active.
      do k=clubbtop,pver
         clubb_s(k) = clubb_s(k) - se_dis*gravit
      enddo
