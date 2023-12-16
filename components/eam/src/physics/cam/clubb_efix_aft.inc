
      ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
      ! after CLUBB is called.  This is for energy conservation purposes.

      se_a = 0._r8
      ke_a = 0._r8
      wv_a = 0._r8
      wl_a = 0._r8

      do k=1,pver
         enthalpy = cpair*((thlm(i,k)+(latvap/cpair)*rcm(i,k))/exner_clubb(i,k))
         clubb_s(k) = enthalpy + gravit*state1%zm(i,k)+state1%phis(i)
        !se_a = se_a + clubb_s(k)*state1%pdel(i,k)*invrs_gravit
         se_a = se_a + enthalpy * state1%pdel(i,k)*invrs_gravit
         ke_a = ke_a + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)*invrs_gravit
         wv_a = wv_a + (rtm(i,k)-rcm(i,k))*state1%pdel(i,k)*invrs_gravit
         wl_a = wl_a + (rcm(i,k))*state1%pdel(i,k)*invrs_gravit
      enddo

      ! Based on these integrals, compute the total energy before and after CLUBB call
      ! TE as in Williamson2015, E= \int_{whole domain} (K+c_p*T) +
      ! \int_{surface} p_s\phi_s (up to water forms), but we ignore surface term
      ! under assumption that CLUBB does not change surface pressure
      te_a = se_a + ke_a + (latvap+latice)*wv_a + latice*wl_a

      ! Limit the energy fixer to find highest layer where CLUBB is active
      ! Find first level where wp2 is higher than lowest threshold
      clubbtop = 1
      do while (wp2(i,clubbtop) .eq. w_tol_sqd .and. clubbtop .lt. pver-1)
         clubbtop = clubbtop + 1
      enddo

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis = (te_a - te_b)/(state1%pint(i,pverp)-state1%pint(i,clubbtop))

      ! Apply this fixer throughout the column evenly, but only at layers where
      ! CLUBB is active.
      do k=clubbtop,pver
         clubb_s(k) = clubb_s(k) - se_dis*gravit
      enddo
