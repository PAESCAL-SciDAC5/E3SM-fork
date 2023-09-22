module mo_spitfire_transport

!---------------------------------------------------------------------------------
! Purpose:
!
! Contains routines to compute tendencies from sedimentation
!
! Author: Phil Rasch
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, rair
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun

  private

  public :: getflx

contains

!===============================================================================
! Calculate tracer fluxes across cell boundaries using the 1D SPITFIRE
! (SPlit Implementation of Transport using Flux Integral REpresentation) 
! algorithm of Rasch and Lawrence (1998):
!   Rasch, P. J., and M. Lawrence, Recent development in transport methods
!   at NCAR, MPI-Rep. 265, pp. 65 – 75, Max-Planck-Inst. fuer Meteorol.,
!   Hamburg, Germany, 1998. 
!===============================================================================
  subroutine getflx(ncol, xw, phi, vel, deltat, flux)

!-----------------------------------------------------------------
! Assumed grid staggering:
!
! Input:
!
!     xw1.......xw2.......xw3.......xw4.......xw5.......xw6
!    vel1......vel2......vel3......vel4......vel5......vel6
!     ....phi1......phi2.......phi3.....phi4.......phi5....
!
! Work arrays:
!
!    psi1......psi2......psi3......psi4......psi5......psi6
!
! Output:
!
!   flux1.....flux2.....flux3.....flux4.....flux5.....flux6
!-----------------------------------------------------------------

    implicit none

    integer, intent(in) :: ncol     ! number of colums to process

    real(r8),intent(in) ::  xw(pcols,pverp)     ! coordinate variable, values at layer interfaces. In EAM this is pint.
    real(r8),intent(in) :: phi(pcols,pverp-1)   ! grid cell mean tracer mixing ratio
    real(r8),intent(in) :: vel(pcols,pverp)     ! velocity in the xw coordinate. In EAM this is grav * rho * v
                                                ! where v is velocity in the height (z) coordinate 
    real(r8),intent(in) :: deltat

    real(r8),intent(out) :: flux(pcols,pverp)

    real(r8) psi(pcols,pverp)    ! integral of phi along the xw coordinate
    real(r8) xxk(pcols,pver)     ! departure point in the xw coordinate
    real(r8) psistar(pcols)      ! integral of phi at departure point
    real(r8) fdot(pcols,pverp)   ! derivative of interpolating polynomial

    integer :: kk  ! vertical layer index

    ! Set fluxes at boundaries to zero

    flux(:ncol,1) = 0
    flux(:ncol,pverp) = 0._r8

    ! Get the vertical integral of phi. 
    ! See Rasch and Lawrence (1998), Eq (3) but note we are using a pressure coordinate here.

    psi(:ncol,1) = 0._r8
    do kk = 2,pverp
       psi(:ncol,kk) = phi(:ncol,kk-1)*( xw(:ncol,kk)-xw(:ncol,kk-1) ) + psi(:ncol,kk-1)
    end do

    ! Calculate the derivatives for the interpolating polynomial

    call cfdotmc_pro (ncol, xw, psi, &! in
                      fdot           )! out

    ! Calculate fluxes at interior interfaces

    do kk = 2,pver

       ! Find departure point. Rasch and Lawrence (1998), Eq (4)

       xxk(:ncol,kk) = xw(:ncol,kk)-vel(:ncol,kk)*deltat   

       ! Calculate the integral, psistar, at the departure point xxk.

       call cfint2(ncol, xw, psi, fdot, xxk(:,kk), &! in
                   psistar                         )! out

       ! Calculate the flux at interface kk. Rasch and Lawrence (1998), Eq (5)

       flux(:ncol,kk) = psi(:ncol,kk)-psistar(:ncol)
    end do

  end subroutine getflx



!##############################################################################
! Given a coordinate xw, an interpolating polynomial ff and its derivative fdot,
! calculate the value of the polynomial (psistar) at xin. 
!##############################################################################

  subroutine cfint2 (ncol, xw, ff, fdot, xin, psistar)

    implicit none

    integer ncol                      ! number of colums to process

    real(r8),intent(in) ::   xw(pcols, pverp)
    real(r8),intent(in) ::   ff(pcols, pverp)
    real(r8),intent(in) :: fdot(pcols, pverp)
    real(r8),intent(in) ::  xin(pcols)

    real(r8),intent(out) :: psistar(pcols)

    real(r8) fxdot(pcols)
    real(r8) fxdd(pcols)

    integer ii
    integer kk
    integer intz(pcols)
    real (r8) dx
    real (r8) ss
    real (r8) c2
    real (r8) c3
    real (r8) xx
    real (r8) xinf
    real (r8) psi1, psi2, psi3, psim
    real (r8) cfint
    real (r8) cfnew
    real (r8) xins(pcols)

    do ii = 1,ncol
       xins(ii) = medan(xw(ii,1), xin(ii), xw(ii,pverp))
       intz(ii) = 0
    end do

    ! first find the interval 
    do kk =  1,pverp-1
    do ii = 1,ncol
       if ((xins(ii)-xw(ii,kk))*(xw(ii,kk+1)-xins(ii)).ge.0._r8) then
          intz(ii) = kk
       endif
    end do
    end do

    do ii = 1,ncol
       if (intz(ii).eq.0) then
          write(iulog,*) ' interval was not found for col i ', ii
          write(iulog,*) ' interval was not found for value ', xins(ii) !added by Jianfeng Li
          call endrun('mo_spitfire_transport: cfint2 -- interval was not found ')
       endif
    end do

    ! now interpolate

    do ii = 1,ncol
       kk = intz(ii)
       dx = (xw(ii,kk+1)-xw(ii,kk))
       ss = (ff(ii,kk+1)-ff(ii,kk))/dx

       !cubic Hermite interpolant coefficient, Eq 2.2 on Page 5 of Huynh (1991)
       c2 = (3*ss-2*fdot(ii,kk)-fdot(ii,kk+1))/dx
       c3 = (fdot(ii,kk)+fdot(ii,kk+1)-2*ss)/dx**2
       xx = (xins(ii)-xw(ii,kk))
       fxdot(ii) =  (3*c3*xx + 2*c2)*xx + fdot(ii,kk) !first-order derivative of q(x)
       fxdd (ii) = 6*c3*xx + 2*c2 !second-order derivative of q(x)
       cfint = ((c3*xx + c2)*xx + fdot(ii,kk))*xx + ff(ii,kk) !q(x)

       ! limit the interpolant

       psi1 = ff(ii,kk)+(ff(ii,kk+1)-ff(ii,kk))*xx/dx

       if (kk.eq.1) then
          psi2 = ff(ii,1)
       else
          psi2 = ff(ii,kk) + (ff(ii,kk)-ff(ii,kk-1))*xx/(xw(ii,kk)-xw(ii,kk-1))
       endif

       if (kk+1.eq.pverp) then
          psi3 = ff(ii,pverp)
       else
          psi3 = ff(ii,kk+1) - (ff(ii,kk+2)-ff(ii,kk+1))*(dx-xx)/(xw(ii,kk+2)-xw(ii,kk+1))
       endif

       psim = medan(psi1, psi2, psi3)
       cfnew = medan(cfint, psi1, psim)

       !if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36_r8)  .gt..03_r8) then
       !   write(iulog,*) ' cfint2 limiting important ', cfint, cfnew 
       !endif

       psistar(ii) = cfnew

    end do

  end subroutine cfint2


!##############################################################################
! Calculate the derivative for the interpolating polynomial.
! Multi column version.
! Reference: Hung T. Huync, Accurate Monotone Cubic Interpolation (1991). NASA Technical Memorandum 103789, 61 pages.
!##############################################################################

  subroutine cfdotmc_pro (ncol, xw, ff, fdot)

    implicit none

    integer ncol                            ! number of colums to process

    real(r8),intent(in)  ::   xw(pcols,pverp)   ! coordinate variable
    real(r8),intent(in)  ::   ff(pcols,pverp)   ! value at notes 
    real(r8),intent(out) :: fdot(pcols,pverp)   ! derivative at nodes

! Assumed variable distribution (staggering)
!     xw1.....xw2......xw3......xw4......xw5.......xw6    1,pverp points
!     ff1.....ff2......ff3......ff4......ff5.......ff6    1,pverp points
!     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
!     ........dd2......dd3......dd4......dd5.........     2,pver points
!     ........ss2......ss3......ss4......ss5.........     2,pver points
!     .............dh2......dh3......dh4.............     2,pver-1 points
!     .............eh2......eh3......eh4.............     2,pver-1 points
!     .................ee3.......ee4.................     3,pver-1 points
!     .................ppl3......ppl4................     3,pver-1 points
!     .................ppr3......ppr4................     3,pver-1 points
!     .................tt3.......tt4.................     3,pver-1 points
!     ................fdot3.....fdot4................     3,pver-1 points


    integer ii
    integer kk

    real (r8) aa                   ! work var
    real (r8) bb                   ! work var
    real (r8) cc                   ! work var
    real (r8) ss(pcols,pverp)      ! first divided differences at nodes
    real (r8) sh(pcols,pverp)      ! first divided differences between nodes
    real (r8) dd(pcols,pverp)      ! second divided differences at nodes
    real (r8) dh(pcols,pverp)      ! second divided differences between nodes
    real (r8) ee(pcols,pverp)      ! third divided differences at nodes
    real (r8) eh(pcols,pverp)      ! third divided differences between nodes
    real (r8) pp                   ! p prime
    real (r8) ppl(pcols,pverp)     ! p prime on left
    real (r8) ppr(pcols,pverp)     ! p prime on right
    real (r8) qpl
    real (r8) qpr
    real (r8) ttlmt
    real (r8) tt
    real (r8) tmin
    real (r8) tmax
    real (r8) delxh(pcols,pverp)


    !-----------------
    do kk = 1,pver

       ! First divided differences between nodes

       do ii = 1, ncol
          delxh(ii,kk) = xw(ii,kk+1) - xw(ii,kk)

          !first divided difference, Eq (2.1) on Page 4 of Huynh (1991)
          sh(ii,kk) = (ff(ii,kk+1)-ff(ii,kk))/delxh(ii,kk)
       end do

       ! First and second divided differences at nodes

       if (kk.ge.2) then
          do ii = 1,ncol
             !second divided difference, Eq (3.1) on Page 19 of Huynh (1991)
             dd(ii,kk) = (sh(ii,kk)-sh(ii,kk-1))/(xw(ii,kk+1)-xw(ii,kk-1))

             !Eq (2.8) on Page 8 of Huynh (1991)
             ss(ii,kk) = minmod(sh(ii,kk),sh(ii,kk-1))
          end do
       endif

    end do

    ! Second and third divided diffs between nodes

    do kk = 2,pver-1
    do ii = 1,ncol
          !third divided difference, Eq (4.1) on Page 36 of Huynh (1991)
          eh(ii,kk) = ( dd(ii,kk+1)-dd(ii,kk) )/( xw(ii,kk+2)-xw(ii,kk-1) )

          !Eq (3.4) on Page 19 of Huynh (1991)
          dh(ii,kk) = minmod( dd(ii,kk), dd(ii,kk+1) )
    end do
    end do

    ! Treat the boundaries

    do ii = 1,ncol

       ee(ii,2)    = eh(ii,2)
       ee(ii,pver) = eh(ii,pver-1)

       ! Outside level
       ! The first equation after Eq 5.1 on Page 42 of Huynh (1991)
       fdot(ii,1) = sh(ii,1) - dd(ii,2)*delxh(ii,1) - eh(ii,2)*delxh(ii,1)*(xw(ii,1)-xw(ii,3))
       
       !Eq 5.2 on Page 43 of Huynh (1991)
       fdot(ii,1) = minmod( fdot(ii,1), 3*sh(ii,1))

       fdot(ii,pverp) = sh(ii,pver) + dd(ii,pver)*delxh(ii,pver)  &
                      + eh(ii,pver-1)*delxh(ii,pver)*(xw(ii,pverp)-xw(ii,pver-1))
       !Eq 5.2 on Page 43 of Huynh (1991)
       fdot(ii,pverp) = minmod( fdot(ii,pverp), 3*sh(ii,pver) )

       ! One in from boundary
       ! The second equation after Eq 5.1 on Page 42 of Huynh (1991)
       fdot(ii,2) = sh(ii,1) + dd(ii,2)*delxh(ii,1) - eh(ii,2)*delxh(ii,1)*delxh(ii,2)

       fdot(ii,2) = minmod(fdot(ii,2),3*ss(ii,2))

       fdot(ii,pver) = sh(ii,pver) - dd(ii,pver)*delxh(ii,pver)   &
                     - eh(ii,pver-1)*delxh(ii,pver)*delxh(ii,pver-1)
       fdot(ii,pver) = minmod( fdot(ii,pver), 3*ss(ii,pver) )

    end do


    do kk = 3,pver-1
    do ii = 1,ncol
       !Eq (4.4) on Page 37 of Huynh (1991)
       ee(ii,kk) = minmod( eh(ii,kk), eh(ii,kk-1) )
    end do
    end do


    do kk = 3,pver-1
    do ii = 1,ncol

       ! p prime at k-0.5
       ! The first Equation after Eq 3.16' on Page 24 of Huynh (1991)
       ppl(ii,kk)=sh(ii,kk-1) + dh(ii,kk-1)*delxh(ii,kk-1)  

       ! The second Equation after Eq 3.16' on Page 24 of Huynh (1991)
       ! p prime at k+0.5
       ppr(ii,kk)=sh(ii,kk)   - dh(ii,kk)  *delxh(ii,kk)

       ! The third Equation after Eq 3.16'  on Page 24 of Huynh (1991)
       tt = minmod(ppl(ii,kk),ppr(ii,kk))

       ! derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
       ! The second Equation of Eq 3.5b on Page 20 of Huynh (1991)
       pp = sh(ii,kk-1) + dd(ii,kk)*delxh(ii,kk-1) 

       ! quartic estimate of fdot
       ! Eq 4.21 on Page 42 of Huynh (1991)
       fdot(ii,kk) = pp                                      &
                   - delxh(ii,kk-1)*delxh(ii,kk)             &
                   *(  eh(ii,kk-1)*(xw(ii,kk+2)-xw(ii,kk  )) &
                     + eh(ii,kk  )*(xw(ii,kk  )-xw(ii,kk-2)) &
                   )/(xw(ii,kk+2)-xw(ii,kk-2))

       ! now limit it
       ! Eq 4.8b on Page 38 of Huynh (1991)
       qpl = sh(ii,kk-1)       &
           + delxh(ii,kk-1)*minmod( dd(ii,kk-1)+ee(ii,kk-1)*(xw(ii,kk)-xw(ii,kk-2)), &
                                    dd(ii,kk)  -ee(ii,kk  )*delxh(ii,kk)             )

       ! Eq 4.8a on Page 38 of Huynh (1991)
       !qpr = sh - delxh ... ? by Jianfeng Li                     
       qpr = sh(ii,kk)         &
           + delxh(ii,kk  )*minmod(dd(ii,kk)  +ee(ii,kk)*delxh(ii,kk-1),        &
                                   dd(ii,kk+1)+ee(ii,kk+1)*(xw(ii,kk)-xw(ii,kk+2)))


       ! Eq 4.22 on Page 42 of Huynh (1991)
       fdot(ii,kk) = medan( fdot(ii,kk), qpl, qpr )

       ! Eq 4.11 on Page 38 of Huynh (1991)
       ttlmt = minmod(qpl, qpr)

       ! Eq 4.19a on Page 41 of Huynh (1991)
       tmin = min(0._r8,3*ss(ii,kk),1.5_r8*tt,ttlmt)

       ! Eq 4.19b on Page 41 of Huynh (1991)
       tmax = max(0._r8,3*ss(ii,kk),1.5_r8*tt,ttlmt)

       ! Eq 4.20 on Page 41 of Huynh (1991)
       fdot(ii,kk) = fdot(ii,kk) + minmod( tmin-fdot(ii,kk), tmax-fdot(ii,kk))

    end do
    end do

  end subroutine cfdotmc_pro

  !##############################################################################
  ! The minmod function 
  !##############################################################################
  real(r8) function minmod(aa,bb)

    implicit none
    real(r8),intent(in) :: aa, bb

    minmod = 0.5_r8*(sign(1._r8,aa) + sign(1._r8,bb))*min(abs(aa),abs(bb))

  end function minmod

  !##############################################################################
  ! The medan function
  !##############################################################################
  real(r8) function medan(aa,bb,cc)

    implicit none
    real(r8),intent(in) :: aa, bb, cc

    medan = aa + minmod(bb-aa,cc-aa)

  end function medan

end module mo_spitfire_transport
