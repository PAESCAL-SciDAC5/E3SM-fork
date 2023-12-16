

   !  Determine CLUBB time step and make it sub-step friendly
   !  For now we want CLUBB time step to be 5 min since that is
   !  what has been scientifically validated.  However, there are certain
   !  instances when a 5 min time step will not be possible (based on
   !  host model time step or on macro-micro sub-stepping

   dtime = clubb_timestep
   hdtime_core_rknd = real(hdtime, kind = core_rknd)

   !  Now check to see if dtime is greater than the host model
   !    (or sub stepped) time step.  If it is, then simply
   !    set it equal to the host (or sub step) time step.
   !    This section is mostly to deal with small host model
   !    time steps (or small sub-steps)

   if (dtime .gt. hdtime_core_rknd) then
     dtime = hdtime_core_rknd
   endif

   !  Now check to see if CLUBB time step divides evenly into
   !    the host model time step.  If not, force it to divide evenly.
   !    We also want it to be 5 minutes or less.  This section is
   !    mainly for host model time steps that are not evenly divisible
   !    by 5 minutes

   if (mod(hdtime_core_rknd,dtime) .ne. 0) then
     dtime = hdtime_core_rknd/2._core_rknd
     do while (dtime .gt. 300._core_rknd)
       dtime = dtime/2._core_rknd
     end do
   endif

   !  If resulting host model time step and CLUBB time step do not divide evenly
   !    into each other, have model throw a fit.

   if (mod(hdtime_core_rknd,dtime) .ne. 0) then
     call endrun('clubb_tend_cam:  CLUBB time step and HOST time step NOT compatible'//errmsg(__FILE__,__LINE__))
   endif

   !  determine number of timesteps CLUBB core should be advanced,
   !  host time step divided by CLUBB time step
   nadv = max(hdtime_core_rknd/dtime,1._core_rknd)
