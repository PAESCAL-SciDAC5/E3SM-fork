  subroutine clubb_stats_output(pcols,lchnk)

    use cam_history,      only: outfld
    use clubb_api_module, only: l_stats, stats_tout, stats_zt, stats_sfc, &
                                stats_zm, stats_rad_zt, stats_rad_zm, l_output_rad_files

    ! The following variables module variables in clubb_intr.F90:
    !  out_zt, out_zm, out_sfc, out_radzt, out_radzm.

    integer, intent(in) :: pcols, lchnk

    integer :: i
    character(len=200) :: temp1, sub                  ! Strings needed for CLUBB output

    !-------- 
    if (l_stats) then

      do i=1,stats_zt%num_output_fields

         temp1 = trim(stats_zt%file%var(i)%name)
         sub   = temp1
         if (len(temp1) .gt. 16) sub = temp1(1:16)

         call outfld(trim(sub), out_zt(:,:,i), pcols, lchnk )
      enddo

      do i=1,stats_zm%num_output_fields

         temp1 = trim(stats_zm%file%var(i)%name)
         sub   = temp1
         if (len(temp1) .gt. 16) sub = temp1(1:16)

         call outfld(trim(sub),out_zm(:,:,i), pcols, lchnk)
      enddo

      if (l_output_rad_files) then
         do i=1,stats_rad_zt%num_output_fields
            call outfld(trim(stats_rad_zt%file%var(i)%name), out_radzt(:,:,i), pcols, lchnk)
         enddo

         do i=1,stats_rad_zm%num_output_fields
            call outfld(trim(stats_rad_zm%file%var(i)%name), out_radzm(:,:,i), pcols, lchnk)
         enddo
      endif

      do i=1,stats_sfc%num_output_fields
         call outfld(trim(stats_sfc%file%var(i)%name), out_sfc(:,:,i), pcols, lchnk)
      enddo

    end if
    !-------- 

  end subroutine clubb_stats_output
