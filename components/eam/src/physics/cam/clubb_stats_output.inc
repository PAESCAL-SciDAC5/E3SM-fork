
   !  Output CLUBB history here
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

   endif
