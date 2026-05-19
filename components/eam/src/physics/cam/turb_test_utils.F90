module turb_test_utils

  use shr_kind_mod,  only: r8=>shr_kind_r8

  implicit none

  public

contains

subroutine set_switches_for_txt_output( single_column, masterproc, l_turb_standalone, &! in
                                        l_first_step, macmic_it,                      &! in
                                        l_open_txt_output, l_clse_txt_output          )! out

   logical, intent(in) :: single_column
   logical, intent(in) :: masterproc
   logical, intent(in) :: l_turb_standalone
   logical, intent(in) :: l_first_step
   integer, intent(in) :: macmic_it

   logical, intent(out) :: l_open_txt_output
   logical, intent(out) :: l_clse_txt_output

   ! Txt output is handled only by the master MPI process in an SCM run.

   if (single_column.and.masterproc) then

      ! If the SCM is configured to be turb standalone, in principle
      ! we only need a single call of this subroutine (shoc_tend_eam)
      ! and a single txt file from the call, despite the fact
      ! that a 1-timestep EAM run will call this subroutine 3 times
      ! because of model initialization and finalization.

      if (l_turb_standalone) then

         l_open_txt_output = .true.   ! open a new file and write header info everytime we get here.
         l_clse_txt_output = .true.   ! close the file and free i/o unit in this subroutine after shoc is called.

      ! If the SCM is not confugured to be turb standalone,
      ! open a new file and write header info at the very beginning of the simulation,
      ! i.e., in the very first EAM time step during the first macmic substep;
      ! keep the file open during the entire simulation.
      else

         l_open_txt_output = l_first_step .and. (macmic_it.eq.1)
         l_clse_txt_output = .false.

      end if

   else
   ! In global simulations, we won't have txt output; set both switches to .false.
   ! In case an SCM simulation contains multiple grid columns for some reason,
   ! do not produce txt output from MPI processes other than the master.

      l_open_txt_output = .false.
      l_clse_txt_output = .false.

   end if

end subroutine set_switches_for_txt_output


subroutine txt_file_open_and_init( prefix_in, nstep_current, macmic_it, nvarnames, varnames, &! in
                                   txtout_unit                                               )! out

   use units, only: getunit

   character(len=*),intent(in) :: prefix_in
   integer,intent(in)          :: nstep_current
   integer,intent(in)          :: macmic_it
   integer,intent(in)          :: nvarnames
   character(len=*),intent(in) :: varnames(nvarnames)
   integer,intent(out)         :: txtout_unit

   ! Local variables

   character(len=256)  :: txt_output_prefix = ''
   character(len=256)  :: outfname
   integer :: iv

   ! Determine name of txt output file

   txt_output_prefix = 'turb_output'                                ! set a default
   if (len_trim(prefix_in) > 0) txt_output_prefix = trim(prefix_in) ! overwrite by nml input
   write(outfname, '(A,2(A,I0),A)') trim(txt_output_prefix), &
                                    '_nstep',nstep_current, '_macmicsub',macmic_it,'.txt'

   ! Open txt file and write a header

   txtout_unit = getunit()
   open(txtout_unit, file=trim(outfname), status='replace')
   write(txtout_unit,'(A10,A5,20A23)') (varnames(iv), iv=1,nvarnames) 
 
end subroutine txt_file_open_and_init


subroutine txt_write_one_column(txtout_unit,modeltime,nlev,nv,merged_array )

  integer, intent(in) :: txtout_unit
  real(r8),intent(in) :: modeltime
  integer, intent(in) :: nlev
  integer, intent(in) :: nv    ! # of physical quantities, not counting time and level index

  real(r8),intent(in) :: merged_array(:,:)

  integer :: kk, iv
  character(len=*),parameter :: fmt = '(F8.2,I5,20E23.14)'

  do kk = nlev, 1, -1
     write(txtout_unit,fmt) modeltime,&!
                            kk - 1,   &! 0 = TOM, pver - 1 = sfc
                            (merged_array(kk,iv),iv=1,nv)
  end do

end subroutine txt_write_one_column

end module turb_test_utils
