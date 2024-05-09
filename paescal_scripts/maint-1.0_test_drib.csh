#!/bin/csh
date
set echo verbose

# Note: below, we source an older version of the E3SM Unified environment because the "latest"
# does not work (in spring 2024)

source /global/common/software/e3sm/anaconda_envs/load_e3sm_unified_1.8.1_pm-cpu.csh

setenv PROJECT  m4359

set compile_model = 0   # 0 = No, >0 = Yes
set run_model     = 1   # 0 = No, >0 = Yes

####################################################################
# code
####################################################################
setenv CODEROOT /global/cfs/projectdirs/${PROJECT}/${USER}/dcape/code
setenv CCSMTAG  maint-1.0_dcape
setenv CCSMROOT ${CODEROOT}/${CCSMTAG}

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv COMPSET FC5AV1C-L
setenv MRES ne30
setenv RESOLUTION ${MRES}_${MRES}
setenv MACH      pm-cpu 
setenv PTMP      /pscratch/sd/h/${USER}/dcape/${CCSMTAG}

setenv ntasks 512 
setenv nthrds 1

setenv CASE     ${CCSMTAG}_${COMPSET}_dcape_off_drib_on

setenv CASEROOT  $PTMP/cases/$CASE
setenv RUNDIR    $PTMP/$CASE/run

####################################################################
# Compile model
####################################################################
if ($compile_model > 0 || $run_model > 0) then

   rm -rf $CASEROOT
   cd  $CCSMROOT/cime/scripts

###usage: create_newcase [-h] [-d] [-v] [-s] --case CASE --compset COMPSET --res
###                      RES [--mach MACH] [--compiler COMPILER] [--ninst NINST]
###                      [--mpilib MPILIB] [--project PROJECT]
###                      [--pecount PECOUNT] [--mach-dir MACH_DIR]
###                      [--user-mods-dir USER_MODS_DIR] [--user-compset]
###                      [--pesfile PESFILE] [--user-grid] [--gridfile GRIDFILE]
###                      [--srcroot SRCROOT] [--test] [--walltime WALLTIME]
###                      [-q QUEUE]


   ./create_newcase --case $CASEROOT --project $PROJECT --mach $MACH \
                    --res $RESOLUTION --compset $COMPSET

#====================================================================
# set up case
#====================================================================

   ###./create_newcase -list grids

   cd $CASEROOT

   ./xmlchange -file env_run.xml      -id RUNDIR     -val $RUNDIR
   ./xmlchange -file env_build.xml    -id EXEROOT    -val $PTMP/$CASE/bld/

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'

   ./case.setup

#====================================================================
# Compile 
#====================================================================
   cd $CASEROOT

#  ./xmlchange -file env_build.xml -id DEBUG       -val 'TRUE'
   ./xmlchange -file env_build.xml -id PIO_VERSION -val '1'

   if ($compile_model > 0) then # Build the model

      ./case.build

   else # mark as already built

      ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

   endif

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then
#------------------
## set environment
#------------------

cd $CASEROOT

./xmlchange  -file env_run.xml  -id  RUN_STARTDATE      -val '2009-10-01'
./xmlchange  -file env_run.xml  -id  RESUBMIT           -val '0'
./xmlchange  -file env_run.xml  -id  CONTINUE_RUN       -val 'FALSE'
./xmlchange  -file env_run.xml  -id  STOP_N             -val '5'
./xmlchange  -file env_run.xml  -id  STOP_OPTION        -val 'ndays'
./xmlchange  -file env_run.xml  -id  REST_N             -val '5'
./xmlchange  -file env_run.xml  -id  REST_OPTION        -val 'ndays'
./xmlchange  -file env_run.xml  -id  DOUT_S             -val 'FALSE'

./xmlchange  -file env_run.xml  -id  PIO_TYPENAME       -val 'netcdf'

./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val '2:00:00'
#./xmlchange JOB_QUEUE=short --force

cat <<EOF >! user_nl_cam
&camexp

cld_cpl_opt = 2  ! turn on dribbling

!.......................................................
! history files
!.......................................................
 nhtfrq          =  0,
 mfilt           =  1,
 avgflag_pertape = 'A',

 history_aero_optics = .false.
 history_amwg        = .true.
 history_aerosol     = .true.
 history_verbose     = .true.
!
!...................
! change init data
!...................
ncdata = '/global/cfs/projectdirs/m4359/huiwan/e3sm_input/v1_init/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.cam.i.0000-10-01-00000.nc'
/
EOF

cat <<EOB >! user_nl_clm
&clm_inparm
 finidat = '/global/cfs/projectdirs/m4359/huiwan/e3sm_input/v1_init/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.clm2.r.0000-10-01-00000.nc'
 check_finidat_fsurdat_consistency = .false.
/
EOB

./case.submit

endif
