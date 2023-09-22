#!/bin/bash -fe

# E3SM Water Cycle v2 run_e3sm script template.
#
# Inspired by v1 run_e3sm script as well as SCREAM group simplified run script.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

# TO DO:
# - custom pelayout

main() {

# For debugging, uncomment line below
#set -x

# --- Configuration flags ----

# Machine and project
readonly E3SM_version="1.0"
readonly MACHINE=pm-cpu
readonly PROJECT="m4359"      #need to change
readonly COMPILER=intel
readonly email_address="jianfeng.li@pnnl.gov"             #need to change

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# atmospheric coupling option
readonly cflx_cpl_opt=42
readonly nverlvl=72

# Simulation
readonly COMPSET="F20TRC5-CMIP6"
readonly RESOLUTION="ne30_ne30"
# BEFORE RUNNING : CHANGE the following CASE_NAME to desired value
readonly date_string=`date +"%Y%m%d%H%M"`
readonly CASE_NAME="E3SMv${E3SM_version}.${COMPSET}.${RESOLUTION}.${date_string}_cflx_${cflx_cpl_opt}_L${nverlvl}"
# If this is part of a simulation campaign, ask your group lead about using a case_group label
# readonly CASE_GROUP=""

# Code and compilation
readonly BRANCH="huiwanpnnl/gmd_2020_330_forc+cflx_202305_pm_cpu_gnu"
readonly REPO_name="PAESCAL-SciDAC5/E3SM-fork"
readonly DEBUG_COMPILE=false
readonly CCSMTAG="gmd_2020_330_cflx_202305"

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="2009-10-01"
readonly START_SECOND="0"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFDIR="/global/cscratch1/sd/forsyth/E3SMv2/v2.LR.piControl/init"
readonly RUN_REFCASE="20210625.v2rc3c-GWD.piControl.ne30pg2_EC30to60E2r2.chrysalis"
readonly RUN_REFDATE="1001-01-01"   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
CODE_ROOT="/global/cfs/projectdirs/m4359/jli628/E3SM_v1.0_SciDAC/"     #need to change
if [ "${do_fetch_code,,}" != "false" ]; then
	CODE_ROOT=${CODE_ROOT}${CCSMTAG}     #need to change
fi
echo $CODE_ROOT
readonly CASE_ROOT="${SCRATCH}/E3SMv${E3SM_version}/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Production simulation
readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
readonly CASE_RUN_DIR=${CASE_ROOT}/run
#PELAYOUT doesn't work on pm-cpu
#readonly PELAYOUT="X2"   #Allowed options are  ('S','M','L','X1','X2','[0-9]x[0-9]','[0-9]').
readonly WALLTIME="12:00:00"
readonly QUEUE="regular"
readonly STOP_OPTION="nmonths"
readonly STOP_N="15"
readonly REST_OPTION="nmonths"
readonly REST_N="15"
readonly RESUBMIT="0"
readonly DO_SHORT_TERM_ARCHIVING=false

# Setup processor layout
readonly nnodes_atm=8
readonly nnodes_ocn=8     #be careful, not all parallel configurations are available for maps-seaice 
readonly nthreads=1
readonly mpi_tasks_per_node=128
readonly ntasks_atm=$(expr ${nnodes_atm} \* ${mpi_tasks_per_node})
readonly ntasks_ocn=$(expr ${nnodes_ocn} \* ${mpi_tasks_per_node})
readonly total_tasks_per_node=$(expr ${mpi_tasks_per_node} \* ${nthreads})
if [ ${ntasks_ocn} -ne ${ntasks_atm} ]; then
    nnodes=$(expr ${nnodes_atm} + ${nnodes_ocn})
else
    nnodes=${nnodes_atm}
fi
pelayout=${nnodes}x${mpi_tasks_per_node}x${nthreads}

stridc=1
npcpl=$(expr ${ntasks_atm} \/ ${stridc})
echo "npcpl=${npcpl}"

# Coupler history
readonly HIST_OPTION="nmonths"
readonly HIST_N="1"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Now, do the work ---

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

# Setup
case_setup

# Build
case_build

# Configure runtime options
runtime_options

# Copy script into case_script directory for provenance
copy_script

# Submit
case_submit

# All done
echo $'\n----- All done -----\n'

}

# =======================
# Custom user_nl settings
# =======================

user_nl() {

cat << EOF >> user_nl_cam
&camexp
!
cflx_cpl_opt = ${cflx_cpl_opt}
!

dtime           = 1800

inithist = 'MONTHLY'
inithist_all = .true.
!
 rad_diag_1 = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2',
             'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4',
             'N:CFC11:CFC11', 'N:CFC12:CFC12',
 docosp    = .true.,
 cosp_lite = .true.,
 cosp_ncolumns        = 10
 cosp_nradsteps       = 3
 cosp_lmisr_sim       = .true.
 cosp_lisccp_sim      = .true.
 cosp_lmodis_sim      = .true.
 cosp_llidar_sim      = .true.
 history_amwg         = .true.
 history_aero_optics  = .true.
 history_aerosol      = .true.
 history_clubb        = .true.
 history_budget       = .true.
 history_verbose      = .true.
 hist_hetfrz_classnuc = .true.
 do_aerocom_ind3 = .true.

 Nudge_Model          = .True.
 Nudge_Path           = '/pscratch/sd/k/kaizhang/ndata/eraint_ne30L72/'
 Nudge_File_Template  = 'interim_se_%y%m%d00.nc'
 Nudge_Times_Per_Day  = 4  !! nudging input data frequency
 Model_Times_Per_Day  = 48 !! should not be larger than 48 if dtime = 1800s
 Nudge_Uprof          = 2
 Nudge_Ucoef          = 1.
 Nudge_Vprof          = 2
 Nudge_Vcoef          = 1.
 Nudge_Tprof          = 0
 Nudge_Tcoef          = 0.
 Nudge_Qprof          = 0
 Nudge_Qcoef          = 0.
 Nudge_PSprof         = 0
 Nudge_PScoef         = 0.
 Nudge_Beg_Year       = 0001
 Nudge_Beg_Month      = 1
 Nudge_Beg_Day        = 1
 Nudge_End_Year       = 9999
 Nudge_End_Month      = 1
 Nudge_End_Day        = 1
 Nudge_Vwin_Lindex    = 5.
 Nudge_Vwin_Hindex    = 58.
 Nudge_Vwin_Ldelta    = 0.1
 Nudge_Vwin_Hdelta    = 0.1
 Nudge_Vwin_lo        = 0.
 Nudge_Vwin_hi        = 1.
 Nudge_Method         = 'Linear'
 Nudge_Loc_PhysOut    = .True.
 Nudge_Tau            = 6.        !! relaxation time scale, unit: 6h
 Nudge_CurrentStep    = .False.
 Nudge_File_Ntime     = 4

!.......................................................
! emissions
!.......................................................
ext_frc_specifier              = 'SO2         -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_so2_elev_1850-2014_c180205.nc',
        'SOAG        -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_soag_elev_1850-2014_c180205.nc',
        'bc_a4       -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
        'num_a1      -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
        'num_a2      -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
        'num_a4      -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
        'pom_a4      -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
        'so4_a1      -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
        'so4_a2      -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'
ext_frc_type           = 'CYCLICAL'
ext_frc_cycle_yr       = 2010
srf_emis_specifier             = 'DMS       -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/DMSflux.2010.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20190220.nc',
        'SO2       -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_so2_surf_1850-2014_c180205.nc',
        'bc_a4     -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_bc_a4_surf_1850-2014_c180205.nc',
        'num_a1    -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_num_a1_surf_1850-2014_c180205.nc',
        'num_a2    -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_num_a2_surf_1850-2014_c180205.nc',
        'num_a4    -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_num_a4_surf_1850-2014_c180205.nc',
        'pom_a4    -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_pom_a4_surf_1850-2014_c180205.nc',
        'so4_a1    -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_so4_a1_surf_1850-2014_c180205.nc',
        'so4_a2    -> /global/cfs/projectdirs/m4359/jli628/emissions_2010/cmip6_mam4_so4_a2_surf_1850-2014_c180205.nc'
srf_emis_type          = 'CYCLICAL'
srf_emis_cycle_yr      = 2010

 nhtfrq          =  0,
 mfilt           =  1,
 avgflag_pertape = 'A',

 fincl1  =   'PS',
            'CAPE',
            'CIN',
            'WSUBI',       !! characteristic vertical velocity for ice nucleation
            'WLARGE',      !! large scale vertical velocity
            'WSIG',        !! subgrid scale vertical velocity
            'WSUBI2',      !! Mean updraft, with stddev=f(TKE)
            'DCST',        !! dcs
            'RHICE',       !! RH ice
            'TOT_CLD_VISTAU',
            'TOT_ICLD_VISTAU',
            'LIQ_ICLD_VISTAU',
            'ICE_ICLD_VISTAU',
            'SNOW_ICLD_VISTAU',
            'AODVIS',  !! AeroCom indirect diagnostics
            'angstrm',
            'cod',
            'cdr',
            'cdnc',
            'cdnum',
            'icnum',
            'clt',
            'lcc',
            'lwp',
            'iwp',
            'icc',
            'icnc',
            'icr',
            'LHFLX',
            'SHFLX',
            'OMEGA500',
            'rh700',
            'colrv',
            'ccn',
            'ccn.1bl',
            'ccn.3bl',
            'ptop',
            'ttop',
            'rwp',
            'lwp2',
            'iwp2',
            'autoconv',
            'accretn',
            'FSUTOA_d1',
            'FSUTOAC_d1',
            'FSUTOA',
            'FSUTOAC',
            'FLUTC',
            'FLUT',
            'PRECC',
            'PRECL',
            'PRECT',
            'TH7001000',

!...................
! change init data
!...................
ncdata = '/global/cfs/projectdirs/m4359/jli628/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.cam.i.0000-10-01-00000.nc'
/
EOF

cat << EOF >> user_nl_clm
&clm_inparm
 finidat = '/global/cfs/projectdirs/m4359/jli628/constance_FC5AV1C-04P2_ne30_ne30_E3SMv1.0_CLIMO_PD_TUNE.clm2.r.0000-10-01-00000.nc'
 check_finidat_fsurdat_consistency = .false.
/
EOF

}

patch_mpas_streams() {

echo

}

######################################################
### Most users won't need to change anything below ###
######################################################

#-----------------------------------------------------
fetch_code() {

    if [ "${do_fetch_code,,}" != "true" ]; then
        echo $'\n----- Skipping fetch_code -----\n'
        return
    fi

    echo $'\n----- Starting fetch_code -----\n'
    local path=${CODE_ROOT}
    local repo=${REPO_name}

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    git clone -b $BRANCH --recursive git@github.com:${repo}.git .

    # Setup git hooks
    #rm -rf .git/hooks
    #git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
    #git config commit.template .git/hooks/commit.template

    # Check out desired branch
    #git checkout ${BRANCH}

    # Custom addition
    if [ "${CHERRY}" != "" ]; then
        echo ----- WARNING: adding git cherry-pick -----
        for commit in "${CHERRY[@]}"
        do
            echo ${commit}
            git cherry-pick ${commit}
        done
        echo -------------------------------------------
    fi

    # Bring in all submodule components
    #git submodule update --init --recursive

    popd
    cp ${CODE_ROOT}/paescal_scripts/addtional_files/ODEMod.F90 ${CODE_ROOT}/components/clm/src/external_models/sbetr/src/betr/betr_math/
}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
        echo $'\n----- Skipping create_newcase -----\n'
        return
    fi

    echo $'\n----- Starting create_newcase -----\n'

	if [[ -z "$CASE_GROUP" ]]; then
		${CODE_ROOT}/cime/scripts/create_newcase \
			--case ${CASE_NAME} \
			--output-root ${CASE_ROOT} \
			--script-root ${CASE_SCRIPTS_DIR} \
			--handle-preexisting-dirs u \
			--compset ${COMPSET} \
			--res ${RESOLUTION} \
			--machine ${MACHINE} \
			--compiler ${COMPILER} \
			--project ${PROJECT}
			#--walltime ${WALLTIME} \
			#--queue ${QUEUE} \
			#--pecount ${PELAYOUT}
	else
		${CODE_ROOT}/cime/scripts/create_newcase \
			--case ${CASE_NAME} \
			--case-group ${CASE_GROUP} \
			--output-root ${CASE_ROOT} \
			--script-root ${CASE_SCRIPTS_DIR} \
			--handle-preexisting-dirs u \
			--compset ${COMPSET} \
			--res ${RESOLUTION} \
			--machine ${MACHINE} \
			--compiler ${COMPILER} \
			--project ${PROJECT}
			#--walltime ${WALLTIME} \
			#--queue ${QUEUE} \
			#--pecount ${PELAYOUT}
	fi
	

    if [ $? != 0 ]; then
      echo $'\nNote: if create_newcase failed because sub-directory already exists:'
      echo $'  * delete old case_script sub-directory'
      echo $'  * or set do_newcase=false\n'
      exit 35
    fi

}

#-----------------------------------------------------
case_setup() {

    if [ "${do_case_setup,,}" != "true" ]; then
        echo $'\n----- Skipping case_setup -----\n'
        return
    fi

    echo $'\n----- Starting case_setup -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Setup some CIME directories
    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}

    # Short term archiving
    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING^^}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

    # Set processor layout
    if [ ${ntasks_ocn} -ne ${ntasks_atm} ]; then
        ./xmlchange NTASKS=${ntasks_ocn}
        ./xmlchange NTASKS_ATM=${ntasks_atm}
        ./xmlchange ROOTPE_ATM=${ntasks_ocn}
    else
        ./xmlchange NTASKS=${ntasks_atm}
    fi
    ./xmlchange MAX_MPITASKS_PER_NODE=${mpi_tasks_per_node}
    ./xmlchange NTASKS_CPL=${npcpl} # change tasks in CPL if using strid
    ./xmlchange PSTRID_CPL=${stridc}
    ./xmlchange NTHRDS=1 # set all to 1 first
    ./xmlchange NTHRDS_ATM=${nthreads}
    ./xmlchange MAX_TASKS_PER_NODE=${total_tasks_per_node}

    #./xmlchange EPS_AGRID=1e-9

    #./xmlchange ATM_NCPL=96
    
    ./xmlchange --id CAM_CONFIG_OPTS --val='-phys cam5 -clubb_sgs -microphys mg2 -chem linoz_mam4_resus_mom_soag -rain_evap_to_coarse_aero -bc_dep_to_snow_updates -nlev '${nverlvl}

    # Build with COSP, except for a data atmosphere (datm)
    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Custom user_nl
    user_nl

    # Finally, run CIME case.setup
    ./case.setup --reset

    popd
}

#-----------------------------------------------------
case_build() {

    pushd ${CASE_SCRIPTS_DIR}

    # do_case_build = false
    if [ "${do_case_build,,}" != "true" ]; then

        echo $'\n----- case_build -----\n'

        if [ "${OLD_EXECUTABLE}" == "" ]; then
            # Ues previously built executable, make sure it exists
            if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
                echo 'Skipping build because $do_case_build = '${do_case_build}
            else
                echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
                exit 297
            fi
        else
            # If absolute pathname exists and is executable, reuse pre-exiting executable
            if [ -x ${OLD_EXECUTABLE} ]; then
                echo 'Using $OLD_EXECUTABLE = '${OLD_EXECUTABLE}
                cp -fp ${OLD_EXECUTABLE} ${CASE_BUILD_DIR}/
            else
                echo 'ERROR: $OLD_EXECUTABLE = '$OLD_EXECUTABLE' does not exist or is not an executable file.'
                exit 297
            fi
        fi
        echo 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
        ./xmlchange BUILD_COMPLETE=TRUE

    # do_case_build = true
    else

        echo $'\n----- Starting case_build -----\n'

        # Turn on debug compilation option if requested
        if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
            ./xmlchange DEBUG=${DEBUG_COMPILE^^}
        fi

	./xmlchange -file env_build.xml -id PIO_VERSION -val '1'

        # Run CIME case.build
        ./case.build

        # Some user_nl settings won't be updated to *_in files under the run directory
        # Call preview_namelists to make sure *_in and user_nl files are consistent.
        ./preview_namelists

    fi

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Set simulation start date
    ./xmlchange RUN_STARTDATE=${START_DATE}
    ./xmlchange START_TOD=${START_SECOND}

    # Segment length
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

    # Restart frequency
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

    # Coupler history
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

    # Coupler budgets (always on)
    ./xmlchange BUDGETS=TRUE

    # Set resubmissions
    if (( RESUBMIT > 0 )); then
        ./xmlchange RESUBMIT=${RESUBMIT}
    fi

    # Run type
    # Start from default of user-specified initial conditions
    if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
        ./xmlchange RUN_TYPE="startup"
        ./xmlchange CONTINUE_RUN="FALSE"

    # Continue existing run
    elif [ "${MODEL_START_TYPE,,}" == "continue" ]; then
        ./xmlchange CONTINUE_RUN="TRUE"

    elif [ "${MODEL_START_TYPE,,}" == "branch" ] || [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then
        ./xmlchange RUN_TYPE=${MODEL_START_TYPE,,}
        ./xmlchange GET_REFCASE=${GET_REFCASE}
	./xmlchange RUN_REFDIR=${RUN_REFDIR}
        ./xmlchange RUN_REFCASE=${RUN_REFCASE}
        ./xmlchange RUN_REFDATE=${RUN_REFDATE}
        echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE}
	echo '$RUN_REFDIR = '${RUN_REFDIR}
	echo '$RUN_REFCASE = '${RUN_REFCASE}
	echo '$RUN_REFDATE = '${START_DATE}

    else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    # Patch mpas streams files
    patch_mpas_streams

    popd
}

#-----------------------------------------------------
case_submit() {

    if [ "${do_case_submit,,}" != "true" ]; then
        echo $'\n----- Skipping case_submit -----\n'
        return
    fi

    echo $'\n----- Starting case_submit -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    ./xmlchange JOB_QUEUE=${QUEUE}
    ./xmlchange JOB_WALLCLOCK_TIME=${WALLTIME}
    # Run CIME case.submit
    ./case.submit --batch-args="--mail-type=ALL --mail-user=${email_address}"

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=`basename $0`
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp ${this_script_name} ${script_provenance_dir}/${script_provenance_name}

}

#-----------------------------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Now, actually run the script
#-----------------------------------------------------
main
