#!/bin/bash -fe

# E3SM Coupled Model Group run_e3sm script template.
#
# Bash coding style inspired by:
# https://web.archive.org/web/20200620202413/http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

main() {

# For debugging, uncomment line below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE="dane"
readonly PROJECT="paescal"

# Simulation
readonly COMPSET="F2010"
readonly RESOLUTION="ne30pg2_r05_IcoswISC30E3r5"
# readonly RESOLUTION="ne30pg2_r05_oECv3"
# BEFORE RUNNING : CHANGE the following CASE_NAME to desired value
readonly CASE_NAME="EAMv3_test_canopyfluxes_reducedpelayout_1year_sfc_inout_"$RESOLUTION
# readonly COMPSET="GMPAS-JRA1p4"
# readonly RESOLUTION="TL319_WC14to60E2r3"
# BEFORE RUNNING : CHANGE the following CASE_NAME to desired value
# readonly CASE_NAME="MPASv3_1yeardata_reducedtimestep2_sfc_inout_"$RESOLUTION
# If this is part of a simulation campaign, ask your group lead about using a case_group label
# readonly CASE_GROUP=""

# Code and compilation
# BEFORE RUNNING: CHANGE CHECKOUT to date string like 20240301
readonly CHECKOUT="20260409"
readonly BRANCH="ELM-flux-analysis"
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="2010-01-01"

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
#readonly RUN_REFDIR=""
#readonly RUN_REFCASE=""
#readonly RUN_REFDATE=""   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
readonly CASE_ROOT="/p/lustre2/dong9/e3sm_scratch/dane/cases/${CASE_NAME}"
readonly CODE_ROOT="/p/lustre2/dong9/E3SM-maint-3.0"


# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays',
#               'M_1x10_ndays', 'M2_1x10_ndays', 'M80_1x10_ndays', 'L_1x10_ndays'
#               * can replace XS, M, etc. with custom-XY with XY being the node count
#  or 'production' for full simulation
readonly run='production'
if [ "${run}" != "production" ]; then
  echo "setting up Short test simulations: ${run}"
  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}

  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  readonly PELAYOUT="custom-8"
  readonly WALLTIME="1:00:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false

else

  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="custom-8"
  readonly WALLTIME="3:00:00"
  readonly STOP_OPTION="nmonths"
  readonly STOP_N="15"
  readonly REST_OPTION="nmonths"
  readonly REST_N="4"
  readonly RESUBMIT="0"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# --- Now, do the work ---

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

# Custom PE layout
custom_pelayout

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

cat << EOF >> user_nl_eam
 cosp_lite = .true.

 empty_htapes = .true.

 ! -- chemUCI settings ------------------
 history_chemdyg_summary = .true.
 history_gaschmbudget_2D = .false.
 history_gaschmbudget_2D_levels = .false.
 history_gaschmbudget_num = 6 !! no impact if  history_gaschmbudget_2D = .false.

 ! -- MAM5 settings ------------------
 is_output_interactive_volc = .true.
EOF

cat << EOF >> user_nl_elm
hist_empty_htapes = .true.
hist_dov2xy = .true.,.true.
 hist_fincl1 = 'SNOWDP','COL_FIRE_CLOSS','NPOOL','PPOOL','TOTPRODC'
 hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA',
               'UBOT', 'VBOT', 'RHOBOT', 'UGUST', 'THM', 'TG', 'QBOT', 'QG', 'THBOT', 'FORC_HGT_U', 'Z0MG', 'Z0HG', 'Z0QG', 'THV', 'ZII',
               'DISPLA', 'ELAI', 'ESAI', 'Z0MG', 'Z0MV', 'DLEAF', 'HTOP', 'PBOT', 'TV', 'FDRY', 'LAISUN', 'LAISHA', 'RSSUN', 'RSSHA', 
               'RH', 'FWET', 'QVEGT', 'H2OCAN', 'FRAC_VEG_NOSNO', 'SNOW_DEPTH', 'SoilBeta', 'SNO_T', 'TSOI', 'FH2OSFC', 'TH2OSFC' 'SABV', 'EMV', 'EMG', 'QFLX_EVAP_VEG',
               'C3PSN', 'SLATOP', 'LEAFCN', 'FLNR', 'FNITR', 'I_VC', 'S_VC', 'LEAFN', 'LEAFP', 'T10', 'LAISUN_Z', 'LAISHA_Z',
               'NRAD', 'VCMAXCINTSUN', 'VCMAXCINTSHA', 'TLAI_Z', 'PARSUN_Z', 'PARSHA_Z',
               'C3FLAG', 'KC', 'KO', 'CP', 'QE', 'THETA_CJ', 'BBB', 'MBB', 'GB_MOL', 'GS_MOL', 'AC', 'AJ', 'AP', 'AG', 'AN', 'RH_LEAF',
               'VCMAX_Z', 'TPU_Z', 'KP_Z', 'DAYL', 'MAX_DAYL', 'SNL', 
               'watsat', 'SOILICE', 'EFF_POROSITY', 'SOILLIQ', 'H2OSOI_LIQVOL', 'SMPSO', 'SMPSC', 'TC_STRESS', 'sucsat',
               'bsw', 'ROOTFR', 'ROOTR', 'BTRAN', 'BTRAN2', 'RRESIS', 'Z0HV', 'Z0QV'
 hist_mfilt = 1,30
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','I'

 finidat = '/p/lustre2/dong9/e3sm_inputdata/v3.LR.amip_0101.elm.r.2010-01-01-00000.nc'
 fsurdat = '/p/vast1/e3sm/ccsm3data/inputdata/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c200609_with_TOP.nc'
 check_finidat_pct_consistency = .false.
EOF

}

#  hist_fincl3 = 'C3PSN', 'SLATOP', 'LEAFCN', 'FLNR', 'FNITR', 'I_VC', 'S_VC', 'LEAFN', 'LEAFP', 'T10'

#  hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA',
#                'UBOT', 'VBOT', 'RHOBOT', 'UGUST', 'THM', 'TG', 'QBOT', 'QG', 'THBOT', 'FORC_HGT_U', 'Z0MG', 'Z0HG', 'Z0QG', 'THV', 'ZII',
#                'DISPLA', 'ELAI', 'ESAI', 'Z0MG', 'Z0MV', 'DLEAF', 'HTOP', 'PBOT', 'TV', 'FDRY', 'LAISUN', 'LAISHA', 'RSSUN', 'RSSHA', 
#                'RH', 'FWET', 'QVEGT', 'H2OCAN', 'FRAC_VEG_NOSNO', 'SNOW_DEPTH', 'SoilBeta', 'TSOI', 'FH2OSFC', 'SABV', 'EMV', 'EMG', 'QFLX_EVAP_VEG'
#  fsurdat = '/p/vast1/e3sm/ccsm3data/inputdata/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c200609_with_TOP.nc'
#  finidat = '/p/lustre2/dong9/e3sm_inputdata/v3.LR.amip_0101.elm.r.2010-01-01-00000.nc'
#  check_finidat_pct_consistency = .false.

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
    local repo=e3sm

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    # This will put repository, with all code
    git clone git@github.com:E3SM-Project/${repo}.git .

    # Setup git hooks
    rm -rf .git/hooks
    git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
    git config commit.template .git/hooks/commit.template

    # Check out desired branch
    git checkout ${BRANCH}

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
    git submodule update --init --recursive

    popd
}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
        echo $'\n----- Skipping create_newcase -----\n'
        return
    fi

    echo $'\n----- Starting create_newcase -----\n'


    if [[ ${PELAYOUT} == custom-* ]];
    then
        layout="M" # temporary placeholder for create_newcase
    else
        layout=${PELAYOUT}
    fi

	if [[ -z "$CASE_GROUP" ]]; then
		${CODE_ROOT}/cime/scripts/create_newcase \
			--case ${CASE_NAME} \
			--output-root ${CASE_ROOT} \
			--script-root ${CASE_SCRIPTS_DIR} \
			--handle-preexisting-dirs u \
			--compset ${COMPSET} \
			--res ${RESOLUTION} \
			--machine ${MACHINE} \
			--project ${PROJECT} \
			--walltime ${WALLTIME} \
			--pecount ${layout}
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
			--project ${PROJECT} \
			--walltime ${WALLTIME} \
			--pecount ${layout}
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
custom_pelayout() {

if [[ ${PELAYOUT} == custom-* ]];
then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    # Number of cores per node (machine specific)
    if [ "${MACHINE}" == "pm-cpu" ]; then
        ncore=128
    elif [ "${MACHINE}" == "chrysalis" ]; then
        ncore=64
    elif [ "${MACHINE}" == "compy" ]; then
        ncore=40
    elif [ "${MACHINE}" == "anvil" ]; then
        ncore=36
    elif [ "${MACHINE}" == "dane" ]; then
        ncore=112
    else
        echo 'ERROR: MACHINE = '${MACHINE}' is not supported for custom PE layout.' 
        exit 400
    fi

    # Extract number of nodes
    tmp=($(echo ${PELAYOUT} | tr "-" " "))
    nnodes=${tmp[1]}

    # Customize
    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=$(( $nnodes * $ncore ))
    ./xmlchange NTHRDS=1
    ./xmlchange LND_NTASKS=560
    ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
    ./xmlchange MAX_TASKS_PER_NODE=$ncore
    popd

fi

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

        # Run CIME case.build
        ./case.build

        # Some user_nl settings won't be updated to *_in files under the run directory
        # Call preview_namelists to make sure *_in and user_nl files are consistent.
	 echo $'\n----- Preview namelists -----\n'
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

    # Run CIME case.submit
    ./case.submit

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=$( basename -- "$0"; )
    local this_script_dir=$( dirname -- "$0"; )
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp "${this_script_dir}/${this_script_name}" ${script_provenance_dir}/${script_provenance_name}

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
