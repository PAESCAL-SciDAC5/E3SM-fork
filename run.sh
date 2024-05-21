#!/bin/bash
#SBATCH --time=5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=cpu
#SBATCH --account=m4359

srun CondiDiag_examples/run_maint-2.0_cam_inout_RHNregularization.sh
