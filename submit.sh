#!/bin/bash
#SBATCH --job-name=oceananigans
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --nodes 1
#SBATCH --mem 15G
#SBATCH --partition gpu_devel
#SBATCH --gpus=1

module reset
module load miniconda
module load Julia/1.11.3-linux-x86_64

###---------------------------------------------------------
### if true, this will instantiste the project
### set to true only if the project is not instantiated yet
### this will only need to be done once
###---------------------------------------------------------

INSTANTIATE=false

###---------------------------------------------------------
### path to simulation you want to run
###---------------------------------------------------------

SIMULATION=/home/${USER}/scratch_pi_mt477/${USER}/cross-section-model/simulations/simulation.jl

###---------------------------------------------------------
### this is path where where Project.toml is
### note: do not put Project.toml at end of the path
###---------------------------------------------------------

PROJECT=/home/${USER}/scratch_pi_mt477/${USER}/cross-section-model/

###---------------------------------------------------------
### this is where all downloaded file and packages will live
### will make scratch directory if does not already exist
###---------------------------------------------------------

DEPOT_PATH=/home/${USER}/scratch_pi_mt477/${USER}/julia-depot/

# this just makes DEPOT_PATH directory if non-existent
mkdir -p ${DEPOT_PATH}

export JULIA_DEPOT_PATH=${DEPOT_PATH}

###-------------------------------------------
### instantiates packages
### should only have to run this once
###-------------------------------------------

if ${INSTANTIATE}; then
    julia --project="${PROJECT}" -e "using Pkg; Pkg.instantiate()"
fi

wait

###-------------------------------------------
### runs the actual simulation
###-------------------------------------------

julia --project=${PROJECT} ${SIMULATION}

