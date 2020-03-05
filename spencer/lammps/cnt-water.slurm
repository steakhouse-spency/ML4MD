#!/bin/bash
#SBATCH --export=none
#SBATCH --ntasks=16 --constraint=xeon
#SBATCH --nodes=-1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sportega@usc.edu


cd "$SLURM_SUBMIT_DIR"

source /usr/usc/intel/16.0/setup.sh
source /usr/usc/openmpi/1.8.8/setup.sh.intel 
srun --mpi=pmi2 /auto/rcf-proj2/jgc/cobenare/lammps/lammps-12Dec18/src/lmp_foo < in-sicnt-wat-prod-000-nm-10_3-0-a.in
