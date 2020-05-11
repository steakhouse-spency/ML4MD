#!/bin/bash
#SBATCH --job-name=f
#SBATCH --export=none
#SBATCH --ntasks=16 --constraint=xeon
#SBATCH --nodes=-1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sportega@usc.edu
#SBATCH --output=C_29_4_11-270_1_f.out

source /usr/usc/intel/16.0/setup.sh
source /usr/usc/openmpi/1.8.8/setup.sh.intel
cd /staging/an2/sportega/ML4MD/spencer/sim/C_29_4_11/270/run_1
srun --mpi=pmi2 /home/rcf-40/sportega/disk/lammps/lammps/src/lmp_foo < C_29_4_11-270_1_f.in
