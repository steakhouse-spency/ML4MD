#!/bin/bash
#SBATCH --export=none
#SBATCH --ntasks=16 --constraint=xeon
#SBATCH --nodes=-1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sportega@usc.edu
#SBATCH --output=C_29_4_11-rebo-therm.out

source /usr/usc/intel/16.0/setup.sh
source /usr/usc/openmpi/1.8.8/setup.sh.intel

cd /staging/an2/sportega/ML4MD/spencer/sim/C_29_4_11-rebo/therm
srun --mpi=pmi2 /home/rcf-40/sportega/disk/lammps/lammps/src/lmp_foo < C_29_4_11-rebo_therm.in
