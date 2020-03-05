#!/bin/bash
#SBATCH --export=none
#SBATCH --ntasks=16 --constraint=xeon
#SBATCH --nodes=-1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sportega@usc.edu


lmp_src=~/disk/lammps/lammps/src/ 
wd=~/disk/ML4MD/spencer/lammps/
cd $wd
source /usr/usc/intel/16.0/setup.sh
source /usr/usc/openmpi/1.8.8/setup.sh.intel 

srun --mpi=pmi2 $lmp_src/lmp_foo < cnt_water.in
