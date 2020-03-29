#!/bin/bash




# must set path of src

# HPC
# wd="/home/rcf-40/sportega/disk/ML4MD/spencer"
# lmpexec="/home/rcf-40/sportega/disk/lammps/lammps/src/lmp_foo"

# Local
wd="/home/sportega/Desktop/dr/spencer"
lmpexec="/home/rcf-proj/an2/sportega/lammps/lammps/src/lmp_foo"

# lammps mpi execution command
lmprun="srun --mpi=pmi2 $lmpexec"

material="C"

cd $wd
#datafiles=($(ls data_file))
 datafiles=("C_25_7_14")

for file in ${datafiles[*]}; do
	echo $file

	# remove file extension for label
	label="${file%.*}"

	# mkdir and cd into new dir
	mkdir -p sim/$label/therm && cd "$_" 

	# copy input file and slurm script to new dir
	cp $wd/slurm/$material/therm.slurm $label.slurm
	cp $wd/input_file/$material/therm.in $label.in

	# edit files
	sed -i "s/variable label string/variable label string $label/1" $label.in
	sed -i "s/<label>/$label/1" $label.slurm
	echo "cd $wd/sim/$label/therm" >> $label.slurm
	echo "${lmprun} < $label.in" >> $label.slurm

	# submit to slurm
	#sbatch $label.slurm
	#sleep 2

	# reset to wd
	cd $wd

done
