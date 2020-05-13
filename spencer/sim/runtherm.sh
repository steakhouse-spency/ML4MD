#!/bin/bash

# must set path of src

# HPC
 wd="/staging/an2/sportega/ML4MD/spencer"
 lmpexec="/home/rcf-40/sportega/disk/lammps/lammps/src/lmp_foo"

# lammps mpi execution command
lmprun="srun --mpi=pmi2 $lmpexec"
material="C"
datafiles=("C_29_4_11")

cd $wd
#datafiles=($(ls data_file))

for file in ${datafiles[*]}; do
	echo $file

	# remove file extension for label
	label="${file%.*}"

	# get nanotube material
	material=$(echo $label |  cut -d'_' -f1)

	# mkdir and cd into new dir
	mkdir -p sim/$label/therm && cd "$_" 
	
	infile="${label}-therm.in"
	slfile="${label}-therm.slurm"

	# copy input file and slurm script to new dir
	cp $wd/input_file/$material/therm.in $infile
	cp $wd/slurm/$material/therm.slurm $slfile
	
	# copy other files
	cp $wd/filled_tube/$label.pdb .
	cp $wd/data_file/$label.data .


	# substitute values for templat slurm/in files
	sed -i "s/variable label string/variable label string $label/1" $infile
	sed -i "s/<label>/$label/1" $slfile
	echo "cd $wd/sim/$label/therm" >> $slfile
	echo "${lmprun} < ${infile}" >> $slfile

	# submit to slurm
    sbatch $slfile
	sleep 2

	# reset to wd
	cd $wd

done
