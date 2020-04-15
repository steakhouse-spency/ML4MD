#!/bin/bash




# must set path of src

# HPC
 wd="/home/rcf-40/sportega/disk/ml4md/spencer"
 lmpexec="/home/rcf-40/sportega/disk/lammps/lammps/src/lmp_foo"

# Local
#wd="/home/sportega/Desktop/dr/spencer"
#lmpexec="/home/rcf-proj/an2/sportega/lammps/lammps/src/lmp_foo"

# lammps mpi execution command
lmprun="srun --mpi=pmi2 $lmpexec"

material="C"

cd $wd
#datafiles=($(ls data_file))
 #datafiles=("C_25_7_14-a" "C_25_7_14-b")
 datafiles=("C_25_7_14-b")

for file in ${datafiles[*]}; do
	echo $file

	# remove file extension for label
	label="${file%.*}"

	# mkdir and cd into new dir
	mkdir -p sim/$label/therm && cd "$_" 
	
	infile="${label}-therm.in"
	slfile="${label}-therm.slurm"

	# copy input file and slurm script to new dir
	cp $wd/input_file/$material/therm.in $infile
	cp $wd/slurm/$material/therm.slurm $slfile

	# copy pdb for reference
	cp $wd/filled_tube/$label.pdb .
	cp $wd/input_file/$material/CH.rebo .

	# edit files
	sed -i "s/variable label string/variable label string $label/1" $infile
	sed -i "s/<label>/$label/1" $slfile
	echo "cd $wd/sim/$label/therm" >> $slfile
	echo "${lmprun} < ${infile}" >> $slfile

	# submit to slurm
	sbatch $label-therm.slurm
	#sleep 2

	# reset to wd
	cd $wd

done
