#!/bin/bash

# must set path of src

# HPC
 wd="/staging/an2/sportega/ML4MD/spencer"
 lmpexec="/home/rcf-40/sportega/disk/lammps/lammps/src/lmp_foo"

# Local
#wd="/home/sportega/Desktop/dr/spencer"
#lmpexec="/home/rcf-proj/an2/sportega/lammps/lammps/src/lmp_foo"

# lammps mpi execution command
lmprun="srun --mpi=pmi2 $lmpexec"

material="C"

cd $wd
#datafiles=($(ls data_file))
# datafiles=("C_25_7_14")
datafiles=("C_25_7_14-rebo" "C_25_7_14-norebo")

for file in ${datafiles[*]}; do
	echo $file

	# remove file extension for label
	label="${file%.*}"

	# get nanotube material
	material=$(echo $label |  cut -d'_' -f1)

	# mkdir and cd into new dir
	mkdir -p sim/$label/therm && cd "$_" 
	
	infile="${label}_therm.in"
	slfile="${label}_therm.slurm"

	opt=$(echo $file | cut -d'-' -f2)
	echo $opt

	# copy input file and slurm script to new dir
	cp $wd/input_file/$material/"${opt}_therm.in" $infile
	cp $wd/slurm/$material/therm.slurm $slfile
	
	# copy other files
	cp $wd/filled_tube/$label.pdb .
	cp $wd/data_file/$label.data .
	cp $wd/input_file/$material/CH.rebo .

	output=($(python3 $wd/sim/get_nonanchors.py $label))
	tube_count=${output[0]}
	water_count=${output[1]}
	non_anchor=("${output[*]:2}")

	#echo $tube_count
	#echo $water_count
	#echo $non_anchor
	#exit 1

	# get the id of the first water atom (H/O)
	water_start=$((tube_count + 1))
	# increment by number of water atoms
	water_end=$((water_start + water_count))
	# create sequence of water ids
	water=($(seq $water_start $water_end))
	
	#echo $water
	#exit 1

	# edit input file
	sed -i "s|<non_anchor>|$non_anchor|g" $infile
	sed -i "s/<water>/${water[*]}/g" $infile

	sed -i "s/variable label string/variable label string $label/1" $infile
	sed -i "s/<label>/$label/1" $slfile
	echo "cd $wd/sim/$label/therm" >> $slfile
	echo "${lmprun} < ${infile}" >> $slfile

	# submit to slurm
    sbatch $slfile
	#sleep 2

	# reset to wd
	cd $wd

done
