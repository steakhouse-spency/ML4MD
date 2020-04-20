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
datafiles=("C_25_7_14")

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

	# copy pdb for reference
	cp $wd/filled_tube/$label.pdb .
	cp $wd/input_file/$material/CH.rebo .


	# # get count of atoms for both nanotube and water
	# # format: (count_1, id_1, ... , count_n, id_n)
	# atom_count=($(cat $label.pdb | awk '{if($1=="ATOM"){print $3}}' | sort | uniq -c ))

	# # for every other index
	# nums=($(seq 1 2 ${#atom_count[*]}))
	# for i in ${nums[*]}; do
	# 	# get values
	# 	count=${atom_count[$((i-1))]}
	# 	id=${atom_count[$((i))]}
		
	# 	# store counts
	# 	if [ $id == $material ]; then
	# 		tube_count=$count
	# 	else
	# 		water_count=$count
	# 	fi
	# done


	# echo "tube:" $tube_count
	# echo "water:" $water_count

	# # percentages represent where rings will be placed
	# # no ring at 1.0
	# rings=(0.25 0.5 0.75 1.0)

	# # atoms per ring
	# ring_len=10

	# # store all non anchor ids
	# non_anchor=( )

	# # for every ring
	# start=$((ring_len + 1))
	# half_ring=$(echo "$ring_len * 0.5 / 1" | bc)
	# for ring in ${rings[*]}; do
	# 	# get id of atom in middle of ring
	# 	# bash cant do float arithmitic -__-
	# 	# divide by 1 to floor decimal value
	# 	mid=$(echo "$tube_count * $ring / 1" | bc)
		
	# 	# get range of atoms in ring
	# 	# using 10 atoms/ring 
	# 	range[0]=$(echo "$mid - $half_ring - 1" | bc)
	# 	range[1]=$((mid + half_ring))

	# 	# debug
	# 	echo ring middle: $mid
	# 	echo ring range: ${range[*]}
	
	# 	# append non-anchor atom ids
	# 	if [ "$ring" == "1.0" ]; then
	# 		# from end of last ring to end of tube
	# 		# -10 for last ring at end of tube
	# 		#echo "last iter"
	# 		end=$((tube_count - 10))
	# 	else
	# 		# from 'start' to the first atom in this ring
	# 		end=$((range[0] - 1))
	# 	fi
	# 	non_anchor+=($(seq $start $end))
	# 	echo non-anchors: [$start, $end]		
	# 	echo --------
	# 	# update start to the last atom in this ring + 1
	# 	start=$((range[1] + 1))
	# done


	output=($(python3 $wd/sim/get_nonanchors.py $label))
	water_count=${output[0]}
	tube_count=${output[1]}
	non_anchor=${output[2]}


	# get the id of the first water atom (H/O)
	water_start=$((tube_count + 1))
	# increment by number of water atoms
	water_end=$((water_start + water_count))
	# create sequence of water ids
	water=($(seq $water_start $water_end))


	# output to file to check
	# echo ${non_anchor[*]} > non_anchor.txt
	# echo ${water[*]} > water.txt


	# edit input file
	sed -i "s|<non_anchor>|${non_anchor[*]}|g" $label-therm.in
	sed -i "s|<water>|${water[*]}|g" $label-therm.in

	sed -i "s/variable label string/variable label string $label/1" $infile
	sed -i "s/<label>/$label/1" $slfile
	echo "cd $wd/sim/$label/therm" >> $slfile
	echo "${lmprun} < ${infile}" >> $slfile

	# submit to slurm
	# sbatch $label-therm.slurm
	#sleep 2

	# reset to wd
	cd $wd

done
