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
 datafiles=("C_25_7_14-a")

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


	# get count of atoms for both nanotube and water
	# format: (count_1, id_1, ... , count_n, id_n)
	atom_count=($(cat $label.pdb | awk '{if($1=="ATOM"){print $3}}' | sort | uniq -c ))

	# for every other index
	nums=($(seq 1 2 ${#atom_count[*]}))
	for i in ${nums[*]}; do
		# get values
		count=${atom_count[$((i-1))]}
		id=${atom_count[$((i))]}
		
		# store counts
		if [ $id == $material ]; then
			tube_count=$count
		else
			water_count=$count
		fi
	done


	echo "tube:" $tube_count
	echo "water:" $water_count

	# percentages represent where rings will be placed
	# no ring at 1.0
	rings=(0.25 0.5 0.75 1.0)

	# store all non anchor ids
	non_anchor=( )

	# for every ring
	start=1
	count=1
	for ring in ${rings[*]}; do
		# get id of atom in middle of ring
		# bash cant do float arithmitic -__-
		i=$(echo "scale=2; $tube_count*$ring" | bc)
		# floor i value to create index
		i=${i%.*}

		# get range of atoms in ring
		# using 10 atoms/ring 
		range[0]=$((i-4))
		range[1]=$((i+5))

		# debug
		echo $ring ":" $i
		echo "range:" ${range[*]}
		
		# append non-anchor atom ids
			# from end of last ring to end of tube
		if [ $count -eq ${#rings[*]} ]; then
			echo "last iter!"
			end=$tube_count
			# from 'start' to the first atom in this ring
		else
			end=${range[0]}
			count=$((count + 1))
		fi
		non_anchor+=($(seq $start $end))
		
		# update start to the last atom in this ring + 1
		start=$((range[1] + 1))
	done

	# get the id of the first water atom (H/O)
	water_start=$((tube_count+1))
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
	sbatch $label-therm.slurm
	#sleep 2

	# reset to wd
	cd $wd

done
