#!/bin/bash




# must set path of src

# HPC
wd="/staging/an2/sportega/ML4MD/spencer"
lmpexec="/home/rcf-40/sportega/disk/lammps/lammps/src/lmp_foo"

lmpsrc="/home/rcf-proj/an2/sportega/lammps/lammps/src"

# lammps mpi execution command
lmprun="srun --mpi=pmi2 $lmpexec"

# check arguments
# 	ex: ./runlmp.sh C 20 7 15
if [ $# -lt 5 ]; then
	echo "default arguments"
	material="C"
	L=29
	N=4
	M=11
	nruns=1
else
	# nanotube material
	material=$1
	# length of tube in nm
	L=$2
	# diameter n
	N=$3
	# diameter m
	M=$4
	# number of iterations
	nruns=$5
fi

# label for current nanotube
label="${material}_${L}_${N}_${M}"
#name of data file and force field file
# datafile="data_file/$label.data"
# fffname="SiC_1989real.tersoff"

# variables for sim
SLEEP=2
DEBUG=0
dependency=""


# function to generate random number > FLOOR
get_randnum () {
	local FLOOR=1500
	local num=0
	while [ $num -le $FLOOR ]
	do
		num=$RANDOM
	done
	echo $num
}



# create list of temperatures
# T=($(seq 270 20 430))
T=($(seq 270 20 290))
# T=($(seq 270 20 290)) # test case 
T=(270)



# for every temperature
for t in "${T[@]}"; do

	# path to simulation directory for current temperature
	t_path=$wd/sim/$label/$t

	# ?
	trun=$[(t-10)*500]

	# create Temp directory (if dne) and enter
	mkdir -p $t_path && cd "$_"

	# for n = 1 to nruns
	runs=($(seq 1 $nruns))
	for n in "${runs[@]}"; do

		# if nrun id has been used already
		if [ -d run_$n ]; then
			echo "$t_path/run_$n exists. Skipping."
			continue
		fi

		# enter nrun directory
		mkdir run_$n && cd "$_"

		cp $wd/sim/timesteps.times .

		# get 4 random numbers and store in list
		declare -a rand
		for i in {0..3}
		do
			rand[$i]=$(get_randnum)
		done

		# copy data file to sim directory
		# cp $datafile $t_path/$label-$t.data

		# for every part of the current sim (a..f)
		for i in {a,b,c,d,e,f}; do
			# create copy of input file and slurm script
			run_label=$label-${t}_${n}_$i
			echo $run_label

			# infile=$t_path/$run_label.in
			# slfile=$t_path/$run_label.slurm
			infile=$run_label.in
			slfile=$run_label.slurm

			# copy template input file 
			cp $wd/input_file/$material/$i.in 		$infile
			cp $wd/slurm/$material/template.slurm 	$slfile

			#changes in input file, only the one that has the langeving thermostat, that is, the a- file
			sed -i '/variable T equal .*/c\variable T equal '"$t"'' 			$infile
			sed -i '/variable nrun equal .*/c\variable nrun equal '"$n"'' 		$infile
			sed -i '/variable label string/c\variable label string '"$label"'' 	$infile
			sed -i '/run 298000 #change/c\run '"${trun}"' #change' 				$infile

			#changes only in file a:
			if [ $i == "a" ]; then
				sed -i '/fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 .*/c\fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 '"${rand[1]}"'' $infile
				sed -i '/fix addtemp2water 5Water_atoms	langevin 10.0 ${T} 5.0 .*/c\fix addtemp2water 5Water_atoms	langevin 10.0 ${T} 5.0 '"${rand[3]}"'' $infile
			fi

			# set output file
			# sed -i '/#SBATCH --output=/c\#SBATCH --output='"$run_label"'.out' $slfile
			sed -i 's|<i>|'"${i}"'|g' $slfile
			sed -i 's|<label>|'"${run_label}"'|g' $slfile

			# append execution command to slurm script
			echo "cd $t_path/run_$n" >> $slfile
			echo "${lmprun} < ${run_label}.in" >> $slfile 
		done

		# go to sim directory
		# cd $t_path

		# submit jobs with dependencies:
		for i in {a,b,c,d,e,f}; do
			run_label=$label-${t}_${n}_$i
			#last_job=""
			last_job=$(sbatch --parsable ${dependency} ${run_label}.slurm)
			echo $i: $last_job
			if [ "$?" == "0" ]; then
				echo "Info: JOBID: $last_job, $n submitted"
				dependency="--dependency=afterany:$last_job"
				sleep $SLEEP
			else
				echo "Error exiting: sbatch $dependency ${job}${i}.slurm failed to run"
				#exit 1
			fi
		done
		dependency=""
		cd ..
	done
done
