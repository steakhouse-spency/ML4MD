#!/bin/bash




# must set path of src

# HPC
# home="/home/rcf-proj/an2/sportega"
# wd="$home/ML4MD/spencer/"
# lmpsrc="$home/lammps/lammps/src"

# Local
wd="/home/sportega/Desktop/dr/spencer"
lmpsrc="/home/rcf-proj/an2/sportega/lammps/lammps/src"

# lammps mpi execution command
lmprun="srun --mpi=pmi2 $lmpsrc/lmp_foo"

# check arguments
# 	ex: ./runlmp.sh C 20 7 15
if [ $# -lt 5 ]; then
	echo "default arguments"
	material="C"
	L=20
	N=7
	M=15
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
datafile="data_file/$label.data"
# fffname="SiC_1989real.tersoff"

# variables for sim
SLEEP=4
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
# T=($(seq 270 20 290)) # test case (270 290)

# for every temperature
for t in "${T[@]}"; do

	# go to working directory
	cd $wd

	# path to simulation directory for current temperature
	t_path=sim/$label/$t

	# ?
	trun=$[(t-10)*500]

	# create directories (if dne)
	mkdir -p $t_path

	# for n = 1 to nruns
	runs=($(seq 1 $nruns))
	for n in "${runs[@]}"; do

		# get 4 random numbers and store in list
		declare -a rand
		for i in {0..3}
		do
			rand[$i]=$(get_randnum)
		done

		# copy data file to sim directory
		cp $datafile $t_path/$label-$t.data

		# for every part of the current sim (a..f)
		for i in {a,b,c,d,e,f}; do
			# create copy of input file and slurm script
			run_label=$label-${t}_${n}_$i
			echo $run_label
			infile=$t_path/$run_label.in
			slfile=$t_path/$run_label.slurm

			# copy template input file 
			cp input_file/$material/$i.in 		$infile
			cp slurm/$material/template.slurm 	$slfile

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
			sed -i '/#SBATCH --output=/c\#SBATCH --output='"$run_label"'.out' $slfile

			# append execution command to slurm script
			echo "cd ${t_path}" >> $slfile
			echo "${lmprun} < ${run_label}.in" >> $slfile 
		done

		# go to sim directory
		cd $t_path

		# submit jobs with dependencies:
		for i in {a,b,c,d,e}; do
			last_job=""
			# last_job=$(sbatch --parsable ${dependency} ${run_label}.slurm)
			if [ "$?" == "0" ]; then
				[[ $DEBUG == 1 ]] && echo "Info: JOBID: $last_job, $n submitted"
				dependency="--dependency=afterany:$last_job"
				# sleep $SLEEP
			else
				echo "Error exiting: sbatch $dependency ${job}${i}.slurm failed to run"
				#exit 1
			fi
		done
		dependency=""
	done
done
