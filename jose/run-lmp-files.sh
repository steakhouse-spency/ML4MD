#!/bin/bash
m="0"
u="000"
#n="1"
SLEEP=4
DEBUG=1
dependency=""

lmpsrc="~/disk/lammps/"

#route to lammps executable
lmprt="srun --mpi=pmi2 /auto/rcf-proj2/jgc/cobenare/lammps/lammps-12Dec18/src/lmp_foo"

#name of data file and force field file
dfilename="data-file-sicnt-water-0.6.data"
# fffname="SiC_1989real.tersoff"

#simulations variables:
infile="in-sicnt-wat-prod"   #"-nm-10_3-0"
runfile="r-sicwat"           #"-298-nm-10_3-0-a"
timefile="time-steps-isf.times"
FLOOR=1500


nruns="1"
nametube="nm-10_3"

declare -a T=("220" "230")
declare -a T=("270" "290" "310" ... "370")
#declare -a T=("50" "60" "70" "80" "90" "100" "110" "120" "130" "140" "150" "160" "170" "180" "190" "200" "210" "220" "230" "240" "250" "260" "270" "280" "290" "300" "310" "320" "330" "340" "350" "360" "370" "380" "390" "400" "410" "420" "430" "440" "450" "460")


#here it starts 

#create files for each value of velocity in vel
for v in "${T[@]}"; do
	mkdir "${v}"
	
	n="1"
	#create replicates files as indicate in nruns
	while [ $n -le $nruns ]
	do
	
		# to get random numbers step 1
		declare -a arr
		for j in {0,1,2,3}; do
		arr[$j]=0
		done
		
		#to get random numbers step 2
		for i in {0,1,2,3}
		do
		while [ "${arr[$i]}" -le $FLOOR ]
		do
		arr[$i]=$RANDOM
		done
		done

		#---temp problem:
		trun=$[(v-10)*500]
		
		#copy data file and ff file and timefiles
		cp $dfilename $v/$dfilename
		cp $fffname $v/$fffname
		cp $timefile $v/$timefile

		#copy input file and the slurm file for running 
		#cp $infile-$u-${nametube}-$m.in  $v/$infile-$v-${nametube}-$n.in
		#cp $runfile-$u-${nametube}-$m.slurm  $v/$runfile-$v-${nametube}-$n.slurm


		#copy the input files and the continuation input files:
		for i in {a,b,c,d,e}; do
			cp $infile-$u-${nametube}-$m-${i}.in  $v/$infile-$v-${nametube}-$n-$i.in
			cp $runfile-$u-${nametube}-$m-${i}.slurm  $v/$runfile-$v-${nametube}-$n-$i.slurm
		done
		cd $v
		
		#changes in input file, only the one that has the langeving thermostat, that is, the a- file
		for i in {a,b,c,d,e}; do
		sed -i '/variable T equal .*/c\variable T equal '"$v"'' ${infile}-$v-${nametube}-$n-${i}.in
		sed -i '/variable nrun equal .*/c\variable nrun equal '"$n"'' ${infile}-$v-${nametube}-$n-${i}.in
		sed -i '/run 298000 #change/c\run '"${trun}"' #change' ${infile}-$v-${nametube}-$n-${i}.in
		done

        #changes only in file a:
		sed -i '/fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 .*/c\fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 '"${arr[1]}"'' ${infile}-$v-${nametube}-$n-a.in
		sed -i '/fix addtemp2water 5Water_atoms	langevin 10.0 ${T} 5.0 .*/c\fix addtemp2water 5Water_atoms	langevin 10.0 ${T} 5.0 '"${arr[3]}"'' ${infile}-$v-${nametube}-$n-a.in
		#done
		
		#changes in slurm files
		for i in {a,b,c,d,e}; do
        sed -i 's|'"${lmprt}"' < '"${infile}"'-'"${u}"'-'"${nametube}"'-'"$m"'-'"${i}"'.in|'"${lmprt}"' < '"${infile}"'-'"${v}"'-'"${nametube}"'-'"$n"'-'"${i}"'.in|g' $runfile-$v-${nametube}-$n-${i}.slurm
		done

		#sed -i 's|'"${lmprt}"' < '"${infile}"'-'"${u}"'-'"${nametube}"'-'"$m"'-a.in|'"${lmprt}"' < '"${infile}"'-'"${v}"'-'"${nametube}"'-'"$n"'-a.in|g' $runfile-$v-${nametube}-$n-a.slurm
		#sed -i 's|'"${lmprt}"' < '"${infile}"'-'"${u}"'-'"${nametube}"'-'"$m"'-b.in|'"${lmprt}"' < '"${infile}"'-'"${v}"'-'"${nametube}"'-'"$n"'-b.in|g' $runfile-$v-${nametube}-$n-b.slurm
		#sed -i 's|'"${lmprt}"' < '"${infile}"'-'"${u}"'-'"${nametube}"'-'"$m"'-c.in|'"${lmprt}"' < '"${infile}"'-'"${v}"'-'"${nametube}"'-'"$n"'-c.in|g' $runfile-$v-${nametube}-$n-c.slurm
		#sed -i 's|'"${lmprt}"' < '"${infile}"'-'"${u}"'-'"${nametube}"'-'"$m"'-d.in|'"${lmprt}"' < '"${infile}"'-'"${v}"'-'"${nametube}"'-'"$n"'-d.in|g' $runfile-$v-${nametube}-$n-d.slurm
		#sed -i 's|'"${lmprt}"' < '"${infile}"'-'"${u}"'-'"${nametube}"'-'"$m"'-e.in|'"${lmprt}"' < '"${infile}"'-'"${v}"'-'"${nametube}"'-'"$n"'-e.in|g' $runfile-$v-${nametube}-$n-e.slurm


		#last_job="0"
		#submitting the jobs with dependencies:
		for i in {a,b,c,d,e}; do
		#for i in {a,b}; do		
			last_job=$(sbatch --parsable ${dependency} "$runfile-$v-${nametube}-$n-${i}.slurm")
			if [ "$?" == "0" ]; then
				[[ $DEBUG == 1 ]] && echo "Info: JOBID: $last_job, $n submitted"
				dependency="--dependency=afterany:$last_job"
				sleep $SLEEP
			else
				echo "Error exiting: sbatch $dependency ${job}${i}.slurm failed to run"
				#exit 1
			fi
		done

		dependency=""
		#last_job="0"
		#last_job=$(sbatch --parsable "$runfile-$T-$v-$n.slurm")
		#if [ "$?" == "0" ]; then
		#[[ $DEBUG == 1 ]] && echo "Info: JOBID: $last_job, $n submitted"
		
		#sleep $SLEEP
		
		n=$[n+1]
		
		#else
		#echo "there is a problem"
		#fi
		cd ..
	done
	
done
exit 1
