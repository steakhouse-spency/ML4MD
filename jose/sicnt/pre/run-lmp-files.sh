#!/bin/bash
m="0"
n="1"
SLEEP=4
DEBUG=1

#route to lammps executable
lmprt="srun --mpi=pmi2 /auto/rcf-proj2/jgc/cobenare/lammps/lammps-12Dec18/src/lmp_foo"

#name of data file and force field file
dfilename="data-file-sicnt-water-0.65.data"
fffname="SiC_1989real.tersoff"

#simulations variables:
infile="in-vz"
runfile="r-sicwat"
tstepsf="time-steps-isf.times"
FLOOR=1500
T="230"
nruns="20"
v0="0.0000"

declare -a vel=("0.06" "0.05" "0.04")
#declare -a vel=("0.03" "0.00" "0.01" "0.02" "0.005")
#declare -a vel=("0.00005" "0.0001" "0.0005" "0.001" "0.0015" "0.002" "0.0025" "0.003")


#here it starts 

#create files for each value of velocity in vel
for v in "${vel[@]}"; do
	mkdir "${v}"
	
	n="1"
	#create replicates files as indicate in nruns
	while [ $n -le $nruns ]
	do
	
		# to get random numbers 1
		declare -a arr
		for j in {0,1,2,3}; do
		arr[$j]=0
		done
		
		#to get random numbers 2
		for i in {0,1,2,3}
		do
		while [ "${arr[$i]}" -le $FLOOR ]
		do
		arr[$i]=$RANDOM
		done
		done
		
		#copy data file and ff file
		
		cp $dfilename $v/$dfilename
		cp $fffname $v/$fffname
		cp ${tstepsf} $v/${tstepsf}

		#copy files to a new directory (one per velocity value)
		cp $infile-$T-${v0}-$m.in  $v/$infile-$T-$v-$n.in 
		cp $runfile-$T-${v0}-$m.slurm  $v/$runfile-$T-$v-$n.slurm
		
		cd $v
		
		#changes in input file		
		sed -i '/variable fvaluez equal .*/c\variable fvaluez equal '"$v"'' $infile-$T-$v-$n.in 
		sed -i '/variable nrun equal .*/c\variable nrun equal '"$n"'' $infile-$T-$v-$n.in 
		sed -i '/fix addtemp1sicnt 9SiCNT_no_anchors langevin 0.0 10.0 5.0 .*/c\fix addtemp1sicnt 9SiCNT_no_anchors langevin 0.0 10.0 5.0 '"${arr[0]}"'' $infile-$T-$v-$n.in 
		sed -i '/fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 .*/c\fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 '"${arr[1]}"'' $infile-$T-$v-$n.in 
		sed -i '/fix addtemp2sicnt 9SiCNT_no_anchors langevin 10.0 ${T} 5.0 .*/c\fix addtemp2sicnt 9SiCNT_no_anchors langevin 10.0 ${T} 5.0 '"${arr[2]}"'' $infile-$T-$v-$n.in 
		sed -i '/fix addtemp2water 5Water_atoms	langevin 10.0 ${T} 5.0 .*/c\fix addtemp2water 5Water_atoms	langevin 10.0 ${T} 5.0 '"${arr[3]}"'' $infile-$T-$v-$n.in 
		
		
		#changes in slurm file
		sed -i 's|'"${lmprt}"' < '"${infile}"'-'"$T"'-'"${v0}"'-'"$m"'.in|'"${lmprt}"' < '"${infile}"'-'"$T"'-'"$v"'-'"$n"'.in|g' $runfile-$T-$v-$n.slurm

		last_job=$(sbatch --parsable "$runfile-$T-$v-$n.slurm")
		if [ "$?" == "0" ]; then
		[[ $DEBUG == 1 ]] && echo "Info: JOBID: $last_job, $n submitted"
		
		sleep $SLEEP
		
		n=$[n+1]
		
		else
		echo "there is a problem"
		fi
		cd ..
	done
	
done
exit 1
