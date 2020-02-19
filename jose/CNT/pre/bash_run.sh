#!/bin/bash
m="0"
SLEEP=4
DEBUG=1

#file names format:
#in-GK-Temp-run.in
#r-sicw-GK-Temp-run.slurm
#0_visco_from_pab_water_scf-run.py


#route to lammps executable and file names
lmprt="srun --mpi=pmi2 /auto/rcf-proj2/jgc/cobenare/lammps/lammps-22Aug18/src/lmp_foo"
dfilename="data-file-sicnt-water-0.65.data"
fffname="SiC_1989real.tersoff"
infile="in-GK"
runfile="r-sicw-GK"
#pyfile="0_visco_from_pab_water_scf"


#simulations variables:
FLOOR=1500
T="298"
nruns="30"


#start 

#create files for each value of velocity in vel
	
	v="1"
	#create replicates files as indicate in nruns
	while [ $v -le $nruns ]
	do
	
		mkdir "r${v}"
	
		# to get random numbers 1 (zeroing the array)
		declare -a arr
		for j in {0,1,2,3}; do
		arr[$j]=0
		done
		
		#to get random numbers 2 (filling up the array)
		for i in {0,1,2,3}
		do
		while [ "${arr[$i]}" -le $FLOOR ]
		do
		arr[$i]=$RANDOM
		done
		done
		
		#copy data file and ff file
		cp $dfilename r$v/$dfilename
		cp $fffname r$v/$fffname

		#copy files to a new directory (one per velocity value)
		cp $infile-$T-$m.in  r$v/$infile-$T-$v.in 
		cp $runfile-$T-$m.slurm  r$v/$runfile-$T-$v.slurm
		#cp $pyfile${m}.py r$v/$pyfile${v}.py
		
		cd r$v
		
		#changes in input file
		sed -i '/variable nrun equal .*/c\variable nrun equal '"$v"'' $infile-$T-$v.in 
		sed -i '/fix addtemp1sicnt 9SiCNT_no_anchors langevin 0.0 10.0 5.0 .*/c\fix addtemp1sicnt 9SiCNT_no_anchors langevin 0.0 10.0 5.0 '"${arr[0]}"'' $infile-$T-$v.in 
		sed -i '/fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 .*/c\fix addtemp1water 5Water_atoms	langevin 0.0 10.0 5.0 '"${arr[1]}"'' $infile-$T-$v.in 
		sed -i '/fix addtemp3sicnt 9SiCNT_no_anchors langevin 10.0 ${T} 5.0 .*/c\fix addtemp3sicnt 9SiCNT_no_anchors langevin 10.0 ${T} 5.0 '"${arr[2]}"'' $infile-$T-$v.in 
		sed -i '/fix addtemp3water 5Water_atoms	langevin 10.0 ${T} 5.0 .*/c\fix addtemp3water 5Water_atoms	langevin 10.0 ${T} 5.0 '"${arr[3]}"'' $infile-$T-$v.in 
		
		
		#changes in slurm file, including the python file submission
		sed -i 's|'"${lmprt}"' < '"${infile}"'-'"$T"'-'"$m"'.in|'"${lmprt}"' < '"${infile}"'-'"$T"'-'"$v"'.in|g' $runfile-$T-$v.slurm
		#sed -i 's|python3 0_visco_from_pab_water_scf'"$m"'.py|python3 0_visco_from_pab_water_scf'"$v"'.py|g' $runfile-$T-$v.slurm
		
		#changes in py file
		#sed -i 's|run = '"$m"'|run = '"$v"'|g' $pyfile$v.py


		#submit jobs
		last_job=$(sbatch --parsable "$runfile-$T-$v.slurm")
		if [ "$?" == "0" ]; then
		[[ $DEBUG == 1 ]] && echo "Info: JOBID: $last_job, $v submitted"
		
		sleep $SLEEP
		
		v=$[v+1]
		
		else
		echo "there is a problem"
		fi
		cd ..
	done
	
exit 1
