# Setting up confined Water in Nanotube Simulation

## Fill Nanotube with Water
Before we can perform our simulation we must first add water atoms to our empty NTs.

1. Create an empty NT PDB file of your choice via some NT builder.

2. Change PDB name to match our naming scheme.  
Ex. &nbsp;&nbsp;M_L_n_m.pdb  
&nbsp;&nbsp;&nbsp;&nbsp;M:&nbsp;&nbsp;Material Type (C, SiC, etc.)  
&nbsp;&nbsp;&nbsp;&nbsp;L: &nbsp;&nbsp;NT length (nm)  
&nbsp;&nbsp;&nbsp;&nbsp;n: &nbsp;&nbsp;NT 'n' value (diameter)  
&nbsp;&nbsp;&nbsp;&nbsp;m: &nbsp;&nbsp;NT 'm' value (diameter)   


3. Place this PDB inside the directory 'empty_tube/'.

4. Enter 'src' directory. Open 'globalvar.py' and update variables to match your current NT. 

5. Once all parameter are set for the current NT you can now execute 'python3 addwater.py' 

6. Your new PDB file and LAMMPS data file can be found in  directories 'filled_tube/' and 'data_file/'

## Thermalize Nanotube
After we have succefully filled our NTs, we can now thermalize them. One thing to note is that you only need to do this process once for every unique NT. These results are then used to run multiple  simulations with different temperatures for this specific NT.

1. Enter 'sim/' directory. Open 'runtherm.sh' and update variables to match NT and local machine.

2. Execute './runtherm.sh'

3. A working directory will be created in 'sim/' for this NT and a thermalize simulation will be submitted as a batch job to the HPC Slurm scheduler.

## Simulate confined water
Once your desired NT has been thermalized, we can now begin running our long simulations.

1. Enter 'sim' directory. Open 'runlmp.sh' and update array 'T' to contain all desired temperatures to simulate.

2. Execute './runlmp'

3. A new directory will be created within 'sim/<NT>/' and all confined water simulations (a-f) will be submitted as a batch job to the HPC Slurm scheduler for all desired temperatures.
