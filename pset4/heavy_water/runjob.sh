#!/bin/sh
 
#### Torque directives:
 
# Give the job a name.
#PBS -N 211_pset
 
# Use the default queue.
#PBS -q default
 
# Request 2 nodes and 12 processors per node.
### #PBS -l nodes=2:ppn=12
 
# Request walltime to run this job.
### #PBS -l walltime=00:30:00
 
 
#### Set the working directory and path.
 
# This cd command tells Torque to execute your code in the same directory that
# you submitted the job in.  So if the input files for your code are in
# /home/bob/openmc_sims/sim5, make sure you cd into that directory before
# running qsub.
cd $PBS_O_WORKDIR
 
# This command will tell Torque to use the same PATH that you were using when
# you ran qsub.  It's commented out because it is not necessary in this case,
# but you may need to uncomment it to get your code to run right.
#PATH=$PBS_O_PATH
 
 
#### Run the simulation.
 
module load openmc
python execution_script.py
