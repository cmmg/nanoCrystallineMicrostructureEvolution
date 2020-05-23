#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=compphys # default "univ", if not specified
#SBATCH --time=0-12:00:00 # run time in days-hh:mm:ss
#SBATCH --nodes=2# require 2 nodes
#SBATCH --ntasks-per-node=20            # (by default, "ntasks"="cpus")
#SBATCH --mem-per-cpu=4000# RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --job-name="DBPG090"
#SBATCH --error=file%j.err
#SBATCH --output=file%j.out
#SBATCH --mail-user=ppandey7@wisc.edu
#SBATCH --mail-type=ALL
#Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

#Now list your executable command (or a string of them).
# Example for non-SLURM-compiled code:
source ~/.bashrc
cd /home/ppandey7/workspace/repos/AROproject/prlGns0-90-tonk
rm -rf CMakeCache.txt cmake_install.cmake Makefile CMakeFiles 
cmake .
make release
mpirun -np 40 ./main
