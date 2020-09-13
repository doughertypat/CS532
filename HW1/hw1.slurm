#!/bin/bash -l
#-----------------------------------
# hw1.slurm
#
# Submit: 
#   sbatch hw1.slurm
#
# Check status:
#   squeue | grep [userid]
#   squeue -u [userid]
#   squeue -j [jobid]
#-----------------------------------

#SBATCH --job-name=hw1
#SBATCH --partition=hodor-cpu
#SBATCH --export=ALL
#SBATCH --exclusive

# Sets the maximum time the job can run (hh:mm:ss).
#SBATCH --time=00:01:00

# Specifies nodes for the job. 
#SBATCH --nodes 1                # Number of nodes
#SBATCH --ntasks 2               # Number of MPI tasks

# Sets the output files.
#SBATCH --output=./HW1.hodor-cpu.%j.txt

# Job events (mail-type): begin, end, fail, all. 
##SBATCH --mail-type=fail
##SBATCH --mail-user=[firstname].[lastname]@und.edu

module list
echo ""

mpirun ./hw1_DoughertyP 5

