#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC164847 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --error=EWAS.err # error file
#SBATCH --output=EWAS.log # output file
#SBATCH --job-name=EWAS


#------------------------------------------------------

# 1. command line argument input is:

# to run sbatch EWASjobSubmission.sh <filepath/to/projectFolder>

#-----------------------------------------------------

## print start date and time
echo Job started on:
  date -u
JOBNAME="EWAS"

echo Job sumbitted from:
  echo $SLURM_SUBMIT_DIR

# Move the user to the project directory
cd $1

# load config file for $RVERS
source config.txt

## load modules
eval "$(conda shell.bash hook)"
conda activate EWAS

# run scripts - CHANGE DEPENDING ON ANALYSIS
Rscript APOE_4_NO_BS_script.r 

## print finish date and time
echo Job finished on:
  date -u
