#!/bin/bash

#SBATCH --job-name=find_fuzzy             # Job name
#SBATCH --nodelist=node02                 # Run all processes on a single node
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --cpus-per-task=72                # Number of CPU cores per task
#SBATCH --mem=10gb                        # Job memory request
#SBATCH --time=10:00:00                   # Time limit hrs:min:sec
#SBATCH --partition=batch                 # Partition name
#SBATCH --output=log_find_fuzzy_%j.log    # Standard output and error log


pwd; hostname; date

echo "Running prime number generator program on $SLURM_CPUS_ON_NODE CPU cores"

###################### Notes ############################
## EXAMPLE slurm script for running find_fuzzy_par.py. ##
#########################################################

module load openmpi/3.1.6, anaconda/3.8

# Script Configuration
TNG_BASE_DIR="/home/DATA/ycli/tng/TNG-100/output/"
TNG_SNAP_NUM="91"
OUTPUT_PATH="/home/dyliu/Filament/Fuzzy_particles/fuzzy_par_indices.h5"

hifigps-find-fuzzy "$TNG_BASE_DIR" "$TNG_SNAP_NUM" -o "$OUTPUT_PATH"

date
