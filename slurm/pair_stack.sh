#!/bin/bash

#SBATCH --job-name=stack                  # Job name
#SBATCH --nodelist=node02                 # Run all processes on a single node
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --cpus-per-task=72                # Number of CPU cores per task
#SBATCH --mem=20gb                        # Job memory request
#SBATCH --time=72:00:00                   # Time limit hrs:min:sec
#SBATCH --partition=batch                 # Partition name
#SBATCH --output=log_stack_pair_%j.log    # Standard output and error log

pwd; hostname; date

echo "Running prime number generator program on $SLURM_CPUS_ON_NODE CPU cores"

###################### Notes ########################
## EXAMPLE slurm script for running pair_stack.py. ##
#####################################################

module load openmpi/4.1.4 anaconda/3.9

# --- Script Configuration ---

# Stacking parameters
export NFS=35                   # Number of frequency slices to extract, equivalent to a frequency width:  2*NFS+1 * freq_resolution
export NWORKER=72               # Number of workers for multiprocessing (should be <= --cpus-per-task)
export SSIZE=500                # Split size for pair catalog processing, results in the same split will be stacked together.
export RANDOM_FLIP="True"       # Randomly flip individual pair map (True/False)
export HALFWIDTH="3.0"          # Stack result map half-width
export NPIX_X="120"             # Stack result map X pixels
export NPIX_Y="120"             # Stack result map Y pixels
# export SAVEKEYS="Signal,Mask"   # Datasets to save in the output file
# export COMPRESSION="gzip"       # Compression method for the output file (gzip/lz4)
# export SKIP_EXIST="False"       # Skip existing output files (True/False)

# Define base paths and prefixes
export INPUT_MAP_BASE="/home/dyliu/data/"
export INPUT_MAP_PREFIX="prepared_map_cube.h5"
export INPUT_MAP_KEYS="T,mask,f_bin_edge,x_bin_edge,y_bin_edge"
export INPUT_MAP_MASKED="True"    # True to read 'mask' dataset in the input map file. Change to false to directly use zore masking.

export INPUT_PAIRCAT_BASE="/home/dyliu/data/sdss_catalog/"
export INPUT_PAIRCAT_PREFIX="pair_catalog"
export INPUT_PAIRCAT_KEYS='is_ra,pos'

export OUTPUT_STACK_BASE="/home/dyliu/data/galaxy_pair_stack/"
export OUTPUT_STACK_PREFIX="stack_result_nfs"$NFS

# Run the Python script
hifigps-stack

date
