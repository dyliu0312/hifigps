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
# NFS: Number of frequency slices to extract, equivalent to a frequency width:  2*NFS+1 * freq_resolution
# NWORKER: Number of workers for multiprocessing (should be <= --cpus-per-task)
# SSIZE: Split size for pair catalog processing, results in the same split will be stacked together.
# RANDOM_FLIP: Randomly flip individual pair map (True/False)
# HALFWIDTH: Stack result map half-width
# NPIX_X/NPIX_Y: Stack result map pixels
# SAVEKEYS: Datasets to save in the output file, default to Signal and Mask
# COMPRESSION: Compression method for the output file (gzip/lz4)
# SKIP_EXIST: Skip existing output files (True/False)

NFS=35
NWORKER=72
SSIZE=500
RANDOM_FLIP="True"
HALFWIDTH="3.0"
NPIX_X="120"
NPIX_Y="120"
# SAVEKEYS="Signal,Mask" 
# COMPRESSION="gzip"
# SKIP_EXIST="False"

# Define base paths and prefixes
INPUT_MAP_BASE="/home/dyliu/data/"
INPUT_MAP_PREFIX="prepared_map_cube.h5"
INPUT_MAP_KEYS="T,mask,f_bin_edge,x_bin_edge,y_bin_edge"
INPUT_MAP_MASKED="True"    # True to read 'mask' dataset in the input map file. Change to false to directly use zore masking.

INPUT_PAIRCAT_BASE="/home/dyliu/data/sdss_catalog/"
INPUT_PAIRCAT_PREFIX="pair_catalog"
INPUT_PAIRCAT_KEYS='is_ra,pos'

OUTPUT_STACK_BASE="/home/dyliu/data/galaxy_pair_stack/"
OUTPUT_STACK_PREFIX="stack_result_nfs"$NFS

# Run the Python script
hifigps-stack \
    --map-base "$INPUT_MAP_BASE" \
    --map-prefix "$INPUT_MAP_PREFIX" \
    --paircat-base "$INPUT_PAIRCAT_BASE" \
    --paircat-prefix "$INPUT_PAIRCAT_PREFIX" \
    --out-base "$OUTPUT_STACK_BASE" \
    --out-prefix "$OUTPUT_STACK_PREFIX" \
    --nfs "$NFS" \
    --ssize "$SSIZE" \
    --nworker "$NWORKER" \
    --random-flip "$RANDOM_FLIP" \
    --halfwidth "$HALFWIDTH" \
    --npix-x "$NPIX_X" \
    --npix-y "$NPIX_Y"

date
