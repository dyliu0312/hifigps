"""
Find INNER fuzzy particles that bound to halos but not to any subhalos.
Those were regraded as the particles of pure filaments.

This script calculates the start and end indices for fuzzy particles 
for a given particle type in an IllustrisTNG simulation snapshot.

The dataset files are assumed to be organized by with TNG default structure, Group/Halo, and then 
by Subhalo within each Group.
"""
import sys
import time
import os
import illustris_python as il
from hifigps.data import save_h5 

def count_part_num(base:str, snap:int=91, part_type:str='gas'):
    """
    Counts the total number of particles for a given type in each halo 
    and the number of those particles contained within its subhalos.
    
    Args:
        base: The base path saving the TNG simulation data (e.g., contains 'groups_xx/').
        snap: The number of the snapshot. Default is 91.
        part_type: The particle type to count (e.g., 'gas', 'dm', 'star'). Default is 'gas'.

    Returns:
        A tuple: (list of total particle count per halo, list of subhalo particle count per halo).
    """
    # Fields needed from the Group Catalog for halos (groups)
    halo_fields = ['GroupLenType','GroupFirstSub','GroupNsubs']
    # Fields needed from the Group Catalog for subhalos
    subhalo_fields = ['SubhaloGrNr','SubhaloParent','SubhaloLenType']

    # Load Group and Subhalo data
    halos = il.groupcat.loadHalos(base,snap,fields=halo_fields)
    subhalos = il.groupcat.loadSubhalos(base,snap,fields=subhalo_fields)
    
    # Get the index corresponding to the particle type (e.g., 0 for 'gas')
    part_num = il.snapshot.partTypeNum(part_type) 

    subhalo_npart = [] # List to store total particle count *within* subhalos for each halo
    halo_npart = []    # List to store total particle count for each halo

    # Iterate through all halos
    for i, h_npart in enumerate(halos['GroupLenType'][:,part_num]):
        # 'GroupLenType' is the total number of particles of 'part_num' in the halo
        
        fsub = halos['GroupFirstSub'][i] # Index of the first subhalo in the group
        
        if (fsub != -1): # Check if the halo has any subhalos
            n_subs = halos['GroupNsubs'][i] # Total number of subhalos in the group
            
            # Sum the particle counts for the given type across all subhalos of this halo
            # SubhaloLenType is an array [N_subhalos, 6 particle types]
            subh_npart = sum(subhalos['SubhaloLenType'][fsub:fsub+n_subs, part_num])
        else:
            subh_npart = 0
            
        subhalo_npart.append(subh_npart)
        halo_npart.append(h_npart)

    return halo_npart, subhalo_npart

def get_fuzzy_indices(halo_npart:list, subhalo_npart:list):
    """
    Calculates the start and end indices of the fuzzy particles for each halo.

    The particle data is assumed to be ordered sequentially:
    ... [Subhalo 1 particles] [Subhalo 2 particles] ... [Fuzz particles] ...
    
    Args:
        halo_npart: A list containing the total number of particles for each halo.
        subhalo_npart: A list containing the total number of particles within subhalos for each halo.

    Returns:
        A tuple: (start indices list, end indices list) of the fuzz particles for each halo.
    """
    st_inds = [] # List for start indices of the fuzz component
    ed_inds = [] # List for end indices of the fuzz component

    temp_sum = 0 # Cumulative sum of particles processed so far (total particles in previous halos)

    # st = temp_sum + j (subhalo particles) -> This is where the fuzz component *starts*
    # ed = temp_sum + i (halo particles)    -> This is where the fuzz component *ends*
    # (Since particles are assumed to be subhalos first, then fuzz)
    for i,j in zip(halo_npart, subhalo_npart):
        # The starting index of the fuzz component for the current halo.
        # It's after all particles of all previous halos, plus all subhalo particles in the current halo.
        st = temp_sum + j 
        
        # The ending index of the fuzz component (which is the end of the current halo's particles).
        ed = temp_sum + i
        
        # Update the cumulative sum for the next halo
        temp_sum = ed
        
        st_inds.append(st)
        ed_inds.append(ed)
    return st_inds, ed_inds


if __name__ == "__main__":
    # --- Read input arguments from environment variables ---
    # Using os.getenv allows execution without command-line arguments if variables are set
    base_dir = os.getenv('TNG_BASE_DIR')
    snap_num_str = os.getenv('TNG_SNAP_NUM')
    output = os.getenv('OUTPUT_PATH')
    
    if not all([base_dir, snap_num_str, output]):
        print('Usage: Set environment variables and run the script:')
        print('export TNG_BASE_DIR="/path/to/TNG/simulation/"')
        print('export TNG_SNAP_NUM="91"')
        print('export OUTPUT_PATH="/path/to/output/fuzz_particles.hdf5"')
        print('python find_fuzz_particles.py')
        sys.exit()
    
    # --- Input Processing and Setup ---
    print(f'Loading data from BASE_DIR: {base_dir}')
    print(f'Using SNAP_NUM: {snap_num_str}')
    print(f'Saving fuzz indices result at OUTPUT_PATH: {output}')
    try:
        snap_num = int(snap_num_str) # type: ignore
    except ValueError:
        print(f"Error: SNAP_NUM must be an integer, got '{snap_num_str}'")
        sys.exit()

    # --- Default Arguments ---
    PTYPE = 'gas' 
    KEYS = ['start_index', 'end_index'] # Keys for the output HDF5 file

    # --- Count Particles Number ---
    print(f'--- Starting particle counting for PTYPE: {PTYPE} ---')
    t0 = time.time()
    halo_npart_list, subhalo_npart_list = count_part_num(base_dir, snap_num, PTYPE) # type: ignore
    print(f'Particle counting finished in {time.time()-t0:.2f} seconds')

    # --- Get Fuzz Indices ---
    print('--- Getting fuzz particle indices ---')
    t1 = time.time()
    st_ids, ed_ids = get_fuzzy_indices(halo_npart_list, subhalo_npart_list)
    print(f'Fuzz particle indices calculation finished in {time.time()-t1:.2f} seconds')

    # --- Save Indices --- 
    print(f'--- Saving fuzz indices to {output} ---')
    t2 = time.time()
    # save_h5 is a custom function, assuming it handles the HDF5 writing
    save_h5(output, KEYS, [st_ids, ed_ids]) # type: ignore
    print(f'Saving fuzz indices finished in {time.time()-t2:.2f} seconds')
    print(f'Total elapsed time: {time.time()-t0:.2f} seconds')