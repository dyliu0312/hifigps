"""
Find INNER fuzzy particles that bound to halos but not to any subhalos.
Those were regraded as the particles of pure filaments.

This script calculates the start and end indices for fuzzy particles
for a given particle type in an IllustrisTNG simulation snapshot.

The dataset files are assumed to be organized by with TNG default structure, Group/Halo, and then
by Subhalo within each Group.
"""
import argparse
import os
import time

import illustris_python as il

from hifigps.data import save_h5


def parse_args():
    parser = argparse.ArgumentParser(
        description="Find inner fuzzy particles in IllustrisTNG simulations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Example:
  hifigps-find-fuzzy /path/to/TNG/simulation/ 91 -o /output/fuzz_particles.hdf5
""",
    )
    parser.add_argument("base", help="Base path to TNG simulation data")
    parser.add_argument("snap", type=int, help="Snapshot number")
    parser.add_argument(
        "-o", "--output", required=True, help="Output file path for fuzz particle indices"
    )
    parser.add_argument(
        "-p",
        "--part-type",
        default="gas",
        choices=["gas", "dm", "star", "bhs", "bhw", "stars"],
        help="Particle type (default: gas)",
    )
    return parser.parse_args() 

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


def main():
    args = parse_args()

    base_dir = args.base
    snap_num = args.snap
    output = args.output
    ptype = args.part_type

    print(f"Loading data from BASE_DIR: {base_dir}")
    print(f"Using SNAP_NUM: {snap_num}")
    print(f"Saving fuzz indices result at OUTPUT_PATH: {output}")
    print(f"Particle type: {ptype}")

    PTYPE = ptype
    KEYS = ["start_index", "end_index"]

    print(f"--- Starting particle counting for PTYPE: {PTYPE} ---")
    t0 = time.time()
    halo_npart_list, subhalo_npart_list = count_part_num(base_dir, snap_num, PTYPE)
    print(f"Particle counting finished in {time.time() - t0:.2f} seconds")

    print("--- Getting fuzz particle indices ---")
    t1 = time.time()
    st_ids, ed_ids = get_fuzzy_indices(halo_npart_list, subhalo_npart_list)
    print(f"Fuzz particle indices calculation finished in {time.time() - t1:.2f} seconds")

    print(f"--- Saving fuzz indices to {output} ---")
    t2 = time.time()
    save_h5(output, KEYS, [st_ids, ed_ids])
    print(f"Saving fuzz indices finished in {time.time() - t2:.2f} seconds")
    print(f"Total elapsed time: {time.time() - t0:.2f} seconds")


if __name__ == "__main__":
    main()