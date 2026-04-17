#!/usr/bin/env python3

from get_contacts import * 
from vips import * 

import os
import glob
import argparse
import multiprocessing

import pandas as pd
import numpy as np
import math

# import datetime
# import json
# from pymol import cmd, stored 
# from pyrosetta import * 
# from pyrosetta.rosetta.core.scoring import *
# from pyrosetta.rosetta.protocols.relax import FastRelax

# init("-mute all")
# scorefxn = create_score_function("ref2015_cart")
# relax = FastRelax() 
# relax.set_scorefxn(scorefxn)

def mean_distance_to_vips(pdbfilename):

    design_name = os.path.splitext(os.path.basename(pdbfilename))[0] # match all 'design_id' from plddt json
    contactsdf, _ = pymol_find_contacts(pdbfilename) # , 100)
    print(f"successfully found contacts for {pdbfilename}")
    contactsdf['chainA_resni'] = contactsdf['chainA_resn'] + contactsdf['chainA_resi'].astype(str)
    
    # don't want the average of every atom, only atoms of receptor residues
    filtered_contactsdf = contactsdf[contactsdf['chainA_resni'].isin(receptor_vip_resis)]
    mean_dist = np.mean(filtered_contactsdf['distance'])

    return design_name, pymol_get_chainB_sequence(pdbfilename), mean_dist

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input-dir", "-i", type=str, help="Path to dir containting input pdbs")
    args = parser.parse_args()

    args.input_dir = args.input_dir.rstrip('/') 

    with open(f"{args.input_dir}/gendist_searchparams_{date_time_string}.txt", 'w') as fout: 
        fout.write(f"receptor residues: {receptor_vip_resis}\n")
        fout.write(f"receptor atoms: {receptor_vip_atoms}\n")
        fout.write(f"receptor hydrophobic residues: {receptor_hydrophobic_resis}\n")
        fout.write(f"receptor hydrophobic atoms: {receptor_hydrophobic_atoms}\n")
        fout.write(f"hydrophibic pairs: {vip_hydrophobic}\n") 

    all_inputs = glob.glob(f'{args.input_dir}/*.pdb')
    # get info abt how many cores avail: 
    npools = min(multiprocessing.cpu_count() - 1, len(all_inputs))

    # calc what chunk size should be for large input list
    chunksize = math.ceil(len(all_inputs) / npools) if npools < len(all_inputs) else 1

    with multiprocessing.Pool(processes=npools) as pool:
        results = pool.map(mean_distance_to_vips, all_inputs, chunksize=chunksize)

    ratingsdf = pd.DataFrame(columns=['design_id',
                                      'peptide_sequence', 
                                      'average_distance', ])

    for result in results: 
        ratingsdf = pd.concat([ratingsdf if not ratingsdf.empty else None, 
                               pd.DataFrame([{'design_id': result[0], 
                                              'peptide_sequence': result[1], 
                                              'average_distance': result[2],}])
                               ], 
                              ignore_index=True)
    ratingsdf.to_csv(f"{args.input_dir}/gendistances_{date_time_string}.csv")
