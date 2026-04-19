#!/usr/bin/env python3

import json
import numpy as np
import glob
import math
import os
import argparse
from vips import * 
from get_contacts import * 

import multiprocessing 

from pyrosetta import init, pose_from_pdb, create_score_function, FastRelax

init()
scorefxn = create_score_function("ref2015_cart")
relax = FastRelax() 
relax.set_scorefxn(scorefxn)

def get_plddts_json(json_path: str) -> dict: 
    plddt_lookup = dict()
    with open(json_path, 'r') as plddt_json: 
        for line in plddt_json: 
            data = json.loads(line.strip())
            plddt_lookup.setdefault(data['design_id'], 
                                    (data['average_plddt'], data['residue_plddt']))
    return plddt_lookup

def score_structure(pdb_filename): 
    """ 
    For pdb_file, 
    return the result of rosetta scoring without relaxing.
    Return score: float
    """
    unrelaxed_pose = pose_from_pdb(pdb_filename)
    return scorefxn(unrelaxed_pose) 

def relax_and_score_structure(pdb_filename, relaxed_dir): 
    """ 
    For pdb_file and directory to save relaxed structures to, 
    return the relaxed pose name and the scores before and after relaxing.
    Return (pose_name: str, score_after: float, score_after: float)
    """

    unrelaxed_pose = pose_from_pdb(pdb_filename)
    score_before_relax = scorefxn(unrelaxed_pose) 

    relaxed_file_name = f"{relaxed_dir}/{pdb_filename.rstrip('.pdb')}_relaxed.pdb"

    if not os.path.exists(relaxed_file_name): 
        relaxed_pose = unrelaxed_pose.clone()
        relax.apply(relaxed_pose) 
        score_after_relax = scorefxn(relaxed_pose) 
        relaxed_pose.dump_pdb(relaxed_file_name)
    else: 
        # Don't relax, just open previously saved file
        relaxed_pose = pose_from_pdb(relaxed_file_name)
        score_after_relax = scorefxn(relaxed_pose) 

    return relaxed_file_name, score_after_relax, score_before_relax

def count_contacts(pdbfilepath:str, cutoff_distance: float = 3.5):
    contactsdf, _ = pymol_find_contacts(pdbfilepath) 

    sequence = pymol_get_chainB_sequence(pdbfilepath)

    print(f"successfully found contacts for {pdbfilepath}") 
    
    # Keep contacts within cutoff distance
    filtereddf = contactsdf[contactsdf['distance'] <= cutoff_distance]

    # create a set of the chain A contact resis and atoms from contactsdf
    chA_contactresis = set(filtereddf['chainA_resn'] + filtereddf['chainA_resi']) # 3 letter code + number in seq
    chA_contactatoms = set(filtereddf['chainA_atomID'].astype(int)) 
    # for hydrophobic pair identification 
    chAchB_contactpairs = set(filtereddf['chainA_resn'] + filtereddf['chainA_resi'] + filtereddf['chainB_resn']) # 3 letter code + number in seq + _peptideres

    total_res_contacts_seen = len(chA_contactresis)
    num_vip_res_contacts_seen = len(chA_contactresis.intersection(receptor_vip_resis))
    num_atom_contacts_seen = len(chA_contactatoms.intersection(receptor_vip_atoms))
    num_hydrophobic_seen = len(chAchB_contactpairs.intersection(vip_hydrophobic))

    result = (sequence,
              total_res_contacts_seen, 
              num_vip_res_contacts_seen, 
              num_atom_contacts_seen,
              num_hydrophobic_seen)
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input-dir", "-i", type=str, help="Path to dir containting input pdbs")
    parser.add_argument("--plddt", "-p", type=str, help="Path to json with plddt data from structure prediction")
    parser.add_argument("--dist", "-d", type=float, default = 3.5, help="Distance cutoff")
    parser.add_argument("--no-relax", "-r", action="store_true", default=False, help="Prevent rosetta relaxing structures")
    args = parser.parse_args()

    args.input_dir = args.input_dir.rstrip('/') 
    input_list = glob.glob(f'{args.input_dir}/*.pdb')

    # get info abt how many cores avail: 
    npools = min(multiprocessing.cpu_count() - 1, len(input_list))

    # calc what chunk size should be for large input list
    chunksize = math.ceil(len(input_list) / npools) if npools < len(input_list) else 1

    with open(f"{args.input_dir}/searchparams_{date_time_string}.txt", 'w') as fout: 
        fout.write(f"cutoff distance: {args.dist} Angstrom(s)\n") 
        fout.write(f"receptor residues: {receptor_vip_resis}\n")
        fout.write(f"receptor atoms: {receptor_vip_atoms}\n")
        fout.write(f"receptor hydrophobic residues: {receptor_hydrophobic_resis}\n")
        fout.write(f"receptor hydrophobic atoms: {receptor_hydrophobic_atoms}\n")
        fout.write(f"hydrophibic pairs: {vip_hydrophobic}\n") 

    ratings = pd.DataFrame(columns=['design_id',
                                    'peptide_sequence', 
                                    'pre-relax_score',
                                    'post-relax_score',
                                    'average_plddt',
                                    'average_chB_plddt',
                                    'total_num_residue_contacts',
                                    'res_contacts', 
                                    'atom_contacts',
                                    'hydrophobic_contacts'])

    pldtt_dict = dict()
    if os.path.exists(args.plddt): 
        plddt_dict = get_plddts_json(args.plddt)
    else:
        pldd_dict = {'empty': (None, None)}

    if args.no_relax: 
        analysis_pdbs = input_list.copy()
        relaxed_scores = [0 for _ in input_list]

        with multiprocessing.Pool(processes=npools) as pool:
            # multiprocess: score structures and store scores in list
            unrelaxed_scores = pool.map(score_structure, input_list, chunksize=chunksize)

    else: 
        relaxed_dir = f"{args.input_dir}/relaxed"
        if not os.path.isdir(relaxed_dir): 
            os.mkdir(relaxed_dir)

        with multiprocessing.Pool(processes=npools) as pool:
            # multiprocess: relax and score structures
            rosetta_scores = pool.starmap(relax_and_score_structure, 
                                            zip(input_list, 
                                                [relaxed_dir for _ in input_list]), 
                                          chunksize=chunksize)
                
        analysis_pdbs = list()
        unrelaxed_scores = list()
        relaxed_scores = list()
        for scores in rosetta_scores: 
            analysis_pdbs.append(scores[0])
            unrelaxed_scores.append(scores[1])
            relaxed_scores.append(scores[2])

    with multiprocessing.Pool(processes=npools) as pool:
        # multiprocess: count contacts and filter by cutoff distance for all files in analysis_pdbs
        count_results = pool.starmap(count_contacts, 
                                     zip(analysis_pdbs, 
                                         [args.cutoff for _ in analysis_pdbs]),
                                     chunksize=chunksize)
        pass

    ratingsdf = pd.DataFrame(columns=['design_id', 
                                      'peptide_sequence', 
                                      'pre-relax_score', 
                                      'post-relax_score', 
                                      'average_plddt', 
                                      'average_chB_plddt', 
                                      'total_num_residue_contacts', 
                                      'res_contacts', 
                                      'atom_contacts', 
                                      'hydrophobic_contacts'])

    for all_data in zip(input_list, count_results, unrelaxed_scores, relaxed_scores):
        # match all 'design_id' from plddt json
        design_name = os.path.splitext(os.path.basename(all_data[0]))[0] 

        try:
            # pyright thinks plddt_dict may be unbound
            mean_plddt = plddt_dict[design_name][0] # pyright: ignore
            perres_plddt = plddt_dict[design_name][1] # pyright: ignore
            chB_plddt = np.mean(perres_plddt[-8:])
        except:
            mean_plddt, perres_plddt, chB_plddt = 0, 0, 0

        ratingsdf = pd.concat([ratingsdf if not ratingsdf.empty else None, 
                               pd.DataFrame([{'design_id': design_name, 
                                              'peptide_sequence': all_data[1][0], 
                                              'pre-relax_score': all_data[2],
                                              'post-relax_score': all_data[3], 
                                              'average_plddt': mean_plddt, 
                                              'average_chB_plddt': chB_plddt, 
                                              'total_num_residue_contacts': all_data[1][1], 
                                              'res_contacts': all_data[1][2], 
                                              'atom_contacts': all_data[1][3], 
                                              'hydrophobic_contacts': all_data[1][4] }])
                               ], ignore_index=True)
    ratingsdf.to_csv(f"{args.input_dir}/ratings_{date_time_string}.csv")
