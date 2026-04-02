#!/usr/bin/env python3

import os
import glob
import json
import datetime
import argparse
import pandas as pd
import numpy as np
from pymol import cmd, stored 
from pyrosetta import * 
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.relax import FastRelax

init("-mute all")
scorefxn = create_score_function("ref2015_cart")
relax = FastRelax() 
relax.set_scorefxn(scorefxn)

date_time_string = datetime.datetime.now().strftime("%Y%m%d-%H%M%S") 

receptor_vip_resis = {'MET6', 
                      'MET54', 
                      'VAL57', 
                      'ARG99', 
                      'ARG100', 
                      'GLY112',
                      'ILE119',
                      'THR122', 
                      }

receptor_hydrophobic_resis = {'MET6',
                              'MET54',
                              'ILE119',}

receptor_vip_atoms = {894, 802, 791,
                      941, 961, 935,
                      425, 446, 444, 
                      43, 45, 46, 
                      24, 54, 74, 
                      75, 
                      }
receptor_hydrophobic_atoms = {}

hydrophobic_resis = {'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 
                     'PRO', 'PHE', 'MET', 'TRP'}

# all combinations to count hydrophobic contacts   
vip_hydrophobic = {f"{r1}{r2}" for r1 in receptor_hydrophobic_resis for r2 in hydrophobic_resis}

def pymol_find_contacts(model_path: str,
                        cutoff_len: float = 3.5):

    cmd.load(model_path, "mymodel")

    natoms_nresi = [cmd.count_atoms("chain A"),
                    cmd.count_atoms("not chain A"),
                    cmd.count_atoms("chain A and name CA"),
                    cmd.count_atoms("not chain A and name CA") ]

    cmd.select("protein_chain", f"chain A and mymodel")
    cmd.select("peptide_ligand", f"not chain A and mymodel")

    pair_list = cmd.find_pairs("protein_chain", "peptide_ligand", 
                                state1=1,
                                state2=1,
                                cutoff=cutoff_len,
                                mode=1,
                                angle=45)

    stored.res1list = []
    stored.res2list = []
    distances = [] 
    for pair in pair_list: # Get residue info for each atom 
        cmd.select("atom1", f"id {pair[0][1]}")
        cmd.iterate("atom1", "stored.res1list.append((name, resn, resi))")
        cmd.select("atom2", f"id {pair[1][1]}")
        cmd.iterate("atom2", "stored.res2list.append((name, resn, resi))")
        distances.append(cmd.get_distance("atom1", "atom2"))

    df = pd.DataFrame({"chainA_atomID": [pair[0][1] for pair in pair_list],
                       "chainA_name": [x[0] for x in stored.res1list],
                       "chainA_resn": [x[1] for x in stored.res1list],
                       "chainA_resi": [x[2] for x in stored.res1list],
                       "chainB_atomID": [pair[1][1] for pair in pair_list],
                       "chainB_name": [x[0] for x in stored.res2list],
                       "chainB_resn": [x[1] for x in stored.res2list],
                       "chainB_resi": [x[2] for x in stored.res2list],
                       "distance": distances
                       })

    cmd.delete("mymodel")

    return df, natoms_nresi

def pymol_get_chainB_sequence(model_path: str):
    cmd.load(model_path, "mymodel")
    fasta = cmd.get_fastastr('mymodel and chain B')
    peptide_sequence = "".join(fasta.split('\n')[1:])
    cmd.delete("mymodel")
    return peptide_sequence

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input-dir", "-i", type=str, help="Path to dir containting input pdbs")
    parser.add_argument("--plddt", "-p", type=str, help="Path to json with plddt data from structure prediction")
    args = parser.parse_args()

    args.input_dir = args.input_dir.rstrip('/') 

    with open(f"{args.input_dir}/searchparams_{date_time_string}.txt", 'w') as fout: 
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
                                    'res_contacts', 
                                    'atom_contacts',
                                    'hydrophobic_contacts'])

    # read in json objects from args.plddt
    # store in dict: "design_id": (mean_plddt, residue_plddts) ## sep dicts? 
    plddt_lookup = dict()
    with open(args.plddt, 'r') as plddt_json: 
        for line in plddt_json: 
            data = json.loads(line.strip())
            plddt_lookup.setdefault(data['design_id'], 
                                    (data['average_plddt'], 
                                     data['residue_plddt']))

    relaxed_dir = f"{args.input_dir}/relaxed"
    if not os.path.isdir(relaxed_dir): 
        os.mkdir(relaxed_dir)

    for pdbfilename in glob.glob(f'{args.input_dir}/*.pdb'):
        design_name = os.path.splitext(os.path.basename(pdbfilename))[0] # match all 'design_id' from plddt json

        unrelaxed_pose = pose_from_pdb(pdbfilename)
        unrelaxed_score = scorefxn(unrelaxed_pose) 

        # look for relaxed, if not exist use rosetta to create
        relaxed_path = f"{relaxed_dir}/{design_name}_relaxed.pdb"
        if not os.path.exists(relaxed_path): 
            print(f"relaxing design \"{design_name}\"")
            relaxed_pose = unrelaxed_pose.clone()
            relax.apply(relaxed_pose) 
            relaxed_pose.dump_pdb(relaxed_path)
        else: 
            print(f"relaxed design \"{design_name}\" found at {relaxed_path}")
            relaxed_pose = pose_from_pdb(relaxed_path)

        relaxed_score = scorefxn(relaxed_pose) 

        mean_plddt = plddt_lookup[design_name][0]
        perres_plddt = plddt_lookup[design_name][1]
        chB_plddt = np.mean(perres_plddt[-8:]) # TODO verify this is right

        try: 
            contactsdf, _ = pymol_find_contacts(relaxed_path, 3.5)
            # create a set of the chain A contact resis and atoms from contactsdf
            chA_contactresis = set(contactsdf['chainA_resn'] + contactsdf['chainA_resi']) # 3 letter code + number in seq
            chA_contactatoms = set(int(contactsdf['chainA_atomID']))
            # for hydrophobic pair identification
            chAchB_contactpairs = set(contactsdf['chainA_resn'] + contactsdf['chainA_resi'] + contactsdf['chainB_resn']) # 3 letter code + number in seq + _peptideres
        except: # leave empty
            contactsdf = pd.DataFrame()
            chA_contactresis = set() 
            chA_contactatoms = set() 
            chAchB_contactpairs = set()

        num_res_contacts_seen = len(chA_contactresis.intersection(receptor_vip_resis))
        num_atom_contacts_seen = len(chA_contactatoms.intersection(receptor_vip_atoms))
        num_hydrophobic_seen = len(chAchB_contactpairs.intersection(vip_hydrophobic))

        ratings = pd.concat([ratings if not ratings.empty else None, 
                             pd.DataFrame([{'design_id': design_name,
                                            'peptide_sequence': pymol_get_chainB_sequence(relaxed_path),
                                            'pre-relax_score': unrelaxed_score,
                                            'post-relax_score': relaxed_score,
                                            'average_plddt': mean_plddt,
                                            'average_chB_plddt': chB_plddt,
                                            'res_contacts': num_res_contacts_seen, 
                                            'atom_contacts': num_atom_contacts_seen,
                                            'hydrophobic_contacts': num_hydrophobic_seen}])
                            ], 
                            ignore_index=True)

    ratings.to_csv(f"{args.input_dir}/ratings{date_time_string}.csv")

if __name__ == "__main__":
    main()
