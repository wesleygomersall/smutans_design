#!/usr/bin/env python3

from vips import * 

import datetime
import pandas as pd
from pymol import cmd, stored 

date_time_string = datetime.datetime.now().strftime("%Y%m%d-%H%M%S") 

def collect_atom_pairs_AB(pdbfile):
    A_list = []
    B_list = []
    pair_list = []
    
    with open(pdbfile, 'r') as fin: 
        for line in fin: 
            if line.split()[0] == 'ATOM' and line.split()[4] == 'A': # in chain
                A_list.append(('mymodel', line.split()[1]))
            if line.split()[0] == 'ATOM' and line.split()[4] == 'B': 
                B_list.append(('mymodel', line.split()[1]))

    for a in A_list: 
        for b in B_list: 
            pair_list.append((a, b))
    
    return pair_list
    

def pymol_find_contacts(model_path: str,
                        # cutoff_len: float = 3.5
                        ):

    cmd.load(model_path, "mymodel")

    natoms_nresi = [cmd.count_atoms("chain A"),
                    cmd.count_atoms("chain B or chain X"),
                    cmd.count_atoms("chain A and name CA"),
                    cmd.count_atoms("(chain B) and name CA") ]

    cmd.select("protein_chain", "chain A and mymodel")
    cmd.select("peptide_ligand", "not chain A and mymodel")

    pair_list = collect_atom_pairs_AB(model_path)
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
