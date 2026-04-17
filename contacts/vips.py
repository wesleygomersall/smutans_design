#!/usr/bin/env python3

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
