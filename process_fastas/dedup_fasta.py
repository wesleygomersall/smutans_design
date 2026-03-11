#!/usr/bin/env python3

import argparse
import numpy as np

def main(args):

    seqs = dict()
    
    with open(args.input, 'r') as fin: 
        while True: 
            name = fin.readline()
            sequence = fin.readline() 
            if name == '' or sequence == '':
                break
            elif sequence not in seqs.keys():
                seqs.setdefault(sequence, [name])
            else: 
                seqs[sequence] = seqs[sequence] + [name]

    with open('deduped.fasta', 'w') as fout:
        for sequence in seqs.keys(): 
            noccurs = len(seqs[sequence]) 

            # average the scores and add meanscore to newname
            scores = []
            for name in seqs[sequence]: 
                scorestr = name.split()[2] # "score=3.14,"
                scoreflt = float(scorestr.lstrip('score=').strip(','))
                scores.append(scoreflt)
            
            newstring = f", ndups={noccurs}, meanscore={np.mean(scores)}\n"
            newname = seqs[sequence][0].strip() + newstring

            fout.write(newname)
            fout.write(sequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, default = "", help="fasta file to filter")
    args = parser.parse_args()
    main(args)
