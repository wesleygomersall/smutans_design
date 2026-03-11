#!/usr/bin/env python3

import re
import argparse

def main(args):
    
    pattern = r"^.F..PFF." # get with arg? 

    with open('filtered.fasta', 'w') as fout:
        with open(args.input, 'r') as fin: 
            while True: 
                name = fin.readline()
                sequence = fin.readline()
                if name == '' or sequence == '':
                    break
                elif re.match(pattern, sequence):
                    print(f"{sequence} matches {pattern}")
                    fout.write(name)
                    fout.write(sequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, default = "", help="fasta file to filter")
    # parser.add_argument("--reg", "-r", type=str, default = "", help="")
    args = parser.parse_args()
    main(args)
