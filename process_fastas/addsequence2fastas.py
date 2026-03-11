#!/usr/bin/env python3

import argparse
import numpy as np


            
prepend = "MNEALMILSNGLLTYLTVLFLLFLFSKVSNVTLSKKELTLFSISNFLIMIAVTMVNVNLFYPAEPLYFIALSIYLNRQNSLSLNIFYGLLPVASSDLFRRAIIFFILDGTQGIVMGSSIITTYMIEFAGIALSYLFLSVFNVDIGRLKDSLTKMKVKKRLIPMNITMLLYYLLIQVLYVIESYNVIPTLKFRKFVVIVYLILFLILISFLSQYTKQKVQNEIMAQKEAQI:"
append = ""

def main(args):

    with open('seqprepended.fasta', 'w') as fout:
        with open(args.input, 'r') as fin: 
            while True: 
                name = fin.readline()
                sequence = fin.readline() 
                if name == '' or sequence == '':
                    break
                else: 
                    fout.write(name)
                    newsequence = prepend + sequence.strip() + append + '\n'
                    fout.write(newsequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, default = "", help="fasta file to filter")
    args = parser.parse_args()
    main(args)
