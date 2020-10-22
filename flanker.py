#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Creates fasta for gene flanks
"""

import pandas as pd
import numpy as np
import subprocess
import argparse
from Bio import SeqIO
from pathlib import Path


__author__ = "Samuel Lipworth, William Matlock"


def get_arguments():
    parser = argparse.ArgumentParser(description = 'flanker',
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta_file', action = 'store',
                        help = 'fasta file')
    parser.add_argument('-g', '--goi', action = 'store',
                        help = 'gene of interest')
    parser.add_argument('-w', '--window', action = 'store',
                        help = 'length of flanking sequences')
    parser.add_argument('-c', '--circ', action = 'store_true',
                        help = 'sequence is circularised')
    parser
    return parser.parse_args()


def run_abricate(file):
    abricate_command = ["abricate", "--db", "resfinder", file] # shell commands
    p = subprocess.Popen(abricate_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run abricate
    out, _ = p.communicate() # read stdout data 
    out = out.decode() # decode from unicode
    o = open(str(file + '_resfinder'),'w') # create output file
    o.write(out) # write output file
    o.close() # close output

    
def flank_positions(file, gene_):
    args = get_arguments()

    data = pd.read_csv(file, sep='\t', header = 0)
    gene = data.query('GENE == @gene_')

    # check if gene is found
    if len(gene) != 0:

        # LHS flank
        lhs_end = int(gene['START'].iloc[0]) # start of gene
        lhs_end -= 1 # end of LHS flank
        w = int(args.window)
        lhs_start = lhs_end - w # start of LHS flank

        # RHS flank
        rhs_start = int(gene['END'].iloc[0]) # end of gene/start of RHS flank
        rhs_end = rhs_start + w # end of RHS flank

        return(lhs_start, lhs_end, rhs_start, rhs_end)

    else:
        return True

    
def flank_fasta_file_circ(file):
    args = get_arguments() 

    abricate_file = str(file + '_resfinder') # name of abricate output for fasta

    pos = flank_positions(abricate_file, args.goi)

    if pos != True:

        for record in SeqIO.parse(file, "fasta"):

            w = int(args.window)
            l = len(record.seq)

            # if window is too long for sequence length
            if (w - pos[1] > l - pos[3]) | (w - (l - pos[2]) > pos[0]):
                print('Window too long for sequence length')

            # if window exceeds sequence length after gene
            elif pos[2] + w > l:

                # loop to start
                record.seq = record.seq[pos[0]:pos[1]] + record.seq[pos[2]:l] + record.seq[0:(w-(l-pos[2]))] 

                record.description = f"{record.description} | {args.goi} | {w}bp window"

                with open(f"{Path(file).stem}_{args.goi}_flank.fasta", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()

            # if window exceeds sequence length before gene
            elif pos[1] - w < 0:

                # loop to end
                record.seq = record.seq[(l-(w-pos[1])):l] + record.seq[0:pos[1]] + record.seq[pos[2]:pos[3]]

                record.description = f"{record.description} | {args.goi} | {w}bp window"

                with open(f"{Path(file).stem}_{args.goi}_flank.fasta", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()

            else:
            
                record.seq = record.seq[pos[0]:pos[1]] + record.seq[pos[2]:pos[3]]

                record.description = f"{record.description} | {args.goi} | {w}bp window"

                with open(f"{Path(file).stem}_{args.goi}_flank.fasta", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()

    else:
        print('Gene not found')


def flank_fasta_file_lin(file):
    args = get_arguments() 

    abricate_file = str(file + '_resfinder') # name of abricate output for fasta

    pos = flank_positions(abricate_file, args.goi)

    if pos != True:

        for record in SeqIO.parse(file, "fasta"):

            w = int(args.window)
            l = len(record.seq)

            record.seq = record.seq[min(0, pos[0]):pos[1]] + record.seq[pos[2]:max(len(record.seq), pos[3])]

            record.description = f"{record.description} | {args.goi} | {w}bp window"

            with open(f"{Path(file).stem}_{args.goi}_flank.fasta", "w") as f:
                SeqIO.write(record, f, "fasta")
                print(f"{f.name} sucessfully created!")
                f.close()
          
def main():
    args = get_arguments()
    run_abricate(args.fasta_file)
    if args.circ == True:
        flank_fasta_file_circ(args.fasta_file)
    else:
        flank_fasta_file_lin(args.fasta_file)
    
    

if __name__ == '__main__':
    main()

