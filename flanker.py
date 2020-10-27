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
                        help = 'sequence is circularised'),
    parser.add_argument('-i', '--include_gene', action = 'store_true',
                        help = 'include the gene of interest')
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
    gene = data[data["GENE"].str.contains(gene_)]

    # check if gene is found
    if len(gene) != 0:

        # gene found
        g = gene['GENE'].iloc[0]

        # LHS flank
        start = int(gene['START'].iloc[0]) # start of gene
        start -= 1 # end of LHS flank

        # RHS flank
        end = int(gene['END'].iloc[0]) # end of gene/start of RHS flank

        return(start, end, g)

    else:
        return True

    
def flank_fasta_file_circ(file):
    args = get_arguments() 

    abricate_file = str(file + '_resfinder') # name of abricate output for fasta

    pos = flank_positions(abricate_file, args.goi)

    if pos != True:

        for record in SeqIO.parse(file, "fasta"):

            print(pos[2] + ' found!')

            w = int(args.window)
            l = len(record.seq)

            # if window is too long for sequence length
            if w > 0.5 * (pos[0] - pos[1] + l):
                print('Window too long for sequence length')

            # if window exceeds sequence length after gene
            elif pos[1] + w > l:
                
                #include the gene if desired
                if args.include_gene == True:
                    record.seq = record.seq[(pos[0]-w):(pos[1]+w)] + record.seq[pos[1]:l]
                else:
                # loop to start
                    record.seq = record.seq[(pos[0]-w):pos[0]] + record.seq[pos[1]:l] + record.seq[0:((pos[1]+w)-l)] 

                record.description = f"{record.description} | {pos[2]} | {w}bp window"

                with open(f"{Path(file).stem}_{pos[2]}_flank.fasta_{w}", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()

            # if window exceeds sequence length before gene
            elif pos[0] - w < 0:

                #include the gene if desired
                if args.include_gene == True:
                    record.seq = record.seq[0:pos[0]:(pos[1]+w)]
                else:

                # loop to end
                    record.seq = record.seq[0:pos[0]] + record.seq[pos[1]:(pos[1]+w)] + record.seq[(l-(w-pos[0])):l]

                record.description = f"{record.description} | {pos[2]} | {w}bp window"

                with open(f"{Path(file).stem}_{pos[2]}_{w}_flank.fasta", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()

            else:
            
                #include the gene if desired
                if args.include_gene == True:
                    record.seq = record.seq[(pos[0]-w):(pos[1]+w)]
                else:

                    record.seq = record.seq[(pos[0]-w):pos[0]] + record.seq[pos[1]:(pos[1]+w)]

                record.description = f"{record.description} | {pos[2]} | {w}bp window"

                with open(f"{Path(file).stem}_{pos[2]}_{w}_flank.fasta", "w") as f:
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

            print(pos[2] + ' found')

            w = int(args.window)
            l = len(record.seq)

            #include the gene if desired
            if args.include_gene == True:
                record.seq = record.seq[max(0,pos[0]-w):min(len(record.seq), pos[1]+w)]
            else:

                record.seq = record.seq[max(0, pos[0]-w):pos[0]] + record.seq[pos[1]:min(len(record.seq), pos[1]+w)]

            record.description = f"{record.description} | {pos[2]} | {w}bp window"



            with open(f"{Path(file).stem}_{pos[2]}_{w}_flank.fasta", "w") as f:
                SeqIO.write(record, f, "fasta")
                print(f"{f.name} sucessfully created!")
                f.close()

    else:
        print('Gene not found')
          
def main():
    args = get_arguments()
    run_abricate(args.fasta_file)
    if args.circ == True:
        flank_fasta_file_circ(args.fasta_file)
    else:
        flank_fasta_file_lin(args.fasta_file)
    
    

if __name__ == '__main__':
    main()

