#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Creates fasta for gene flanks
"""

import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from pathlib import Path


__author__ = "Samuel Lipworth, William Matlock"


# arguments for the script
def get_arguments():
    parser = argparse.ArgumentParser(description = 'flanker',
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group('required arguments')

    # input fasta file
    required.add_argument('-i', '--fasta_file', action = 'store',
                        required = True,
                        help = 'Input fasta file')

    # flanks desired
    parser.add_argument('-f','--flank', action='store',
                        help='Choose which side(s) of the gene to extract (left/right/both)',
                        default='both')

    # window arguments
    parser.add_argument('-w', '--window', action = 'store', type=int,
                        help = 'Length of flanking sequence/first window length',
                        default = 1000)
    parser.add_argument('-wstop', '--window_stop', action='store',type=int,
                        help = 'Final window length',
                        default = None)
    parser.add_argument('-wstep', '--window_step', action='store',type=int,
                        help = 'Step in window sequence',
                        default = None)

    # is sequence circularised?
    parser.add_argument('-circ', '--circ', action = 'store_true',
                        help = 'Is sequence circularised')

    # include gene in output sequence?
    parser.add_argument('-inc', '--include_gene', action = 'store_true',
                        help = 'Include the gene of interest')

    # speciify abricate database
    parser.add_argument('-db', '--database', action = 'store',
                        help = 'Choose Abricate database e.g. NCBI, resfinder',
                        default = 'resfinder')

    # gene(s) to annotate
    gene_group = parser.add_mutually_exclusive_group(required = True)
    gene_group.add_argument('-g', '--goi', action = 'store',
                        help = 'Gene of interest (escape any special characters)')
    gene_group.add_argument('-lg', '--list_of_genes', action= 'store',
                        help = 'Takes a .txt /n list of genes to process')


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    return args


def run_abricate(file):
    args=get_arguments()
    abricate_command = ["abricate", "--db", args.database, file] # shell commands
    p = subprocess.Popen(abricate_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # run abricate
    out, _ = p.communicate() # read stdout data
    out = out.decode() # decode from unicode
    o = open(str(file + '_resfinder'),'w') # create output file
    o.write(out) # write output file
    o.close() # close output


def flank_positions(file, gene_):
    data = pd.read_csv(file, sep='\t', header = 0)
    gene = data[data["GENE"].str.contains(gene_, regex=False)]

    # check if gene is found
    if len(gene) == 0:
        return True

    # gene found
    g = gene['GENE'].iloc[0]

    # LHS flank
    start = int(gene['START'].iloc[0]) # start of gene
    start -= 1 # end of LHS flank

    # RHS flank
    end = int(gene['END'].iloc[0]) # end of gene/start of RHS flank

    return(start, end, g)


def flank_fasta_file_circ(file, window,gene):
    args = get_arguments()

    abricate_file = str(file + '_resfinder') # name of abricate output for fasta

    pos = flank_positions(abricate_file, gene)

    if (pos == True):
        return print(f"Gene {args.goi} not found in {args.fasta_file}")
        
    for record in SeqIO.parse(file, "fasta"):

            print(pos[2] + ' found!')

            w = int(window)
            l = len(record.seq)
            x = args.flank
            
            # if window is too long for sequence length
            if w > 0.5 * (pos[0] - pos[1] + l):
                return print('Window too long for sequence length')

            # if window exceeds sequence length after gene and gene is included
            if (pos[1] + w > l) and (args.include_gene == True) and (args.flank == 'both'):
                record.seq = record.seq[(pos[0]-w):l] + record.seq[0:(pos[1]+w-l)]
                break

            if (pos[1] + w > l) and (args.include_gene == True) and (args.flank == 'left'):
                record.seq = record.seq[(pos[0]-w):pos[1]]
                break

            if (pos[1] + w > l) and (args.include_gene == True) and (args.flank == 'right'):
                record.seq = record.seq[(pos[0]):l] + record.seq[0:(pos[1]+w-l)] # logic needs checking
                break

            # if window exceeds sequence length after gene and gene is excluded
            if (pos[1] + w > l) and (args.include_gene == False) and (args.flank == 'both'):
                record.seq = record.seq[(pos[0]-w):pos[0]] + record.seq[pos[1]:l] + record.seq[0:((pos[1]+w)-l)]
                break

            if (pos[1] + w > l) and (args.include_gene == False) and (args.flank == 'left'):
                record.seq = record.seq[(pos[0]-w):pos[0]]
                break

            if (pos[1] + w > l) and (args.include_gene == False) and (args.flank == 'right'): 
                record.seq = record.seq[pos[1]:l] + record.seq[0:((pos[1]+w)-l)]
                break

            # if window exceeds sequence length before gene and gene is included
            if (pos[0] - w < 0) and (args.include_gene == True) and (args.flank == 'both'):
                record.seq = record.seq[0:(pos[1]+w)] + record.seq[(l-(w-pos[0])):l]
                break

            if (pos[0] - w < 0) and (args.include_gene == True) and (args.flank == 'left'):
                record.seq = record.seq[0:(pos[1])] + record.seq[(l-(w-pos[0])):l] # logic needs checking
                break

            if (pos[0] - w < 0) and (args.include_gene == True) and (args.flank == 'right'):
                record.seq = record.seq[pos[0]:(pos[1]+w)]
                break
            
            # if window exceeds sequence length before gene and gene is excluded
            if (pos[0] - w < 0) and (args.include_gene == False) and (args.flank == 'both'):
                record.seq = record.seq[0:pos[0]] + record.seq[pos[1]:(pos[1]+w)] + record.seq[(l-(w-pos[0])):l]
                break

            if (pos[0] - w < 0) and (args.include_gene == False) and (args.flank == 'left'):
                record.seq = record.seq[0:pos[0]] + record.seq[(l-(w-pos[0])):l]
                break

            if (pos[0] - w < 0) and (args.include_gene == False) and (args.flank == 'right'):
                record.seq =  record.seq[pos[1]:(pos[1]+w)]
                break
            
            # if window length is short enough and gene is included
            if (args.include_gene == True) and (args.flank == 'both'):
                record.seq = record.seq[(pos[0]-w):(pos[1]+w)]
                break

            if (args.include_gene == True) and (args.flank == 'left'):
                record.seq = record.seq[(pos[0]-w):pos[1]]
                break

            if (args.include_gene == True) and (args.flank == 'right'):
                record.seq = record.seq[pos[0]:(pos[1]+w)]
                break

            # if window length is short enough and gene is excluded
            if (args.include_gene == False) and (args.flank == 'both'):
                record.seq = record.seq[(pos[0]-w):pos[0]] + record.seq[pos[1]:(pos[1]+w)]
                break

            if (args.include_gene == False) and (args.flank == 'left'):
                record.seq = record.seq[(pos[0]-w):pos[0]]
                break

            if (args.include_gene == False) and (args.flank == 'right'):
                record.seq = record.seq[pos[1]:(pos[1]+w)]
                break
            
    record.description = f"{record.description} | {pos[2]} | {w}bp window"
    with open(f"{Path(file).stem}_{pos[2]}_{w}_{x}_flank.fasta", "w") as f:
        SeqIO.write(record, f, "fasta")
        print(f"{f.name} sucessfully created!")
        f.close()


def flank_fasta_file_lin(file, window,gene):
    args = get_arguments()
    abricate_file = str(file + '_resfinder') # name of abricate output for fasta

    pos = flank_positions(abricate_file, gene)

    if pos != True:

        for record in SeqIO.parse(file, "fasta"):

            print(pos[2] + ' found')

            w = int(window)
            l = len(record.seq)

            #take both flanks
            if args.flank == 'both':

            #include the gene if desired
                if args.include_gene == True:
                    record.seq = record.seq[max(0,pos[0]-w):min(len(record.seq), pos[1]+w)]

                else:
                    record.seq = record.seq[max(0, pos[0]-w):pos[0]] + record.seq[pos[1]:min(len(record.seq), pos[1]+w)]

                    record.description = f"{record.description} | {pos[2]} | {w}bp window"

                with open(f"{Path(file).stem}_{pos[2]}_{w}_both_flank.fasta", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()

            #or if desired only go left
            elif args.flank == 'left':
                #include the gene if desired
                if args.include_gene == True:
                    record.seq = record.seq[max(0,pos[0]-w):min(len(record.seq),pos[1])]


                else:
                    record.seq = record.seq[max(0, pos[0]-w):pos[0]]

                    record.description = f"{record.description} | {pos[2]} | {w}bp window"

                with open(f"{Path(file).stem}_{pos[2]}_{w}_left_flank.fasta", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()

            #or if desired only go right
            elif args.flank == 'right':
                #include the gene if desired
                if args.include_gene == True:
                    record.seq = record.seq[pos[0]:min(len(record.seq), pos[1]+w)]


                else:
                    record.seq = record.seq[pos[1]:min(len(record.seq), pos[1]+w)]

                    record.description = f"{record.description} | {pos[2]} | {w}bp window"

                with open(f"{Path(file).stem}_{pos[2]}_{w}_right_flank.fasta", "w") as f:
                    SeqIO.write(record, f, "fasta")
                    print(f"{f.name} sucessfully created!")
                    f.close()



def main():
    args = get_arguments()
    run_abricate(args.fasta_file)
    if args.list_of_genes is not None:
        with open(args.list_of_genes) as gene_list:
           for gene in gene_list:
               print("Working on gene {}".format(gene.strip()))
               if args.window_stop is not None:
                   for i in range(args.window, args.window_stop, args.window_step):
                       if args.circ == True:
                           flank_fasta_file_circ(args.fasta_file, i, gene.strip())
                       else:
                           flank_fasta_file_lin(args.fasta_file, i, gene.strip())
               else:
                   if args.circ == True:
                       flank_fasta_file_circ(args.fasta_file, args.window, gene.strip())
                   else:
                       flank_fasta_file_lin(args.fasta_file, args.window,gene.strip())
    else:
        if args.window_stop is not None:
            for i in range(args.window, args.window_stop, args.window_step):
                if args.circ == True:
                    flank_fasta_file_circ(args.fasta_file, i, args.goi)
                else:
                    flank_fasta_file_lin(args.fasta_file, i, args.goi)
        else:
            if args.circ == True:
                flank_fasta_file_circ(args.fasta_file, args.window,args.goi)
            else:
                flank_fasta_file_lin(args.fasta_file, args.window,args.goi)


if __name__ == '__main__':
    main()
