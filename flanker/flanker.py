#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Flanker v1.0
"""
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from pathlib import Path
from cluster import *
from salami import *
import time
import logging as log

start = time.time()

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

    # gene(s) to annotate
    required.add_argument('-g', '--gene', nargs='+', action = 'store', required=True,
                        help = 'Gene of interest (escape any special characters)')

    # flanks desired
    parser.add_argument('-f','--flank', action='store',
                        help='Choose which side(s) of the gene to extract (upstream/downstream/both)',
                        default='both')

    # running mode
    parser.add_argument('-m', '--mode',action='store',
                        help = 'One of "default" - normal mode, "mm" - multi-allelic cluster, or "sm" - salami-mode',
                        default = "default")

    # is sequence circularised?
    parser.add_argument('-circ', '--circ', action = 'store_true',
                        help = 'Is sequence circularised')

    # include gene in output sequence?
    parser.add_argument('-inc', '--include_gene', action = 'store_true',
                        help = 'Include the gene of interest')

    # specify abricate database
    parser.add_argument('-db', '--database', action = 'store',
                        help = 'Choose Abricate database e.g. NCBI, resfinder',
                        default = 'resfinder')

    # output verbosity
    parser.add_argument("-v", "--verbose", const=1, default=0, type=int, nargs="?",
                    help="Increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity.")

    # window arguments
    window=parser.add_argument_group('window options')
    window.add_argument('-w', '--window', action = 'store', type=int,
                        help = 'Length of flanking sequence/first window length',
                        default = 1000)
    window.add_argument('-wstop', '--window_stop', action='store',type=int,
                        help = 'Final window length',
                        default = None)
    window.add_argument('-wstep', '--window_step', action='store',type=int,
                        help = 'Step in window sequence',
                        default = None)

    # clustering options
    cluster=parser.add_argument_group('clustering options')
    cluster.add_argument('-cl','--cluster',help='Turn on clustering mode?',action='store_true')
    cluster.add_argument('-o', '--outfile',action='store',help='Prefix for the clustering file')
    cluster.add_argument('-tr', '--threshold',action='store',help='mash distance threshold for clustering'),
    cluster.add_argument('-p', '--threads',action='store',help='threads for mash to use')

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

# returns the start and end positions of the annotation
def flank_positions(data, gene_):

    gene = data[data["GENE"].str.match(gene_)]

    # check if gene is found
    if len(gene) == 0:
        return True

    # gene foundname=str(recorname=str(record.description)d.description)
    g = gene['GENE'].iloc[0]

    # LHS flank
    start = int(gene['START'].iloc[0]) # start of gene
    start -= 1 # end of LHS flank

    # RHS flank
    end = int(gene['END'].iloc[0]) # end of gene/start of RHS flank

    return(start, end, g)

# writes output fasta
def writer(record, gene, window, isolate, x,gene_sense):
    record.description = f"{record.description} | {gene} | {window}bp window"

    with open(f"{isolate}_{gene}_{window}_{x}_flank.fasta", "w") as f:
        if gene_sense == '+':
            SeqIO.write(record, f, "fasta")
            log.info(f"{f.name} sucessfully created!")
            f.close()
        elif gene_sense == '-':
            record.seq=record.seq.reverse_complement()
            SeqIO.write(record, f, "fasta")
            log.info(f"{f.name} sucessfully created!")
            f.close()

#this function is needed for multi-fasta files
def filter_abricate(data, isolate):

    data = data.loc[data['SEQUENCE'] == isolate]

    return(data)

def flank_fasta_file_circ(file, window,gene):
    args = get_arguments()

    unfiltered_abricate_file = str(file + '_resfinder') # name of abricate output for fasta
    data = pd.read_csv(unfiltered_abricate_file, sep='\t', header = 0)

    guids=data['SEQUENCE'].unique()
    log.debug(guids)
    #can't just use abricate output for whole of muli-fasta
    for guid in guids:

        abricate_file=filter_abricate(data,guid)
        pos = flank_positions(abricate_file, gene)
        if (pos == True):
            log.warning(f"Error: Gene {gene} not found in {guid}")

        else:
            gene_sense=abricate_file.loc[abricate_file['GENE']==gene].filter(items=['STRAND']
                                                                             
            log.info(f"Gene {gene} found in {guid}")

            gene_sense=str(gene_sense['STRAND'].iloc[0])

            log.debug(gene_sense)
            log.debug(pos)

    # initialise dictionaries for sequence splicing functions

            d = {(True, 'both'): lambda record, pos, w, l : record.seq[(pos[0]-w):(pos[1]+w)],
                (True, 'upstream'): lambda record, pos, w, l : record.seq[(pos[0]-w):(pos[1])],
                (True, 'downstream'): lambda record, pos, w, l : record.seq[pos[0]:(pos[1]+w)],
                (False, 'both'): lambda record, pos, w, l : record.seq[(pos[0]-w):pos[0]] + record.seq[pos[1]:(pos[1]+w)],
                (False, 'upstream'): lambda record, pos, w, l : record.seq[(pos[0]-w):pos[0]],
                (False, 'downstream'): lambda record, pos, w, l : record.seq[pos[1]:(pos[1]+w)]}

            d_before = {(True, 'both'): lambda record, pos, w, l : record.seq[0:(pos[1]+w)] + record.seq[(l-(w-pos[0])):l],
                (True, 'upstream'): lambda record, pos, w, l : record.seq[0:(pos[1])] + record.seq[(l-(w-pos[0])):l],
                (True, 'downstream'): lambda record, pos, w, l : record.seq[pos[0]:(pos[1]+w)],
                (False, 'both'): lambda record, pos, w, l : record.seq[0:pos[0]] + record.seq[pos[1]:(pos[1]+w)] + record.seq[(l-(w-pos[0])):l],
                (False, 'upstream'): lambda record, pos, w, l : record.seq[0:pos[0]] + record.seq[(l-(w-pos[0])):l],
                (False, 'downstream'): lambda record, pos, w, l : record.seq[pos[1]:(pos[1]+w)]}

            d_after = {(True, 'both'): lambda record, pos, w, l : record.seq[(pos[0]-w):l] + record.seq[0:(pos[1]+w-l)],
                (True, 'upstream'): lambda record, pos, w, l : record.seq[(pos[0]-w):pos[1]],
                (True, 'downstream'): lambda record, pos, w, l : record.seq[(pos[0]):l] + record.seq[0:(pos[1]+w-l)],
                (False, 'both'): lambda record, pos, w, l : record.seq[(pos[0]-w):pos[0]] + record.seq[pos[1]:l] + record.seq[0:((pos[1]+w)-l)],
                (False, 'upstream'): lambda record, pos, w, l : record.seq[(pos[0]-w):pos[0]],
                (False, 'downstream'): lambda record, pos, w, l : record.seq[pos[1]:l] + record.seq[0:((pos[1]+w)-l)]}

        # loop through records in fasta
            for record in SeqIO.parse(file, "fasta"):
                #select the fasta record of interest
                if record.description == guid:

                    name=str(record.description)

                    log.info(pos[2] + ' found!')

                    w = int(window)
                    l = len(record.seq)
                    x = args.flank

                # if window is too long for sequence length
                    if w > 0.5 * (pos[0] - pos[1] + l):
                        log.warning(f"Error: Window length {w} too long for sequence length {l}")
                        continue

                # if window exceeds sequence length after gene

                    if (pos[1] + w > l):
                        log.debug("Window exceeds seq length after gene")
                        record.seq = d_after[(args.include_gene, args.flank)](record, pos, w, l)
                        writer(record, pos[2], w, guid, x, gene_sense)
                        continue

                # if window exceeds sequence length before gene

                    if (pos[0] - w < 0):
                        log.debug("Window excees seq length before gene")
                        record.seq = d_before[(args.include_gene, args.flank)](record, pos, w, l)
                        writer(record, pos[2], w, guid, x, gene_sense)
                        continue

                    else:
                        log.debug("Window is all good")

                        record.seq = d[(args.include_gene, args.flank)](record, pos, w, l)
                        writer(record, pos[2], w, guid, x, gene_sense)
                        continue


def flank_fasta_file_lin(file, window,gene):
    args = get_arguments()
    unfiltered_abricate_file = str(file + '_resfinder') # name of abricate output for fasta
    data = pd.read_csv(unfiltered_abricate_file, sep='\t', header = 0)

    guids=data['SEQUENCE'].unique()

    for guid in guids:
        abricate_file=filter_abricate(data,guid)
        pos = flank_positions(abricate_file, gene)
        if pos == True:
            log.error(f"Error: Gene {gene} not found in {guid}")

        else:
             gene_sense=abricate_file.loc[abricate_file['GENE']==gene].filter(items=['STRAND'])
             gene_sense=str(gene_sense['STRAND'].iloc[0])
                                                                             
             d_lin = {(True, 'both'): lambda record, pos, w, l: record.seq[max(0,pos[0]-w):min(l, pos[1]+w)],
             (True, 'upstream'): lambda record, pos, w, l : record.seq[max(0,pos[0]-w):min(l,pos[1])],
             (True, 'downstream'): lambda record, pos, w, l : record.seq[pos[0]:min(l, pos[1]+w)],
             (False, 'both'): lambda record, pos, w, l : record.seq[max(0, pos[0]-w):pos[0]] + record.seq[pos[1]:min(l, pos[1]+w)],
             (False, 'upstream'): lambda record, pos, w, l : record.seq[max(0, pos[0]-w):pos[0]],
             (False, 'downstream'): lambda record, pos, w, l : record.seq[pos[1]:min(l, pos[1]+w)]}

             w = int(window)
             x = args.flank

             for record in SeqIO.parse(file, "fasta"):
                 if record.description == guid:
                     name=str(record.description)

                     log.info(f"{gene} found in {record.description}")

                     l = len(record.seq)

                     record.seq = d_lin[(args.include_gene, args.flank)](record, pos, w, l)
                     writer(record, pos[2], w, guid, x,gene_sense)
                     continue


def flanker_main():
    args = get_arguments()

    run_abricate(args.fasta_file)

    gene_list=args.gene
    log.debug(gene_list)


    if args.window_stop is not None:
        for i in range(args.window, args.window_stop, args.window_step):
            for gene in gene_list:

                if args.circ == True:
                    flank_fasta_file_circ(args.fasta_file, i, gene.strip())
                else:
                    flank_fasta_file_lin(args.fasta_file, i, gene.strip())

                if args.cluster ==True and args.mode =='default':

                    define_clusters(gene,i,args.threads,args.threshold,args.outfile)
                    flank_scrub()

            if args.cluster==True and args.mode=='mm':

                define_clusters(gene,i,args.threads,args.threshold,args.outfile)
                log.info("Cleaning up")
                flank_scrub()

    else:
        for gene in gene_list:
            if args.circ == True:
                flank_fasta_file_circ(args.fasta_file, args.window, gene.strip())
            else:

                flank_fasta_file_lin(args.fasta_file, args.window,gene.strip())
            if args.cluster ==True and args.mode =='default':
                log.info("Performing clustering")
                define_clusters(gene,i,args.threshold,args.outfile)
                log.info("Cleaning up")
                flank_scrub()

        if args.cluster==True and args.mode=='mm':
            log.info("Performing clustering")
            define_clusters(gene,"mm",args.threads,args.threshold,args.outfile)
            log.info("Cleaning up")
            flank_scrub()

def main():
    args=get_arguments()

    logger = log.getLogger()

    log.basicConfig(format="%(message)s")

    if args.verbose == 0:
        logger.setLevel(log.WARNING)
    elif args.verbose == 1:
        logger.setLevel(log.INFO)
    elif args.verbose == 2:
        logger.setLevel(log.DEBUG)

    log.info(args)
    if args.mode =="default" or args.mode == "mm":
        flanker_main()
    elif args.mode =="sm":
        salami_main(args.gene,args.fasta_file,args.window,args.window_step,args.window_stop,args.outfile,args.threshold,args.cluster)

    end = time.time()
    log.info(f"All done in {round(end - start, 2)} seconds")


if __name__ == '__main__':
    main()
