#!/usr/bin/env python3

import sys
import time
import argparse
import subprocess
import logging as log
from pathlib import Path

import pandas as pd
from Bio import SeqIO

import jellyfish as jf

from flanker import cluster, salami


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
    genes=parser.add_mutually_exclusive_group(required=True)
    genes.add_argument('-g', '--gene', nargs='+', action = 'store',
                        help = 'Gene(s) of interest (escape any special characters). Use space seperation for multipe genes')
    genes.add_argument('-log','--list_of_genes',action='store',default=False,
                        help = 'Line separated file containing genes of interest')

    # closest match mode
    parser.add_argument('-cm', '--closest_match', action = 'store_true',
                        help = 'Find closest match to query')

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
                        default = 'ncbi')

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
    cluster.add_argument('-cl','--cluster',help='Turn on clustering mode?',action='store_true'),
    cluster.add_argument('-o', '--outfile',action='store',help='Prefix for the clustering file',default='out'),
    cluster.add_argument('-tr', '--threshold',action='store',help='mash distance threshold for clustering',default="0.001"),
    cluster.add_argument('-p', '--threads',action='store',help='threads for mash to use',default='1'),
    cluster.add_argument('-k', '--kmer_length',action='store', help='kmer length for Mash',default='21'),
    cluster.add_argument('-s', '--sketch_size',action='store', help='sketch size for mash',default='1000')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    return args

# validate input
def check_input(fasta_file):
    fasta_records = list(SeqIO.parse(fasta_file, 'fasta'))
    assert len(fasta_records) >= 1, 'No records found in fasta file'

# annotate for gene(s) of interest
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
    args = get_arguments()
    
    if args.closest_match == False:
        gene = data[data["GENE"].str.match(gene_)]
        # check if gene is found
        if len(gene) == 0:
            return True
    else:
        data["dist"] = [jf.levenshtein_distance(gene_, x) for x in data["GENE"]]
        gene = data.sort_values(by="dist", ascending=True)

    print(gene)

    g = gene['GENE'].iloc[0]

    print(type(g))

    # LHS flank
    start = int(gene['START'].iloc[0]) # start of gene
    start -= 1 # end of LHS flank

    # RHS flank
    end = int(gene['END'].iloc[0]) # end of gene/start of RHS flank

    return(start, end, g)

# writes output fasta
def writer(record, gene, window, isolate, x,gene_sense):
    record.description = f"{record.description} | {gene} | {window}bp window"

    (gene,window,isolate,x)

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

# for processing multi-fasta files
def filter_abricate(data, isolate):

    data = data.loc[data['SEQUENCE'] == isolate]

    return(data)

# generates flanks for circularised sequences
def flank_fasta_file_circ(file, window,gene):
    args = get_arguments()

    unfiltered_abricate_file = file + '_resfinder' # name of abricate output for fasta
    data = pd.read_csv(unfiltered_abricate_file, sep='\t', header = 0)

    guids=data['SEQUENCE'].unique()
    log.debug(guids)

    for guid in guids:

        abricate_file=filter_abricate(data,guid)
        pos = flank_positions(abricate_file, gene)

        if (pos == True):
            log.warning(f"Error: Gene {gene} not found in {guid}")

        else:
            pos=list(pos)
            gene_sense=abricate_file.loc[abricate_file['GENE']==gene].filter(items=['STRAND'])

            log.info(f"Gene {gene} found in {guid}")

            gene_sense=str(gene_sense['STRAND'].iloc[0])

            log.debug(gene_sense)
            log.debug(pos)

            # initialise dictionaries of sequence splicing functions

            d = {(True, 'both'): lambda record, pos, w, l : record.seq[(pos[0]-w):(pos[1]+w)],
                (True, 'upstream'): lambda record, pos, w, l : record.seq[(pos[0]-w):(pos[1])],
                (True, 'downstream'): lambda record, pos, w, l : record.seq[pos[0]:(pos[1]+w)],
                (False, 'both'): lambda record, pos, w, l : record.seq[(pos[0]-w):pos[0]] + record.seq[pos[1]:(pos[1]+w)],
                (False, 'upstream'): lambda record, pos, w, l : record.seq[(pos[0]-w):pos[0]],
                (False, 'downstream'): lambda record, pos, w, l : record.seq[pos[1]:(pos[1]+w)]}

            d_before = {(True, 'both'): lambda record, pos, w, l : record.seq[(l-(w-pos[0])):l] + record.seq[0:(pos[1]+w)],
                (True, 'upstream'): lambda record, pos, w, l : record.seq[(l-(w-pos[0])):l] + record.seq[0:(pos[1])] ,
                (True, 'downstream'): lambda record, pos, w, l : record.seq[pos[0]:(pos[1]+w)],
                (False, 'both'): lambda record, pos, w, l : record.seq[(l-(w-pos[0])):l] + record.seq[0:pos[0]] + record.seq[pos[1]:(pos[1]+w)] ,
                (False, 'upstream'): lambda record, pos, w, l : record.seq[(l-(w-pos[0])):l] + record.seq[0:pos[0]],
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
                w = int(window)
                l = len(record.seq)
                x = args.flank
                if record.description == guid:
                    if gene_sense == '-':

                        #record.seq = record.seq.reverse_complement()
                        if args.flank == 'upstream':
                            x = 'downstream'
                        else:
                            x = 'upstream'

                    name=record.description

                    log.info(pos[2] + ' found!')

                    # if window is too long for sequence length
                    if w > 0.5 * (pos[0] - pos[1] + l):
                        log.info(f"Error: Window length {w} too long for sequence length {l}")
                        continue

                    # if window exceeds sequence length after gene
                    if (pos[1] + w > l):
                        log.info("Window exceeds seq length after gene")
                        record.seq = d_after[(args.include_gene, x)](record, pos, w, l)
                        writer(record, pos[2], w, guid, args.flank, gene_sense)
                        continue

                    # if window exceeds sequence length before gene
                    if (pos[0] - w < 0):
                        log.info("Window exceeds seq length before gene")
                        record.seq = d_before[(args.include_gene, x)](record, pos, w, l)
                        writer(record, pos[2], w, guid, args.flank, gene_sense)
                        continue

                    else:
                        log.debug("Window is good")

                        record.seq = d[(args.include_gene, x)](record, pos, w, l)
                        writer(record, pos[2], w, guid, args.flank, gene_sense)
                        continue

# generates flanks for linear sequences
def flank_fasta_file_lin(file, window, gene):
    args = get_arguments()
    unfiltered_abricate_file = file + '_resfinder' # name of abricate output for fasta
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

             # initialise dictionary of sequence splicing functions

             d_lin = {(True, 'both'): lambda record, pos, w, l: record.seq[max(0,pos[0]-w):min(l, pos[1]+w)],
             (True, 'upstream'): lambda record, pos, w, l : record.seq[max(0,pos[0]-w):min(l,pos[1])],
             (True, 'downstream'): lambda record, pos, w, l : record.seq[pos[0]:min(l, pos[1]+w)],
             (False, 'both'): lambda record, pos, w, l : record.seq[max(0, pos[0]-w):pos[0]] + record.seq[pos[1]:min(l, pos[1]+w)],
             (False, 'upstream'): lambda record, pos, w, l : record.seq[max(0, pos[0]-w):pos[0]],
             (False, 'downstream'): lambda record, pos, w, l : record.seq[pos[1]:min(l, pos[1]+w)]}

             w = int(window)
             x = args.flank

             # loop through records in fasta
             for record in SeqIO.parse(file, "fasta"):
                 if record.description == guid:
                     if gene_sense == '-':

                         if args.flank == 'upstream':
                             x = 'downstream'
                         else:
                             x = 'upstream'
                     name=record.description

                     log.info(f"{gene} found in {record.description}")

                     l = len(record.seq)

                     record.seq = d_lin[(args.include_gene, x)](record, pos, w, l)
                     writer(record, pos[2], w, guid, args.flank, gene_sense)

                     continue

def flanker_main():
    args = get_arguments()


    run_abricate(args.fasta_file)

    if args.list_of_genes == False:
        gene_list=args.gene

    else:
        gene_list=[]
        with open(args.list_of_genes, 'rb') as gl:
            for line in gl:
                line=line.decode('utf-8')
                gene_list.append(line.strip())

    log.debug(gene_list)

    if args.window_stop is not None:
        for i in range(args.window, args.window_stop, args.window_step):
            for gene in gene_list:

                if args.circ == True:
                    flank_fasta_file_circ(args.fasta_file, i, gene.strip())
                else:
                    flank_fasta_file_lin(args.fasta_file, i, gene.strip())

                if args.cluster ==True and args.mode =='default':

                    cluster.define_clusters(gene,i,args.threads,args.threshold,args.outfile,args.kmer_length,args.sketch_size)
                    cluster.flank_scrub()

            if args.cluster==True and args.mode=='mm':

                cluster.define_clusters(gene,i,args.threads,args.threshold,args.outfile,args.kmer_length,args.sketch_size)
                log.info("Cleaning up")
                cluster.flank_scrub()

    else:
        for gene in gene_list:
            if args.circ == True:
                flank_fasta_file_circ(args.fasta_file, args.window, gene.strip())
            else:

                flank_fasta_file_lin(args.fasta_file, args.window,gene.strip())
            if args.cluster ==True and args.mode =='default':
                log.info("Performing clustering")
                cluster.define_clusters(gene,args.window,args.threads,args.threshold,args.outfile,args.kmer_length,args.sketch_size)
                log.info("Cleaning up")
                cluster.flank_scrub()

        if args.cluster==True and args.mode=='mm':
            log.info("Performing clustering")
            cluster.define_clusters(gene,"mm",args.threads,args.threshold,args.outfile,args.kmer_length,args.sketch_size)
            log.info("Cleaning up")
            cluster.flank_scrub()

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


    if check_input(args.fasta_file):
        log.info(f"{args.fasta_file} is valid and not empty")


    if args.mode =="default" or args.mode == "mm":
        flanker_main()
    elif args.mode =="sm":

        if args.list_of_genes != False:
            salami.salami_main(args.list_of_genes,args.fasta_file,args.include_gene,args.window_step,args.window_stop,args.outfile,args.flank,args.threads,args.threshold,args.cluster,args.kmer_length,args.sketch_size)
        if args.list_of_genes == False:
            salami.salami_main(args.gene,args.fasta_file,args.include_gene,args.window_step,args.window_stop,args.outfile,args.flank,args.threads,args.threshold,args.cluster,args.kmer_length,args.sketch_size)
    end = time.time()
    log.info(f"All done in {round(end - start, 2)} seconds")


if __name__ == '__main__':
    main()
