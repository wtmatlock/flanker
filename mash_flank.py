#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 17:20:36 2020

@author: sam
"""

import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import os
import argparse
import sys
import subprocess
import datetime
import shutil
import itertools
from Bio import SeqIO

def get_arguments():
    parser = argparse.ArgumentParser(description='mash_flank',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f1','--fasta_file1',action='store',
                        help='first fasta file')
    parser.add_argument('-f2','--fasta_file2',action='store',
                        help='second fasta file')
    parser.add_argument('-g','--goi',action='store',
                        help='gene of interest')
    parser.add_argument('-o','--out',action='store',
                        help='prefix of output file')
    parser.add_argument('-w','--window',action='store',
                        help='length of flanking sequences')
    parser
    return parser.parse_args()


def run_abricate(file):
    abricate_command=["abricate","--db","resfinder",file]
    p=subprocess.Popen(abricate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    out=out.decode()
    o=open(str(file +"_resfinder"),'w')
    o.write(out)
    o.close()
    

def load_start_position(file, gene_):
    args = get_arguments()
    data=pd.read_csv(file,sep='\t',header=0)
    gene=data.query('GENE == @gene_')
    start=gene[['START']]
    start-=1
    w=int(args.window)
    start-= w
    return(start)
    
def load_end_position(file, gene_):
    args = get_arguments()
    data=pd.read_csv(file,sep='\t',header=0)
    gene=data.query('GENE == @gene_')
    start=gene[['END']]
    w=int(args.window)
    start+=w
    return(start)
    
def parse_fasta_file(file):
    args = get_arguments()
    abricate_file=str(args.fasta_file1 + '_resfinder')
    start=load_start_position(abricate_file, args.goi)
    end=load_end_position(abricate_file,args.goi)
    print(start)
    print(end)
    sys.stdout.flush()
    cc1= ["grep", '>', file] 
    cc2= ["sed", 's/ .*//']
    cc3= ["sed", 's/>//']
    pp1=subprocess.Popen(cc1,stdout=subprocess.PIPE)
    pp2=subprocess.Popen(cc2, stdin=pp1.stdout, stdout=subprocess.PIPE)
    pp3=subprocess.Popen(cc3, stdin=pp2.stdout, stdout=subprocess.PIPE)
    
    pp = pp3.communicate()[0]
    pp=pp.decode()
    #pp=pp.decode()
    out=pd.concat([start,end], axis=1)
    
    out['p']=pp.strip()
    out=out[['p','START','END']]
    
    sys.stdout.flush()
    out_name=str(args.fasta_file1 + "_out")
    out.to_csv(out_name,index=0,sep='\t',header=0)
    subseq_command = ["seqtk", "subseq",file,out_name]
    p=subprocess.Popen(subseq_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    out=out.decode()
    o=open(str(file +"_trim"),'w')
    o.write(out)
    o.close()
    

def main():
    args = get_arguments()
    print(args.fasta_file1)
    run_abricate(args.fasta_file1)
    parse_fasta_file(args.fasta_file1)
    
    start_time=datetime.datetime.now()
    trim=str(args.fasta_file1 + '_trim')    
    mash_command= ['mash','sketch','-s','1000000','-k','21',trim]
    p=subprocess.Popen(mash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = p.communicate()
    out=out.decode()
    
    
    
    
    
    
if __name__ == '__main__':
    main()

