import os
import glob
import tempfile
import subprocess
import collections

import pandas as pd
import networkx as nx

"""
functions to cluster output sequences from flanker
build mash sketch and pairwise_mash_distances are functions written by Ryan Wick in his assembly de-replicator repo
https://github.com/rrwick/Assembly-Dereplicator. The other functions are also adapted from functions written by Ryan
for assembly de-replicator.
"""

def find_all_assemblies():
    all_assemblies=[]
    for foldername, subfolders, filenames in os.walk(os.getcwd()):
        for filename in filenames:

            #might need something more sophisticated than 50

            if os.path.getsize(os.path.join(foldername,filename)) > 100:
                if filename.endswith('flank.fasta'):
                    all_assemblies.append(filename)


    print('found {:,} files\n'.format(len(all_assemblies)))

    return all_assemblies


def build_mash_sketch(assemblies, threads, temp_dir, sketch_size):
    mash_command = ['mash', 'sketch', '-p', str(threads), '-o', temp_dir + '/mash',
                        '-s', str(sketch_size)] + assemblies

    
    subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return temp_dir + '/mash.msh'

def pairwise_mash_distances(mash_sketch, threads):
    mash_command = ['mash', 'dist', '-p', str(threads), mash_sketch, mash_sketch]

    mash_out = subprocess.run(mash_command, stdout=subprocess.PIPE).stdout.decode()
    return mash_out.splitlines()



def create_graph_from_distances(pairwise_distances, threshold):

    matrix=[]
    assemblies = set()
    graph = collections.defaultdict(set)
    all_connections = collections.defaultdict(set)
    for line in pairwise_distances:
        parts = line.split('\t')
        assembly_1 = parts[0]
        assembly_2 = parts[1]
        distance = float(parts[2])

        matrix.append([assembly_1,assembly_2,distance])



    df=pd.DataFrame(matrix,columns=['assembly_1','assembly_2','distance'])
    df=df[df.distance <= threshold]
    G=nx.from_pandas_edgelist(df,'assembly_1','assembly_2','distance')
    l=list(nx.connected_components(G))
    L=[dict.fromkeys(y,x) for x, y in enumerate(l)]

    d={k: v for d in L for k, v in d.items()}
    df2=df['assembly_1'].unique()
    df2=pd.DataFrame(df2,columns=['assembly_1'])

    df2['cluster']=df2.assembly_1.map(d)

    return(df2)

def flank_scrub():

    filelist=glob.glob(str(str(os.getcwd()) + '/' + str("*flank.fasta")))

    for filename in filelist:
        os.remove(filename)


#here we build clusters using mash distances
def define_clusters(gene,window,threads,threshold,outfile):


    with tempfile.TemporaryDirectory() as temp_dir:
        all_assemblies=find_all_assemblies()
        mash_sketch = build_mash_sketch(all_assemblies, threads, temp_dir, 1000)

        pairwise_distances= pairwise_mash_distances(mash_sketch,threads)



        clusters=create_graph_from_distances(pairwise_distances, float(threshold))

        clusters.to_csv(str(outfile + '_' + gene.strip() + '_' + str(window)),index=False)
