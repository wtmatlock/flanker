#from flanker import get_arguments
#from flanker import run_abricate
import tempfile
import os
import subprocess
import collections
import pandas as pd
import glob

"""
functions to cluster output sequences from flanker
most are written/adapted from function written by Ryan Wick in his assembly de-replicator repo
https://github.com/rrwick/Assembly-Dereplicator
"""

def find_all_assemblies(in_dir):
    all_assemblies=[]
    for foldername, subfolders, filenames in os.walk(in_dir):
        for filename in filenames:
            #might need something more sophisticated than 50
            if os.path.getsize(os.path.join(foldername,filename)) > 100:
                if filename.endswith('flank.fasta'):
                    #print(filename)
                    #print(os.path.getsize(os.path.join(foldername,filename)))
                    all_assemblies.append(filename)
                    #print(all_assemblies)

    print('found {:,} files\n'.format(len(all_assemblies)))
    return all_assemblies


def build_mash_sketch(assemblies, threads, temp_dir, sketch_size):
    mash_command = ['mash', 'sketch', '-p', str(threads), '-o', temp_dir + '/mash',
                        '-s', str(sketch_size)] + assemblies
    print(mash_command)
    subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return temp_dir + '/mash.msh'

def pairwise_mash_distances(mash_sketch, threads):
    mash_command = ['mash', 'dist', '-p', str(threads), mash_sketch, mash_sketch]
    print(mash_command)
    mash_out = subprocess.run(mash_command, stdout=subprocess.PIPE).stdout.decode()
    return mash_out.splitlines()

def dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited

def create_graph_from_distances(pairwise_distances, threshold):
    """
    Builds an undirected graph where nodes are assemblies and edges connect assemblies which have
    a pairwise Mash distance below the threshold.
    """
    assemblies = set()
    graph = collections.defaultdict(set)
    all_connections = collections.defaultdict(set)
    for line in pairwise_distances:
        parts = line.split('\t')
        assembly_1 = parts[0]
        assembly_2 = parts[1]
        distance = float(parts[2])
        assemblies.add(assembly_1)
        assemblies.add(assembly_2)
        if assembly_1 == assembly_2:
            continue
        all_connections[assembly_1].add(assembly_2)
        all_connections[assembly_2].add(assembly_1)

        if distance <= threshold:
            graph[assembly_1].add(assembly_2)
            graph[assembly_2].add(assembly_1)
    assemblies = sorted(assemblies)
    assembly_count = len(assemblies)
    for assembly in assemblies:  # sanity check: make sure we have all the connections
        assert len(all_connections[assembly]) == assembly_count - 1
    return assemblies, graph


def cluster_assemblies(assemblies, graph):
    visited = set()
    clusters = []
    for assembly in assemblies:
        if assembly in visited:
            continue
        connected = dfs(graph, assembly)
        clusters.append(sorted(connected))
        visited |= connected
    return clusters

def cleanup():
    clean = ['rm','*flank.fasta', '*.msh']
    clean_out = subprocess.run(clean, stdout=subprocess.PIPE).stdout.decode()


#here we build clusters using mash distances
def define_clusters(gene,window,in_dir,threads,threshold,outfile):
    #args=get_arguments()
    print(threshold)

    with tempfile.TemporaryDirectory() as temp_dir:
        all_assemblies=find_all_assemblies(in_dir)
        mash_sketch = build_mash_sketch(all_assemblies, threads, temp_dir, 1000)

        pairwise_distances= pairwise_mash_distances(mash_sketch,threads)
        assemblies, graph = create_graph_from_distances(pairwise_distances, float(threshold))
        clusters = cluster_assemblies(assemblies, graph)

        n_clusters= len(clusters)
        print("\nFound " + str(n_clusters) + " clusters at a " + str(threshold) +" threshold.\n")
        i=1
        all=[]
        for cluster in clusters:
            for c in cluster:
                m=[c,i]

                all.append(m)

            i+=1

        df=pd.DataFrame(all,columns=['Isolate','Cluster'])
        print(df)

        df.to_csv(str(outfile + '_' + gene.strip() + '_' + str(window)),index=False)
