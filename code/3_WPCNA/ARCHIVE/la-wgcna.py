#!/usr/bin/env python3

#------------------------------------------------------------------------------
## Multiresolution partitioning of the protein co-expression network using the 
#  Leiden algorithm.
#------------------------------------------------------------------------------

import os
import sys
from igraph import Graph
from pandas import read_csv, DataFrame

# Specify which network to be analyzed (wt,ko,combined):
data_type = 2
net = ['wtAdjm.csv', 'koAdjm.csv', 'combinedAdjm.csv'][data_type]
print(f'Analyzing the {net.split(".")[0]} network!')

# Read bicor adjacency matrix (no additional soft threshold)..
os.chdir('/mnt/d/projects/Synaptopathy-Proteomics/data')
df = read_csv(net, header = 0, index_col = 0)

# Create edge list.
edges = df.stack().reset_index()
edges.columns = ['protA','protB','weight']

# Define dictionary of nodes.
nodes = dict(zip(df.columns, range(len(df.columns))))

# Create list of edge tuples.
edge_list = list(zip(edges['protA'],edges['protB']))

# Edges need to be referenced by node id (i number). 
el = list(zip([nodes.get(e[0]) for e in edge_list],
    [nodes.get(e[1]) for e in edge_list]))

# Create empty graph.
print("Generating protein co-expression graph!", file = sys.stderr)
g = Graph()

# Add vertices and their labels.
g.add_vertices(len(nodes))
g.vs['label'] = nodes.keys()

# Add edges as list of tuples: (1,2) = node 1 interacts with node 2.
print("Adding edges to graph, this will take several minutes!", file = sys.stderr)
g.add_edges(el)
g.es['weight'] = edges['weight']

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection in the co-expression graph with the Leiden algorithm.
#------------------------------------------------------------------------------

import sys
from leidenalg import Optimiser
from leidenalg import CPMVertexPartition
print('''Generating partition profile for protein co-expression graph!
        This will take several hours!''', file = sys.stderr)
optimiser = Optimiser()
profile = optimiser.resolution_profile(
        g, CPMVertexPartition, weights = 'weight', resolution_range=(0,1))

# Collect key results.
print(f"Examined network at {len(profile)} resolutions!")
results = {
        'Modularity' : [partition.modularity for partition in profile],
        'Membership' : [partition.membership for partition in profile],
        'Summary'    : [partition.summary() for partition in profile],
        'Resolution' : [partition.resolution_parameter for partition in profile]}

# Save membership info.
out = ["wtAdjm_partitions.csv", "koAdjm_partitions.csv", "combinedAdjm_partitions.csv"][data_type]
DataFrame(results['Membership']).to_csv(out[data_type])

# Save other info as csv.
#FIXME: NEED TO remove row index and columns for proper import into self-preservation.r
df = DataFrame.from_dict(results)

df.drop(index='Membership', axis = 1) 

out = ["wtAdjm_partition_profile.csv", "koAdjm_partition_profile.csv", "combinedAdjm_partition_profile.csv"]
df.to_csv(out[data_type])

# ENDOFILE
#------------------------------------------------------------------------------
