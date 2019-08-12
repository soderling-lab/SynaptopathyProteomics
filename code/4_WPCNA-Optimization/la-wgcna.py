#!/usr/bin/env python3

#------------------------------------------------------------------------------
## Multiresolution partitioning of the WPCNetwork using the Leiden algorithm.

#------------------------------------------------------------------------------
# ## Analysis of the Weighted Protein Co-expression Network.
#------------------------------------------------------------------------------

import os
import sys
from igraph import Graph
from pandas import read_csv, DataFrame

# Read bicor adjacency matrix (no additional soft threshold)..
os.chdir('/mnt/d/projects/Synaptopathy-Proteomics/code/4_WPCNA-Optimization')
df = read_csv('wtAdjm.csv', header = 0, index_col = 0)

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

# Power for weighting the network.
beta = 13  # scale free fit ~ 13

# Add edges as list of tuples: (1,2) = node 1 interacts with node 2.
print("Adding edges to graph, this will take several minutes!", file = sys.stderr)
g.add_edges(el)
g.es['weight'] = edges['weight']**beta

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection in the WPCNetwork with the Leiden algorithm..
#------------------------------------------------------------------------------

import sys

from ttictoc import TicToc

from leidenalg import Optimiser
from leidenalg import CPMVertexPartition

# Examine resolution profile:
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

# Save partition info.
DataFrame(results['Membership']).to_csv("wtAdjm_weighted_partitions.csv")

# Save other info as csv.
df = DataFrame.from_dict(results)
df.to_csv("wtAdjm_weighted_partition_profile.csv")

# ENDOFILE
#------------------------------------------------------------------------------
