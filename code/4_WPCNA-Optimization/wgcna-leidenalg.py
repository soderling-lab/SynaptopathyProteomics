#!/usr/bin/env python3

#------------------------------------------------------------------------------
## Analysis of the WPCNetwork with the Leiden algorithm.
#  Use the CPM to evaluate quality of the network at multiple resolutions.


#------------------------------------------------------------------------------
# ## Analysis of the Synaptic PPI graph.
#------------------------------------------------------------------------------

import os
import sys

import numpy as np
import leidenalg as la

from pandas import read_csv, DataFrame
from igraph import Graph

# Load Compiled PPIs
os.chdir('/mnt/d/projects/Synaptopathy-Proteomics/tables')
df = read_csv('3_Compiled_PPIs.csv')

# Define dictionary of nodes.
v = np.unique([df['musEntrezA'], df['musEntrezB']])
nodes = dict(zip(v, range(len(v))))

# Create list of edges.
edges = list(zip(df['musEntrezA'], df['musEntrezB']))

# Edges need to be referenced by node id (i number). 
el = list(zip([nodes.get(e[0]) for e in edges],
    [nodes.get(e[1]) for e in edges]))

# Create empty graph.
g = Graph()

# Add vertices and their labels..
g.add_vertices(len(v))
g.vs['label'] = v

# Add edges as list of tuples: (1,2) ~node 1 interacts with node 2.
g.add_edges(el)

# Simplify to remove loops and multiple edges.
g = g.simplify()

# Find partition that optimizes modularity.
partition = la.find_partition(g, la.ModularityVertexPartition, n_iterations = -1)

# Examine partition:
print(partition.summary())
print(f"Modularity: {partition.modularity:.4f}")

# Fine tuning of the partition.
# If n_iter < 0, then optimiser continues until no improvement.
optimiser = la.Optimiser()
diff = optimiser.optimise_partition(partition, n_iterations=-1) # 0.0

# Try again with CPM algorithm at a given resolution.
partition = la.find_partition(g, la.CPMVertexPartition, resolution_parameter = 0.05)
print(partition.summary())
print(f"Resolution: {partition.resolution_parameter}")
print(f"Modularity: {partition.modularity:.4f}")

# Examine clustering at multiple resolutions:
# This will take some additional time...
print("Generating partition profile of Synaptic proteome...", file=sys.stderr)
optimiser = la.Optimiser()
profile = optimiser.resolution_profile(g, la.CPMVertexPartition, resolution_range=(0,1))

# Collect key results.
print(f"Network profile was examined at {len(profile)} different resolutions!") 
results = {
        'Modularity' : [partition.modularity for partition in profile],
        'Membership' : [partition.membership for partition in profile],
        'Summary'    : [partition.summary() for partition in profile],
        'Resolution' : [partition.resolution_parameter for partition in profile]}

# Save as csv.
df = DataFrame.from_dict(results)
df.to_csv("ppi_graph_partition_profile.csv")

#------------------------------------------------------------------------------
# ## Analysis of the Weighted Protein Co-expression Network.
#------------------------------------------------------------------------------

import os
import sys
from pandas import read_csv, DataFrame

# Read bicor adjacency matrix (not weighted).
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

# Add edges as list of tuples: (1,2) = node 1 interacts with node 2.
print("Adding edges to graph, this will take several minutes!", file = sys.stderr)
g.add_edges(el)
beta = 1.0 # Power for weighting the network.
g.es['weight'] = edges['weight']**beta

# Remove self loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection in the WPCNetwork with the Leiden algorithm..
#------------------------------------------------------------------------------

import sys

from ttictoc import TicToc

from leidenalg import Optimiser
from leidenalg import CPMVertexPartition

# Examine network at single resolution.
t = TicToc()
t.tic()
partition = la.find_partition(g, la.CPMVertexPartition, resolution_parameter = 0.05)
t.toc()
print(f"Time elapsed analysing network at single resolution: {t.elapsed}")

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

# Save as csv.
df = DataFrame.from_dict(results)
df.to_csv("wtAdjm_partition_profile.csv")

# Done!
sys.exit()

# ENDOFILE
#------------------------------------------------------------------------------



df = read_csv('wtDat.csv', header = 0, index_col = 0)

import numpy as np
dm = np.transpse(df.values)
cormat = np.corrcoef(dm)
cormat.shape

from astropy.stats.biweight import biweight_midcorrelation as bicor

x = dm[0]
y = dm[1]
bicor(x, y) # Value is the same as from WGCNA::bicor() in R!

# how to compute for all pairs?

# Remove diagonal. 
df.values[[np.arange(len(df))]*2] = np.nan

