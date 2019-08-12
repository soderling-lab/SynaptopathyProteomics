#!/usr/bin/env python3

#------------------------------------------------------------------------------
## Analysis of the Synaptosome PPI graph with the Leiden algorithm.
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

# Write graph to file.
g.write_pajek("Synaptome.net")
import sys
sys.exit()
# Find partition that optimizes modularity.
partition = la.find_partition(g, la.ModularityVertexPartition, n_iterations = -1)
partition.modularity
partition.summary()

# Find the partition that optimizes CPM.
partitionCPM = la.find_partition(g, la.CPMVertexPartition, n_iterations = -1)
partitionCPM.modularity
partitionCPM.summary() # Why do the two answers disagree so much?

# Examine partition:
print(partition.summary())
print(f"Modularity: {partition.modularity:.4f}")

# Fine tuning of the partition.
# If n_iter < 0, then optimiser continues until no improvement.
optimiser = la.Optimiser()
diff = optimiser.optimise_partition(partition, n_iterations=-1) # 0.0
diff

#FIXME: Why doesn't modularity optimization work?
# Optimize modularity: Examine clustering at multiple resolutions:
# This will take some additional time...
print("Generating partition profile of Synaptic proteome...", file=sys.stderr)
optimiser = la.Optimiser()
profile = optimiser.resolution_profile(g, la.ModularityVertexPartition, resolution_range=(0,1))

# Optimize CPM: Examine clustering at multiple resolutions:
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

# Examine partition:
df = DataFrame.from_dict(results)
bestQ = df.loc[df['Modularity'].idxmax(), 'Modularity']
bestS = df.loc[df['Modularity'].idxmax(), 'Summary']
bestR = df.loc[df['Modularity'].idxmax(), 'Resolution']
print(f"Best Partition : {bestS}")
print(f"Best Modularity: {bestQ}")
print(f"Best Resolution: {bestR}")

# Convert dict to pandas df.
df = DataFrame(results['Membership'])
df.to_csv("ppi_graph_partitions.csv")

# Save as csv.
df = DataFrame.from_dict(results)
df.to_csv("ppi_graph_partition_profile.csv")

# ENDOFILE
#------------------------------------------------------------------------------
