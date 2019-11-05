#!/usr/bin/env python3

''' Multiresolution clustering of the protein co-expression graph. 
Executable script that utilizes the Leiden algorithm to perform multiresolution 
clustering of the protein co-expression graph.

USAGE: 
    $ ./leidenalg-clustering.py

OUTPUT:
    la-partitions.csv: node community (cluster) membership for graph partition(s).
    la-profile.csv: other descriptive statistics for the resolution profile.
'''

#------------------------------------------------------------------------------
## Load the adjacency matrix.
#------------------------------------------------------------------------------

from sys import stderr
from pandas import read_csv

# Which analysis are we doing?
# WT or KO network (0,1)
data_type = 0 
geno = ['WT','KO'][data_type]
print("Performing Leiden algorithm clustering of the " + geno + " protein co-expression network.", file = stderr)

# Read bicor adjacency matrix.
datadir = '~/projects/Synaptopathy-Proteoimcs/data'
myfiles = ['~/projects/Synaptopathy-Proteomics/data/3_WTadjm.csv',
        '~/projects/Synaptopathy-Proteomics/data/3_KOadjm.csv']
adjm = read_csv(myfiles[data_type], header = 0, index_col = 0)

#------------------------------------------------------------------------------
## Create an igraph object to be passed to leidenalg.
#------------------------------------------------------------------------------
# I tried a couple ways of creating an igraph object. Simplier approaches like
# using the igraph.Weighted_Adjacency function didn't work. 

from igraph import Graph

# Create edge list.
edges = adjm.stack().reset_index()
edges.columns = ['protA','protB','weight']

# Define dictionary of nodes.
nodes = dict(zip(adjm.columns, range(len(adjm.columns))))

# Create list of edge tuples.
edge_list = list(zip(edges['protA'],edges['protB']))

# Edges need to be referenced by node id (a number). 
el = list(zip([nodes.get(e[0]) for e in edge_list],
    [nodes.get(e[1]) for e in edge_list]))

# Create empty graph.
print("Generating igraph object:", file = stderr)
g = Graph()

# Add vertices and their labels.
g.add_vertices(len(nodes))
g.vs['label'] = nodes.keys()

# Add edges as list of tuples; ex: (1,2) = node 1 interacts with node 2.
print("Adding edges to graph, this will take several minutes...", file = stderr)
g.add_edges(el)
g.es['weight'] = edges['weight'] 

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection in the PCNetwork with the Leiden algorithm.
#------------------------------------------------------------------------------

from numpy import linspace
from leidenalg import find_partition
from leidenalg import CPMVertexPartition
from progressbar import ProgressBar

# Loop to perform leidenalg community detection at 100 resolutions.
print('''Generating partition profile for protein co-expression graph!
        This will take several hours!''', file = stderr)
pbar = ProgressBar()
resolution_range = linspace(0,1,100)
profile = list()
for res in pbar(resolution_range):
    partition = find_partition(g, CPMVertexPartition, 
            weights='weight', resolution_parameter=res)
    profile.append(partition)

#------------------------------------------------------------------------------
## Save partition profile results.
#------------------------------------------------------------------------------

from pandas import DataFrame

# Collect partition results. 
print(f"Examined network at {len(profile)} resolutions!", file = stderr)
results = {
        'Modularity' : [partition.modularity for partition in profile],
        'Membership' : [partition.membership for partition in profile],
        'Summary'    : [partition.summary() for partition in profile],
        'Resolution' : [partition.resolution_parameter for partition in profile]}

# Save cluster membership.
# FIXME: add correct names!
myfile = datadir + geno + "_partitions.csv"
DataFrame(results['Membership']).to_csv(myfile)

# Save partition profile.
# Add correct names!
df = DataFrame.from_dict(results)
myfile = datadir + geno + "_profile.csv"
df.to_csv(myfile)

# ENDOFILE
#------------------------------------------------------------------------------
