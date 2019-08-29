#!/usr/bin/env python3

#------------------------------------------------------------------------------
## Multiresolution partitioning of the protein co-expression network using the 
#  Leiden algorithm.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
## Parse the user's input.
#------------------------------------------------------------------------------

import os
import sys
from argparse import ArgumentParser
from pandas import read_csv

ap = ArgumentParser(description = ''' Perform clustering of the protein 
        co-expression matrix using the Leiden algorithm. ''')
ap.add_argument('adjm.csv', type = str,
        help = 'File path to a n x n adjaceny matrix.')
args = vars(ap.parse_args())

# Read bicor adjacency matrix.
adjm = read_csv(args['adjm.csv'], header = 0, index_col = 0)

#------------------------------------------------------------------------------
## Create an igraph object to be passed to leidenalg.
#------------------------------------------------------------------------------

from igraph import Graph

# Create edge list.
edges = df.stack().reset_index()
edges.columns = ['protA','protB','weight']

# Add weighted edges.
power = 1 
edges['weight'] = edges['weight'] ** power

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
g.es['weighted'] = edges['weighted'] 

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection in the WPCNetwork with the Leiden algorithm..
#------------------------------------------------------------------------------

import sys
from pandas import DataFrame
from leidenalg import Optimiser
from leidenalg import find_partition
from leidenalg import CPMVertexPartition

# Calculate resolution profile:
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
DataFrame(results['Membership']).to_csv("la-partitions.csv")

# Save other info as csv.
#FIXME: NEED TO remove row index and columns for proper import into self-preservation.r script.
df = DataFrame.from_dict(results)
df.to_csv("la-profile.csv")

# ENDOFILE
#------------------------------------------------------------------------------
