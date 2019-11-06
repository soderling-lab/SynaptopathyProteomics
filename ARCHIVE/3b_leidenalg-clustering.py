#!/usr/bin/env python3
''' Multiresolution clustering of the protein co-expression graph. 
Executable script that utilizes the Leiden algorithm to perform multiresolution 
clustering of the protein co-expression graph.

FIXME: Doesnt work for single resolution parameter...
FIXME: Remove row index and columns for proper function of self-preservation.r

USAGE: 
    $ ./leidenalg-clustering.py adjm.csv --resolution

OUTPUT:
    la-partitions.csv: node community (cluster) membership for graph partition(s).
    la-profile.csv: other descriptive statistics for the resolution profile.
'''

#------------------------------------------------------------------------------
## Parse the user's input.
#------------------------------------------------------------------------------

import os
import sys
from argparse import ArgumentParser
from pandas import read_csv

ap = ArgumentParser(description = '''Perform clustering of the protein 
        co-expression matrix using the Leiden algorithm. ''')
ap.add_argument('adjm.csv', type = str,
        help = 'File path to a n x n adjaceny matrix.')
ap.add_argument('--beta','-b' ,type = int, default = 1,
        help = '''Power (beta) for soft-thresholding of the adjacency matrix
        Defaults to 1, (no soft-threshold).''')
ap.add_argument('--resolution', '-r', type = float, nargs='+', default = (0,1),
        help = '''The resolution range for which to analyze the network. 
        Defaults to (0,1). A single resolution can be be specified by passing a 
        single value, e.g. 0.05.''')
args = vars(ap.parse_args())

# Read bicor adjacency matrix.
adjm = read_csv(args['adjm.csv'], header = 0, index_col = 0)

#------------------------------------------------------------------------------
## Create an igraph object to be passed to leidenalg.
#------------------------------------------------------------------------------

from igraph import Graph

# Create edge list.
edges = adjm.stack().reset_index()
edges.columns = ['protA','protB','weight']

# Soft-thresholding with power, beta.
power = args['beta'] 
edges['weight'] = edges['weight'] ** power

# Define dictionary of nodes.
nodes = dict(zip(adjm.columns, range(len(adjm.columns))))

# Create list of edge tuples.
edge_list = list(zip(edges['protA'],edges['protB']))

# Edges need to be referenced by node id (a number). 
el = list(zip([nodes.get(e[0]) for e in edge_list],
    [nodes.get(e[1]) for e in edge_list]))

# Create empty graph.
print("Generating igraph object:", file = sys.stderr)
g = Graph()

# Add vertices and their labels.
g.add_vertices(len(nodes))
g.vs['label'] = nodes.keys()

# Add edges as list of tuples; ex: (1,2) = node 1 interacts with node 2.
print("Adding edges to graph, this will take several minutes...", file = sys.stderr)
g.add_edges(el)
g.es['weight'] = edges['weight'] 

# Remove any self-loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection in the WPCNetwork with the Leiden algorithm.
#------------------------------------------------------------------------------

import sys
from pandas import DataFrame
from leidenalg import Optimiser
from leidenalg import find_partition
from leidenalg import CPMVertexPartition

# Perform community detection with leidenalg.
print('''Generating partition profile for protein co-expression graph!
        This will take several hours!''', file = sys.stderr)
optimiser = Optimiser()
profile = optimiser.resolution_profile(
        g, 
        CPMVertexPartition, 
        weights = 'weight', 
        resolution_range = args['resolution'])

# Collect partition results. 
print(f"Examined network at {len(profile)} resolutions!", file = sys.stderr)
results = {
        'Modularity' : [partition.modularity for partition in profile],
        'Membership' : [partition.membership for partition in profile],
        'Summary'    : [partition.summary() for partition in profile],
        'Resolution' : [partition.resolution_parameter for partition in profile]}

# Save membership (partition) info.
DataFrame(results['Membership']).to_csv("la-partitions.csv")

# Save other info as csv.
df = DataFrame.from_dict(results)
df.to_csv("la-profile.csv")

# ENDOFILE
#------------------------------------------------------------------------------
