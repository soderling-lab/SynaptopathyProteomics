#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

# FIXME: error when saving data. single resolution methods don't have resoltuion parameter...

## User parameters: 
input_adjm = "3_PPI_Adjm.csv" # Input adjacency matrix.
#input_adjm = "3_GO_Semantic_Similarity_RMS_Adjm.csv"
output_name = "PPI" # Output filename.
method = 'SurpriseVertexPartition' # Best for PPI graph.
#method = 'CPMVertexPartition' # For signed co-expression graph and GO graph.
sft = 1 # Power (soft-threshold) for weighting the network. See notes. 

## For multiresolution methods:
rmin = 1 # Min resolution.
rmax = 1 # Max resolution. 
nsteps = 1 # Number of steps.

## Notes: sft -- This is the power to which the adjacency matrix is raised in
# order to apply a soft-threshold to the network. If this power is even then the
# network will become unsigned. 

# For PPI graph use SupriseVertexPartition -- this minimizes the percent
# unclustered while maximizing the number of modules.
# For GO graph use CPM with resolution parameter = 1.0. This minimizes the
# percent unclustered while retaining a large number of modules.

## Methods: Leidenalg supports the following methods for optimization
#           of community detection:
#
# Methods for positive edge weights, no resolution parameter:
# [1] ModularityVertexPartition -- Implements modularity. Only 
#     well-defined for positive edge weights. No resolution parameter.
# [2] SurpriseVertexPartition -- Implements (asymptotic) Surprise. 
#     This quality function is well-defined only for positive edge 
#     weights. No resolution parameter.
#
# Methods for positive edges with linear resolution parameter:
# [3] RBConfigurationVertexPartition -- Implements Reichardt and 
#     Bornholdt’s Potts model with a configuration null model. 
#     Only well-defined for positive edge weights.
# [4] RBERVertexPartition -- Implements Reichardt and Bornholdt’s 
#     Potts model with a configuration null model. 
#     Only well-defined for positive edge weights. 
#
# Methods for positive and negative edge weights:
# [5] CPMVertexPartition -- Implements CPM. Quality function is 
#     well-defined for both positive and negative edge weights. 
#     Can utilize a linear resolution parameter.
#
# Method for unweighted graphs without resolution parameter:
# [6] SignificanceVertexPartition -- Implements Significance. This 
#     quality function is well-defined only for unweighted graphs.
#     Does not utilize a resolution parameter, but tries to find
#     the best resolution/partition.

#------------------------------------------------------------------------------
## Define some utility functions.
#------------------------------------------------------------------------------

# xstr
def xstr(s):
    ''' Convert NoneType to blank ('') string.'''
    if s is None:
        return ''
    else:
        return str(s)
# EOF

# contains
def contains(mylist,value,return_index=False):
    ''' Check if list contains a value. 
    Like list.index(value) but returns False if the provided list
    does not contain value. 
    '''
    list_as_dict = dict(zip(mylist,range(len(mylist))))
    index = list_as_dict.get(value)
    if not return_index: 
        return type(index) is int # True if in list.
    else:
        return index
# EOF

#------------------------------------------------------------------------------
## Load the input adjacency matrix.
#------------------------------------------------------------------------------

# Imports.
import os
import glob
from os.path import dirname
from sys import stderr
from pandas import read_csv

# Get system variables.
myvars = ['SLURM_JOBID','SLURM_CPUS_PER_TASK']
envars = {var:os.environ.get(var) for var in myvars}
jobID = xstr(envars['SLURM_JOBID'])

# Read bicor adjacency matrix as input.
here = os.getcwd()
root = dirname(dirname(here))
datadir = os.path.join(root,"rdata")
myfile = os.path.join(datadir,input_adjm)
adjm = read_csv(myfile, header = 0, index_col = 0)

# Add rownames.
adjm = adjm.set_index(keys=adjm.columns)

#------------------------------------------------------------------------------
## Create an igraph object.
#------------------------------------------------------------------------------
# I tried a couple ways of creating an igraph object. Simplier approaches like
# using the igraph.Weighted_Adjacency function didn't work for me...

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
g = Graph()

# Add vertices and their labels.
g.add_vertices(len(nodes))
g.vs['label'] = nodes.keys()

# Add edges as list of tuples; ex: (1,2) = node 1 interacts with node 2.
# This will take several minutes.
g.add_edges(el)
g.es['weight'] = edges['weight']**sft

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection with the Leiden algorithm.
#------------------------------------------------------------------------------

from numpy import linspace
from leidenalg import Optimiser
from leidenalg import find_partition
from progressbar import ProgressBar
from importlib import import_module

# Dynamically load the partition_type class--the clusering optimization method.
partition_type = getattr(import_module('leidenalg'), method)

# Methods that don't support multi-resolution clustering:
out = ["ModularityVertexPartition", 
        "SurpriseVertexPartition",
        "SignificanceVertexPartition"]
# Check if users optimization method supports resolution parameter. 
single_resolution = contains(out,method)

# Perform Leidenalg community detection. 
if (single_resolution):
    # Single resolution clustering:
    profile = list()
    if method is "SignificanceVertexPartition":  
        # SignificanceVertexPartition only supports unweighted graphs.
        partition = find_partition(g, partition_type, weights=None)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        profile.append(partition)
    else:
        # Analysis of weighted graph at single resolution.
        partition = find_partition(g, partition_type, weights='weight')
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        profile.append(partition)
else:
    # Loop to perform multi-resolution clustering:
    print("Performing Leiden algorithm clustering of the" 
        " protein co-expression network.\n", file = stderr)
    pbar = ProgressBar()
    resolution_range = linspace(rmin,rmax,nsteps)
    profile = list()
    for resolution in pbar(resolution_range):
        # Perfrom La clustering.
        partition = find_partition(g, partition_type, 
                weights='weight', resolution_parameter=resolution)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        profile.append(partition)
        # Ends loop.
# Ends If/else.

#------------------------------------------------------------------------------
## Save clustering results.
#------------------------------------------------------------------------------

from pandas import DataFrame

# Collect partition results. 
print(f"Complete! Examined network at {len(profile)} resolutions!", 
        file = stderr)
results = {
        'Modularity' : [partition.modularity for partition in profile],
        'Membership' : [partition.membership for partition in profile],
        'Summary'    : [partition.summary() for partition in profile],
        'Resolution' : [partition.resolution_parameter for partition in profile]}

# Save cluster membership vectors.
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "_partitions.csv")
DataFrame(results['Membership']).to_csv(myfile)

# Save partition profile summary data.
df = DataFrame.from_dict(results)
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "_profile.csv")
df.to_csv(myfile)
