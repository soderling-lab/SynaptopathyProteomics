#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

# FIXME: Something isn't working... clustering of ppi graph returns one
# module...

## User parameters: 
input_adjm = "3_PPI_Adjm.csv" # Input adjacency matrix.
#input_adjm = "3_GO_Semantic_Similarity_RMS_Adjm.csv"
output_name = "PPI" # Output filename.
#method = 'SurpriseVertexPartition' # Best for PPI graph.
method = "ModularityVertexPartition" 
#method = 'CPMVertexPartition' # For signed co-expression graph and GO graph.
sft = 1 # Power (soft-threshold) for weighting the network. See notes. 
weighted = False

## Resolution range for multiresolution methods:
rmin = 0 # Min resolution.
rmax = 1 # Max resolution. 
nsteps = 100 # Number of steps.

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
## Prepare the workspace.
#------------------------------------------------------------------------------

# Imports.
import sys
import os
import glob
from os.path import dirname
from sys import stderr

# Directories.
here = os.getcwd()
root = dirname(dirname(here))
datadir = os.path.join(root,"rdata")
funcdir = os.path.join(root,"Py")

# Load functions.
sys.path.append(root)
from Py import myfun

# Get system variables.
myvars = ['SLURM_JOBID','SLURM_CPUS_PER_TASK']
envars = {var:os.environ.get(var) for var in myvars}
jobID = myfun.xstr(envars['SLURM_JOBID'])

#------------------------------------------------------------------------------
## Load input adjacency matrix and create an igraph object.
#------------------------------------------------------------------------------
# I tried a couple ways of creating an igraph object. Simplier approaches like
# using the igraph.Weighted_Adjacency function didn't work for me...

from pandas import read_csv
from igraph import Graph

# Read bicor adjacency matrix as input.
myfile = os.path.join(datadir,input_adjm)
adjm = read_csv(myfile, header = 0, index_col = 0)

# Add rownames.
adjm = adjm.set_index(keys=adjm.columns)

# Create edge list.
edges = adjm.stack().reset_index()
edges.columns = ['protA','protB','weight']

# Remove weight == 0.
edges[edges.weight != 0]

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
if weighted:
    g.es['weight'] = edges['weight']**sft
elif not g.is_weighted():
    print("Analyzing an unweighted graph.")

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

# Get single largest connected component.
# Errors may be caused by la clustering of un-connected components.
#g = g.clusters().giant()

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
single_resolution = myfun.contains(out,method)

# Perform Leidenalg community detection. 
if (single_resolution):
    # Single resolution clustering methods:
    profile = list()
    if weighted:  
        # Analysis of weighted graph at single resolution.
        partition = find_partition(g, partition_type, weights='weight')
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        partition = myfun.filter_modules(partition)
        profile.append(partition)
    else:
        # Analysis of unweighted graph at single resolution.
        partition = find_partition(g, partition_type)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        partition = myfun.filter_modules(partition)
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
        if weighted:
            # Analysis of weighted graph.
            partition = find_partition(g, partition_type, 
                    weights='weight', resolution_parameter=resolution)
            optimiser = Optimiser()
            diff = optimiser.optimise_partition(partition,n_iterations=-1)
            partition = myfun.filter_modules(partition)
            profile.append(partition)
        else:
            # Analysis of unweighted graph.
            partition = find_partition(g, partition_type, 
                    resolution_parameter=resolution)
            optimiser = Optimiser()
            diff = optimiser.optimise_partition(partition,n_iterations=-1)
            partition = myfun.filter_modules(partition)
            profile.append(partition)
        # Ends loop.
# Ends If/else.

print(f"Complete! Examined network at {len(profile)} resolutions!", 
        file = stderr)

#------------------------------------------------------------------------------
## Save clustering results.
#------------------------------------------------------------------------------

from pandas import DataFrame

# Collect partition results and save as csv. 
if len(profile) is 1:
    # Single resolution profile:
    results = {
            'Modularity' : [partition.modularity for partition in profile],
            'Membership' : [partition.membership for partition in profile],
            'Summary'    : [partition.summary() for partition in profile]}
else: 
    # Multi-resolution profile:
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




