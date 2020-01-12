#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

## User parameters: 
adjm_type = "Cortex" # See adjms below.
method = 4 # See methods below.
rmin = 0
rmax = 1
nsteps = 100

#------------------------------------------------------------------------------
## Parse the user provided parameters.
#------------------------------------------------------------------------------

import sys
from sys import stderr

## Input adjacency matrix.
adjms = {"Cortex" : "3_Cortex_Adjm.csv",
        "Striatum" : "3_Striatum_Adjm.csv",
        "Combined" : "3_Combined_Adjm.csv",
        "PPI" : "3_PPI_Adjm.csv",
        "GO" : "3_GO_Semantic_Similarity_RMS_Adjm.csv"}

## Leidenalg supports the following optimization methods:
methods = {
        # Modularity
        0: {'partition_type' : 'ModularityVertexPartition', 
            'weights' : 'positive',
            'resolution_parameter' : None},
        # Surprise
        1: {'partition_type' : 'SurpriseVertexPartition', 
            'weights' : 'positive',
            'resolution_parameter' : None},
        # RBConfig
        2: {'partition_type' : 'RBConfigurationVertexPartition', 
            'weights' : 'positive',
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps}},
        # RBEVertex
        3: {'partition_type' : 'RBERVertexPartition', 
            'weights' : 'positive',
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps}},
        # CPM
        4: {'partition_type' : 'CPMVertexPartition', 
            'weights' : 'positive and negative',
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps}},
        # Significance
        5: {'partition_type' : 'SignificanceVertexPartition', 
            'weights':None,
            'resolution_parameter' : None}}

# Parameters for clustering.
parameters = methods.get(method)
method = parameters.get('partition_type')

# Status report.
print("Performing Leidenalg clustering of the {}".format(adjm_type),
        "graph utilizing the {}".format(method),
        "method to find optimal partition(s)...", file=stderr)

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

import os
from os.path import dirname
import glob

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

from pandas import read_csv, DataFrame
from igraph import Graph

# Read bicor adjacency matrix.
input_adjm = adjms.get(adjm_type)
myfile = os.path.join(datadir,input_adjm)
adjm = read_csv(myfile, header = 0, index_col = 0)
adjm = adjm.set_index(keys=adjm.columns) # Add row names.
adjm = adjm.fillna(0)

# Incorporate weighted, positive, unweighted graph types.
if parameters['weights'] is 'positive':
    A = abs(adjm.values)
else:
    A = adjm.values

# Create igraph object.
g = Graph.Adjacency(A.tolist())
#g = Graph.Adjacency((A > 0).tolist())
g.es['weight'] = A[A.nonzero()]
g.vs['label'] = adjm.columns

# Update weights parameter for weighted graphs.
if parameters.get('weights') is not None:
    parameters['weights'] = 'weight'

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

# Add graph to input parameters for clustering.
parameters['graph'] = g

#------------------------------------------------------------------------------
## Community detection with the Leiden algorithm.
#------------------------------------------------------------------------------

import numpy as np
from numpy import linspace
from leidenalg import Optimiser, find_partition
from progressbar import ProgressBar
from importlib import import_module

# Update partition type parameter.
# Dynamically load the partition_type class. This is the method to be used for
# clusering optimization.
parameters['partition_type'] = getattr(import_module('leidenalg'),method)

# Update n_iterations parameter.
parameters['n_iterations'] = -1

# Remove any None type parameters.
out = [key for key in parameters if parameters.get(key) is None]
for key in out: del parameters[key]

# Perform Leidenalg community detection. 
if parameters.get('resolution_parameter') is None:
    # Single resolution methods.
    profile = list()
    partition = find_partition(**parameters)
    optimiser = Optimiser()
    diff = optimiser.optimise_partition(partition,n_iterations=-1)
    partition = myfun.filter_modules(partition)
    profile.append(partition)
    m = np.array(partition.membership)
    unclustered = sum(m==0)/len(m)
    print(partition.summary())
    print("Percent unclustered: {}".format(unclustered) + " (%).\n")
else:
    # Loop to perform multi-resolution clustering methods.
    pbar = ProgressBar()
    profile = list()
    resolution_range = linspace(**parameters.get('resolution_parameter'))
    for resolution in pbar(resolution_range):
        parameters['resolution_parameter'] = resolution
        partition = find_partition(**parameters)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        partition = myfun.filter_modules(partition)
        profile.append(partition)
        # Ends loop.
# Ends If/else.

profile[0].summary()

#------------------------------------------------------------------------------
## Save clustering results.
#------------------------------------------------------------------------------

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
output_name = input_adjm.split("_")[1]
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "_partitions.csv")
DataFrame(results['Membership']).to_csv(myfile)

# Save partition profile summary data.
df = DataFrame.from_dict(results)
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "_profile.csv")
df.to_csv(myfile)

quit()

#--------------------------------------------------------------------
# Other method of creating a graph...
#--------------------------------------------------------------------

# Get the values as np.array, it's more convenenient.
node_names = adjm.columns
A = adjm.values
# Create graph, A.astype(bool).tolist() or (A / A).tolist() can also be used.
g = Graph.Adjacency((A > 0).tolist())
# Add edge weights and node labels.
g.es['weight'] = A[A.nonzero()]
g.vs['label'] = node_names  # or a.index/a.columns


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
