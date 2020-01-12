#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

## User parameters: 
myfile = 0 # See input_adjm below.
method = 5 # See methods below.

#------------------------------------------------------------------------------
## Parse the user provided parameters.
#------------------------------------------------------------------------------

import sys
from sys import stderr

## Input adjacency matrix.
input_adjm = ["3_PPI_Adjm.csv","3_GO_Semantic_Similarity_RMS_Adjm.csv"][myfile]

## Leidenalg supports the following methods for optimization methods:
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
            'resolution_parameter' : {'start':0,'stop':1,'num':100}},
        # RBEVertex
        3: {'partition_type' : 'RBEVertexPartition', 
            'weights' : 'positive',
            'resolution_parameter' : {'start':0,'stop':1,'num':100}},
        # CPM
        4: {'partition_type' : 'CPMVertexPartition', 
            'weights' : 'positive and negative',
            'resolution_parameter' : {'start':0,'stop':1,'num':100}},
        # Significance
        5: {'partition_type' : 'SignificanceVertexPartition', 
            'weights':None,
            'resolution_parameter' : None}}

# Parameters for clustering.
parameters = methods.get(method)
method = parameters.get('partition_type')

# Status report.
print("Performing Leidenalg clustering utilizing the {}".format(method),
        "method to find optimal partition(s).", file=stderr)

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
myfile = os.path.join(datadir,input_adjm)
adjm = read_csv(myfile, header = 0, index_col = 0)
adjm = adjm.set_index(keys=adjm.columns) # Add row names.

# Incorporate weighted, positive, unweighted graph types.
if parameters['weights'] is 'positive':
    A = abs(adjm.values)
else:
    A = adjm.values

# Create igraph object.
g = Graph.Adjacency(A.tolist())  #g = Graph.Adjacency((A > 0).tolist()) # Unweighted or positive?
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
