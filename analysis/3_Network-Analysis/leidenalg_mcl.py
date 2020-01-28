#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

## FIXME: Significance method doesn't seem to be working.

## User parameters: 
# adjm_type - string specifying input adjacency matrix.
# method - string specifying the optimization method. One of: Modularity, 
#    Surprise, RBConfiguration, RBER, CPM, or Significance.
# rmin - min resolution
# rmax - max resolution
# nsteps - number of steps between rmin and rmax.
adjm_type = 'Cortex' 
method = 'CPM' 
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
        "Modularity": {'partition_type' : 'ModularityVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : -1},
        # Surprise
        "Surprise": {'partition_type' : 'SurpriseVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : -1},
        # RBConfiguration
        "RBConfiguration": {'partition_type' : 'RBConfigurationVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : -1},
        # RBER
        "RBER": {'partition_type' : 'RBERVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : -1},
        # CPM
        "CPM": {'partition_type' : 'CPMVertexPartition', 
            'weights' : True, 'signed' : True,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : -1},
        # Significance
        "Significance": 
        {'partition_type' : 'SignificanceVertexPartition', 
            'weights': None, 'signed' : False,
            'resolution_parameter' : None,
            'n_iterations' : -1}
        }

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
adjm = adjm.set_index(keys=adjm.columns)

# Create igraph graph -- this takes a few seconds.
if parameters.get('weights') is not None:
    # Create a weighted graph.
    g = myfun.graph_from_adjm(adjm,weighted=True,signed=parameters.pop('signed'))
    parameters['weights'] = 'weight' # update weights parameter.
else:
    # Create an unweighted graph.
    g = myfun.graph_from_adjm(adjm,weighted=False,signed=parameters.pop('signed'))

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

# Remove any None type parameters.
out = [key for key in parameters if parameters.get(key) is None]
for key in out: del parameters[key]

# Perform Leidenalg community detection. 
if parameters.get('resolution_parameter') is None:
    # Single resolution methods:
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
    # Loop to perform multi-resolution clustering methods:
    pbar = ProgressBar()
    profile = list()
    resolution_range = linspace(**parameters.get('resolution_parameter'))
    for resolution in pbar(resolution_range):

        resolution_range = linspace(**parameters.get('resolution_parameter'))
        resolution = resolution_range[0]
        parameters['resolution_parameter'] = resolution
        partition = find_partition(**parameters)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)

        # Recursively split large modules with MCL.
        # 1. get big modules.
        # 2. Threshold.
        # 3. MCL

        # Get modules that are too big.
        max_size = 500
        modules = set(partition.membership)
        too_big = [mod for mod in modules if partition.size(mod) > max_size]
        # Threshold graphs.
        subg = partition.subgraphs()
        subg = [apply_best_threshold(subg[i]) for i in too_big]

        # Write to file.
for graph in subg:

    edges = graph.get_edgelist()
    nodes = [graph.vs[edge]['name'] for edge in edges]
    weights = graph.es['weight']
    edge_list = [node + [edge] for node,edge in zip(nodes,weights)]
    DataFrame(edge_list).to_csv("network.csv",sep="\t",header=False,index=False)



# Continue...
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
# Ends if/else

# Save cluster membership vectors.
output_name = input_adjm.split("_")[1]
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "_partitions.csv")
df = DataFrame(results['Membership'])
df.columns = profile[0].graph.vs['name']
df.to_csv(myfile)

# Save partition profile summary data.
df = DataFrame.from_dict(results)
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "_profile.csv")
df.to_csv(myfile)
