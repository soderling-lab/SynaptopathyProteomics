#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

## User parameters: the only thing you need to change.
analysis_type = "Cortex" # Cortex or Striatum

## Optional parameters:
rmin = 0
rmax = 1
nsteps = 100
max_size = 100
recursive = True
n_iterations = -1
method = 'Surprise' 

## Description of parameters:
# rmin - minimum resolution for multi-resolution methods.
# rmax - maximum resolution for multi-resolution methods.
# adjm_type - string specifying input adjacency matrix.
# max_size - number of steps to take between rmin and rmax.
# recursive - boolean indicating whehter or not large modules will be
#    recursively split.
# n_iterations - number indicating the number of iterations to run Leidenalg
#    optimizer. If -1, then it will be run until there is no further 
#    improvement in the partition of the graph.
# method - string specifying the optimization method to be used.
# adjm_type - string specifying the input adjacency matrix, see below.

#--------------------------------------------------------------------
## Parse the user provided parameters.
#--------------------------------------------------------------------

import sys
from sys import stderr

## Input adjacency matrix.
adjms = {"Cortex" : "Cortex_Adjm.csv",
        "Enhanced Cortex" : "Cortex_NE_Adjm.csv",
        "Striatum" : "Striatum_Adjm.csv",
        "Enhanced Striatum" : "Striatum_NE_Adjm.csv",
        "Combined" : "Combined_Adjm.csv",
        "PPI" : "PPI_Adjm.csv"
        }

## Leidenalg supports the following optimization methods:
methods = {
        # Modularity
        "Modularity": {'partition_type' : 'ModularityVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations},
        # Surprise
        "Surprise": {'partition_type' : 'SurpriseVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations},
        # RBConfiguration
        "RBConfiguration": {'partition_type' : 'RBConfigurationVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations},
        # RBER
        "RBER": {'partition_type' : 'RBERVertexPartition', 
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations},
        # CPM
        "CPM": {'partition_type' : 'CPMVertexPartition', 
            'weights' : True, 'signed' : True,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations},
        # Significance
        # FIXME: Significance method doesn't seem to be working.
        "Significance": 
        {'partition_type' : 'SignificanceVertexPartition', 
            'weights': None, 'signed' : False,
            'resolution_parameter' : None,
            'n_iterations' : n_iterations}
        }

# Parameters for clustering.
adjm_type = ' '.join(['Enhanced',analysis_type])
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
import glob
from os.path import dirname

# Directories.
here = os.getcwd()
root = dirname(dirname(here))
datadir = os.path.join(root,"rdata")
funcdir = os.path.join(root,"Py")

# Load functions.
sys.path.append(root)
from Py.myfun import *

# Get system variables.
#myvars = ['SLURM_JOBID','SLURM_CPUS_PER_TASK']
#envars = {var:os.environ.get(var) for var in myvars}
#jobID = xstr(envars['SLURM_JOBID'])

#---------------------------------------------------------------------
## Load input adjacency matrix and create an igraph object.
#---------------------------------------------------------------------

from igraph import Graph
from pandas import read_csv, DataFrame

# Load adjacency matrix.
input_adjm = adjms.get(adjm_type)
myfile = os.path.join(datadir,input_adjm)
adjm = read_csv(myfile, header = 0, index_col = 0)
adjm = adjm.set_index(keys=adjm.columns)

# Create igraph graph. Note, this takes several moments.
if parameters.get('weights') is not None:
    # Create a weighted graph.
    g = graph_from_adjm(adjm,weighted=True,signed=parameters.pop('signed'))
    parameters['weights'] = 'weight'
else:
    # Create an unweighted graph.
    g = graph_from_adjm(adjm,weighted=False,signed=parameters.pop('signed'))

# Add graph to input parameters for clustering.
parameters['graph'] = g

#--------------------------------------------------------------------
## Community detection with the Leiden algorithm.
#--------------------------------------------------------------------

import numpy as np
from numpy import linspace
from importlib import import_module
from progressbar import ProgressBar
from leidenalg import Optimiser, find_partition

# Update partition type parameter.
# Dynamically load the partition_type class. This is the method to be used for
# clusering optimization.
parameters['partition_type'] = getattr(import_module('leidenalg'),method)

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
    profile.append(partition)
    if recursive:
        # Recursively split modules that are too big.
        subgraphs = partition.subgraphs()
        too_big = [subg.vcount() > max_size for subg in subgraphs]
        n_big = sum(too_big)
        print("Recursively spliting {} large module(s)...".format(n_big),
                file=stderr)
        print("... Initial partition: " + partition.summary() + ".",file=stderr)
        while any(too_big):
            # Perform clustering for any subgraphs that are too big.
            idx = [i for i, too_big in enumerate(too_big) if too_big] 
            parameters['graph'] = subgraphs.pop(idx[0])
            part = find_partition(**parameters)
            optimiser = Optimiser()
            diff = optimiser.optimise_partition(part,n_iterations=-1)
            # Add to list.
            subgraphs.extend(part.subgraphs())
            too_big = [subg.vcount() > max_size for subg in subgraphs]
        # Collect subgraph membership as a single partition.
        nodes = [subg.vs['name'] for subg in subgraphs]
        parts = [dict(zip(n,[i]*len(n))) for i, n in enumerate(nodes)]
        new_partition = {k: v for d in parts for k, v in d.items()}
        # Set membership of initial graph.
        new_membership = [new_partition.get(node) for node in partition.graph.vs['name']]
        partition.set_membership(new_membership)
        # Replace partition in profile list.
        profile[0] = partition
        print("... Final partition: " + partition.summary() + ".",file=stderr)
else:
    # Loop to perform multi-resolution clustering methods.
    pbar = ProgressBar()
    profile = list()
    resolution_range = linspace(**parameters.get('resolution_parameter'))
    for resolution in pbar(resolution_range):
        # Update resolution parameter.
        parameters['resolution_parameter'] = resolution
        partition = find_partition(**parameters)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,n_iterations=-1)
        profile.append(partition)
        # Ends loop.
# Ends If/else.

#------------------------------------------------------------------------------
## Save Leidenalg clustering results.
#------------------------------------------------------------------------------

# Collect partition results and save as csv. 
if len(profile) == 1:
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

# Save cluster membership (partition).
output_name = "_".join([analysis_type,"NE",method]) + ".csv"

myfile = os.path.join(datadir,  output_name)
df = DataFrame(results['Membership'])
df.columns = profile[0].graph.vs['name']
df.to_csv(myfile)
