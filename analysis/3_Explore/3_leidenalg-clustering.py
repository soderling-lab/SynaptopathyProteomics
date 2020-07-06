#!/usr/bin/env python3
'''
title: Leidenalg Clustering
description: clustering the protein network with Leidenalg(Q=CPM)
authors: Tyler W A Bradshaw
'''

## Parameters for multiresolution methods:
rmin = 0 # Min resolution for multi-resolution methods.
rmax = 1 # Max resolution for multi-resolution methods.
nsteps = 100 # Number of steps to take between rmin and rmax.
max_size = 100 # Maximum allowable size of a module.

## General optimization methods:
optimization_method = 'CPM'
n_iterations = -1  # The number of optimization iterations.

## Input data:
# Input adjacency matrix should be in root/rdata/
adjm_file = 'Cortex_NE_Adjm.csv'

## Output:
# Saved in root/rdata/
# [output_name]_partitions.csv
# FIXME: should pass input arg cortex and striatum.
# FIXME: should get method (eg CPM) from args.
output_name = 'TMT_Cortex_CPM' # Prefix out output partition, saved as .csv.

#------------------------------------------------------------------------------
## Prepare the workspace.
#------------------------------------------------------------------------------

# Imports.

import os
import sys
import glob
from sys import stderr
from os.path import dirname

import numpy as np
from numpy import linspace
from importlib import import_module
from progressbar import ProgressBar
from leidenalg import Optimiser, find_partition

from igraph import Graph
from pandas import read_csv, DataFrame

# Directories.
here = os.getcwd()
root = dirname(dirname(here))
rdatdir = os.path.join(root,"rdata")
funcdir = os.path.join(root,"Py")

# Load user defined functions.
sys.path.append(root)
from Py.myfun import *

# Leidenalg supports the following optimization methods:
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

# Get method specific parameters for clustering.
parameters = methods.get(optimization_method)
method = parameters.get('partition_type')

# Status report.
print("Performing Leidenalg clustering utilizing the {}".format(method),
        "method to find optimal partition(s)...", file=stderr)

#---------------------------------------------------------------------
## Load input adjacency matrix and create an igraph object.
#---------------------------------------------------------------------

# Load graph adjacency matrix.
myfile = os.path.join(rdatdir,adjm_file) 
adjm = read_csv(myfile, header = 0, index_col = 0) 

# Create igraph graph object and add to parameters dictionary. 
# NOTE: uses a function in pyfun to read a graph as adjacency matrix.
# Note, this can take several minutes.
if parameters.get('weights') is not None:
    # Create a weighted graph.
    g = graph_from_adjm(adjm,weighted=True,signed=parameters.pop('signed'))
    parameters['weights'] = 'weight'
    parameters['graph'] = g
else:
    # Create an unweighted graph.
    g = graph_from_adjm(adjm,weighted=False,signed=parameters.pop('signed'))
    parameters['graph'] = g

#--------------------------------------------------------------------
## Leidenalg Community Detection
#--------------------------------------------------------------------

# Update partition type parameter.
# Dynamically load the partition_type class. 
# This is the method to be used to optimize the clustering.
parameters['partition_type'] = getattr(import_module('leidenalg'),method)

# Remove any None type parameters.
out = [key for key in parameters if parameters.get(key) is None]
for key in out: del parameters[key]

# Multi-resolution methods:
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

#------------------------------------------------------------------------------
## Save clustering results.
#------------------------------------------------------------------------------

# Multi-resolution profile:
results = {
    'Modularity' : [partition.modularity for partition in profile],
    'Membership' : [partition.membership for partition in profile],
    'Summary'    : [partition.summary() for partition in profile],
    'Resolution' : [partition.resolution_parameter for partition in profile]}

# Save cluster membership vectors.
myfile = os.path.join(rdatdir, output_name + "_partition.csv")
df = DataFrame(results['Membership'])
df.columns = profile[0].graph.vs['name']
df.to_csv(myfile)

# Save other information about each partition.
results['Resolution']
results['Modularity']
