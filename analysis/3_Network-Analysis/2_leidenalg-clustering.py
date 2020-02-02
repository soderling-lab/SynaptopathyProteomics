#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

## User parameters: 
# adjm_type - string specifying input adjacency matrix.
# method - string specifying the optimization method. One of: Modularity, 
#    Surprise, RBConfiguration, RBER, CPM, or Significance.
# rmin - min resolution
# rmax - max resolution
# nsteps - number of steps between rmin and rmax.
adjm_type = 'Cortex_NE' 
method = 'CPM' 
n_iterations = -1
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
        "Cortex_NE" : "3_Cortex_NEAdjm.csv",
        "Striatum" : "3_Striatum_Adjm.csv",
        "Combined" : "3_Combined_Adjm.csv",
        "PPI" : "3_PPI_Adjm.csv",
        "GO" : "3_GO_Semantic_Similarity_RMS_Adjm.csv"}

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
myvars = ['SLURM_JOBID','SLURM_CPUS_PER_TASK']
envars = {var:os.environ.get(var) for var in myvars}
jobID = xstr(envars['SLURM_JOBID'])

#------------------------------------------------------------------------------
## Load input adjacency matrix and create an igraph object.
#------------------------------------------------------------------------------

from igraph import Graph
from pandas import read_csv, DataFrame

# Read bicor adjacency matrix.
input_adjm = adjms.get(adjm_type)
myfile = os.path.join(datadir,input_adjm)
adjm = read_csv(myfile, header = 0, index_col = 0)
adjm = adjm.set_index(keys=adjm.columns)

# Create igraph graph -- this takes several moments.
if parameters.get('weights') is not None:
    # Create a weighted graph.
    g = graph_from_adjm(adjm,weighted=True,signed=parameters.pop('signed'))
    parameters['weights'] = 'weight'
else:
    # Create an unweighted graph.
    g = graph_from_adjm(adjm,weighted=False,signed=parameters.pop('signed'))

# Add graph to input parameters for clustering.
parameters['graph'] = g

#------------------------------------------------------------------------------
## Community detection with the Leiden algorithm.
#------------------------------------------------------------------------------

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
    #partition = filter_modules(partition)
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
        #partition = filter_modules(partition)
        profile.append(partition)

        # Ends loop.
# Ends If/else.

#------------------------------------------------------------------------------
## Save Leidenalg clustering results.
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

#--------------------------------------------------------------------
## Decompose large communities with MCL.
#--------------------------------------------------------------------

## Parameters for MCL clustering.
max_size = 500 # maximum allowable size of a module.
i_min = 1.2 # Min inflation parameter.
i_max = 5 # Max inflation parameter.
nsteps = 10 # Number of steps between min and max.

# Collect all leidenalg partitions as list of dicts.
la_partitions = [dict(zip(profile[res].graph.vs['name'],
    profile[res].membership)) for res in range(len(profile))]

# Loop to perform MCL clustering.
inflation = linspace(i_min,i_max,nsteps) # inflation space to explore.
mcl_profile = profile.copy()
for resolution in range(len(mcl_profile)):
    print("Working on resolution {} ...".format(resolution))
    print("... Initial clustering result: " + profile[resolution].summary() + ".")
    ## Get modules that are too big.
    partition = mcl_profile[resolution]
    graph = partition.graph
    modules = set(partition.membership)
    too_big = [mod for mod in modules if partition.size(mod) > max_size]
    print("... Resolving {} large module(s) ...".format(len(too_big)))
    ## Threshold graphs.
    # FIXME: THIS IS SLOW!
    subg = partition.subgraphs()
    subg = [apply_best_threshold(subg[i]) for i in too_big]
    ## Perform modularity optimized MCL clustering.
    # FIXME: speed up by adding cluster parameter to MCL function!
    best_clusters = list()
    for g in subg:
        print("...")
        result = clusterMaxMCL(g, inflation) # result is clusters object.
        best_clusters.append(result)
    ## Combine MCL partitions into single partition.
    nodes = [part.graph.vs['name'] for part in best_clusters]
    parts = [part.membership for part in best_clusters]
    # Fix MCL membership indices.
    n = 1
    while n < len(parts):
        parts[n] = parts[n] + max(parts[n-1])
        n += 1
    # Combine as list of dicts.
    comb_parts = [dict(zip(nodes[i],parts[i])) for i in range(len(nodes))]
    # Flatten list.
    mcl_partition = {k: v for d in comb_parts for k, v in d.items()} 
    ## Combine MCL partition with initial LA parition...
    # First, make sure that membership indices are unique.
    part0 = la_partitions[resolution]
    part1 = mcl_partition
    n = max(set(part0.values()))
    part2 = dict(zip(part1.keys(),[x + n for x in part1.values()]))
    # Update LA partition with MCL partition.
    part0.update(part2)
    ## Set graph membership.
    all_nodes = profile[resolution].graph.vs['name']
    mcl_profile[resolution].set_membership([part0[node] for node in all_nodes])
    # Remove small modules.
    mcl_profile[resolution] = filter_modules(mcl_profile[resolution])
    # Summary:
    nc = sum([x==0 for x in mcl_profile[resolution].membership])/len(all_nodes) 
    k = sum([x != 0 for x in set(mcl_profile[resolution].membership)])
    print("Final partition, number of cluster: {}".format(k))
    print("Percent not-clustered: {}".format(round(nc,3)))
# Ends loop.

# Collect partition results and save as csv. 
if len(mcl_profile) is 1:
    # Single resolution profile:
    results = {
            'Modularity' : [partition.modularity for partition in mcl_profile],
            'Membership' : [partition.membership for partition in mcl_profile],
            'Summary'    : [partition.summary() for partition in mcl_profile]}
else: 
    # Multi-resolution profile:
    results = {
        'Modularity' : [partition.modularity for partition in mcl_profile],
        'Membership' : [partition.membership for partition in mcl_profile],
        'Summary'    : [partition.summary() for partition in mcl_profile],
        'Resolution' : [partition.resolution_parameter for partition in mcl_profile]}
# Ends if/else

# Save cluster membership vectors.
output_name = input_adjm.split("_")[1]
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "MCL_partitions.csv")
df = DataFrame(results['Membership'])
df.columns = mcl_profile[0].graph.vs['name']
df.to_csv(myfile)

# Save partition profile summary data.
df = DataFrame.from_dict(results)
myfile = os.path.join(datadir, jobID + "3_" + output_name + "_" + 
        method + "_MCL_profile.csv")
df.to_csv(myfile)
