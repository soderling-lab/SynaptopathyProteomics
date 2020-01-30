#!/usr/bin/env python3
' Breaking down large communities with MCL.'

## Parameters for MCL clustering.
adjm_type = 'Cortex' 
max_size = 500 # maximum allowable size of a module.
i_min = 1.2 # Min inflation parameter.
i_max = 5 # Max inflation parameter.
nsteps = 10 # Number of steps between min and max.

# Input adjacency matrix.
adjms = {"Cortex" : "3_Cortex_Adjm.csv",
        "Striatum" : "3_Striatum_Adjm.csv",
        "Combined" : "3_Combined_Adjm.csv",
        "PPI" : "3_PPI_Adjm.csv",
        "GO" : "3_GO_Semantic_Similarity_RMS_Adjm.csv"}

# Imports.
import os
import sys
import glob
import igraph
import leidenalg
import numpy as np
from igraph import Graph
from numpy import linspace
from os.path import dirname
from pandas import read_csv, DataFrame

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

# Read bicor adjacency matrix.
input_adjm = adjms.get(adjm_type)
myfile = os.path.join(datadir,input_adjm)
adjm = read_csv(myfile, header = 0, index_col = 0)
adjm = adjm.set_index(keys=adjm.columns)

# Create a weighted graph.
graph = graph_from_adjm(adjm,weighted=True,signed=True)

# Load La partitions.
partition_files = {"Cortex":"147731383","Striatum":"148436673"}
myfile = glob.glob(os.path.join(datadir,
    partition_files.get(adjm_type) + "*"))[0]
partitions = read_csv(myfile,index_col=0)

# Collect all La partitions as list of dicts.
nrows = partitions.shape[0]
la_partitions = [dict(zip(partitions.iloc[i].keys(), # Node names.
    partitions.iloc[i].values)) for i in range(nrows)] # Int - node membership.

####################################
## EXAMPLE MCL.
####################################

# Add la partitions to graph.
la_partition = la_partitions[97]
la_clusters = leidenalg.VertexPartition.CPMVertexPartition(graph)
membership = [la_partition.get(node) for node in la_clusters.graph.vs['name']]
la_clusters.set_membership(membership)

# Loop to check sizes.
for la_partition in la_partitions:
    la_clusters = leidenalg.VertexPartition.CPMVertexPartition(graph)
    membership = [la_partition.get(node) for node in la_clusters.graph.vs['name']]
    la_clusters.set_membership(membership)
    print(sum([x>500 for x in set(la_clusters.sizes())]))

# Get a subgraph.
subg = la_clusters.subgraph(0)
cutoff = 0.2
# Loop to test multiple inflation values.
inflation = linspace(1.2,5,10)
for i in range(len(inflation)):
    # Apply thresholding.
    #subg = thresholdGraph(subg,quiet=False)
    subg.es.select(weight_lt=cutoff).delete()
    # Perform mcl.
    mcl_clusters = clusterMCL(subg,inflation[i],quiet=True)
    # Filter small modules.
    mcl_clusters = filterModules(mcl_clusters)
    # Summary.
    pc = sum([m!=0 for m in mcl_clusters.membership])/n
    print("Inflation: {}.".format(round(inflation[i],3)))
    print("MCL partition: {}".format(mcl_clusters.summary()))
    print("Percent clustered: {}".format(pc))
    print("Modularity: {}".format(mcl_clusters.modularity))
    print("\n")
# Ends loop.

####################################

# Loop to perform MCL clustering.
inflation = linspace(i_min,i_max,nsteps) # inflation space to explore.
mcl_partitions = list()
for resolution in range(len(la_partitions)):
    print("Working on resolution {} ...".format(resolution))
    # Define La partition.
    partition = la_partitions[resolution]
    membership = [partition.get(node) for node in graph.vs['name']]
    la_clusters = leidenalg.VertexPartition.CPMVertexPartition(graph)
    la_clusters.set_membership(membership)
    print("Initial partition: " + la_clusters.summary())
    ## Get modules that are too big.
    modules = set(membership)
    too_big = [mod for mod in modules if la_clusters.size(mod) > max_size]
    print("... Breaking down {} large module(s) with MCL ...".format(len(too_big)))
    ## Threshold graphs.
    subg = la_clusters.subgraphs()
    subg = [apply_best_threshold(subg[i]) for i in too_big]
    ## Perform modularity optimized MCL clustering.
    best_clusters = list()
    for g in subg:
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
    part0.update(part2)
    final_partition = part0
    # Update LA partition with MCL partition.
    mcl_partitions.append(final_partition)
    new_membership = [final_partition.get(node) for node in graph.vs['name']]
    la_clusters.set_membership(new_membership)
    # Remove small modules.
    la_clusters = filter_modules(la_clusters)
    print("Final partition: {}.".format(la_clusters.summary()))
    nc = sum([m != 0 for m in la_clusters.membership])/len(la_clusters.graph.vs)
    print("Fraction not-clustered: {}.\n".format(round(nc,3)))
# Ends loop.
