#!/usr/bin/env python3

#------------------------------------------------------------------------------
## Load the input adjacency matrix.
#------------------------------------------------------------------------------

#input_adjm = "3_PPI_Adjm.csv" # Input adjacency matrix.
input_adjm = "3_GO_Semantic_Similarity_RMS_Adjm.csv"

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
## Explore Community detection methods.
#------------------------------------------------------------------------------

# Methods for positive edge weights, no resolution parameter:
# [1] ModularityVertexPartition -- Implements modularity. Only 
# [2] SurpriseVertexPartition -- Implements (asymptotic) Surprise. 

# Method for unweighted graph.
# [6] SignificanceVertexPartition -- Implements Significance. This 

# Methods for positive edges with linear resolution parameter:
# [3] RBConfigurationVertexPartition -- Implements Reichardt and 
# [4] RBERVertexPartition -- Implements Reichardt and Bornholdt’s 
# [5] CPMVertexPartition -- Implements CPM. Quality function is 

import numpy as np
from leidenalg import Optimiser
from leidenalg import find_partition
from progressbar import ProgressBar
from importlib import import_module

# Clustering methods.
m1 = ["ModularityVertexPartition",
        "SurpriseVertexPartition"]
m2 = ["RBConfigurationVertexPartition",
        "RBERVertexPartition",
        "CPMVertexPartition"]
m3 = ["SignificanceVertexPartition"]
methods = {"single" : m1, "multi" : m2, "unweighted" : m3}

# Loop to explore methods as defined by:
mtype = "unweighted"
resolution = 0
profile = list()
min_size = 5 # Mimimum module size. If smaller then set to 0.
for method in methods[mtype]:
    partition_type = getattr(import_module('leidenalg'), method)
    if mtype is "single":
        partition = find_partition(g, partition_type, weights='weight')
    elif mtype is "multi":
        partition = find_partition(g, partition_type, weights='weight',resolution_parameter=resolution)
    elif mtype is "unweighted":
        partition = find_partition(g, partition_type)
    optimiser = Optimiser()
    diff = optimiser.optimise_partition(partition,n_iterations=-1)
    membership = np.array(partition.membership)
    sizes = np.array(partition.sizes())
    remove = np.array(range(len(sizes)))[sizes < min_size]
    out = [protein in remove for protein in membership]
    membership[out] = 0
    partition.set_membership(membership)
    m = np.array(partition.membership)
    unassigned = sum(m==0)/len(m)
    print("{}: ".format(method) + partition.summary() + ".")
    print("Percent unassigned:{}".format(unassigned) + " (%).\n")
    profile.append(partition)
# Ends loop.


