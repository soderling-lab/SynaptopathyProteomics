#!/usr/bin/env python3
' Clustering of the protein co-expression graph with Leidenalg.'

#------------------------------------------------------------------------------
## Load the adjacency matrix.
#------------------------------------------------------------------------------

# User parameters to change: 
data_type = 0 # Combined, KO, or WT network (0,1,2)
rmin = 0      # Min resolution
step = 1      # Step size
rmax = 100    # Max resolution

# Imports.
import os
import glob
from os.path import dirname
from sys import stderr
from pandas import read_csv

# Define a function to convert NoneType to blank string.
def xstr(s):
	if s is None:
		return ''
	return str(s)

# Get system variables.
vars = ['SLURM_JOBID','SLURM_CPUS_PER_TASK']
envars = {var:os.environ.get(var) for var in vars}
jobID = xstr(envars['SLURM_JOBID'])

# Which analysis are we doing?
geno = ['Combined','KO','WT'][data_type]
print("Performing Leiden algorithm clustering of the " + geno + 
        " protein co-expression network.", file = stderr)

# Read bicor adjacency matrix.
here = os.getcwd()
root = dirname(dirname(here))
datadir = os.path.join(root,"rdata")
myfiles = glob.glob(os.path.join(datadir,"*Adjm.csv"))
adjm = read_csv(myfiles[data_type], header = 0, index_col = 0)

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
g.add_edges(el)
g.es['weight'] = edges['weight'] 

# Remove self-loops.
g = g.simplify(multiple = False, loops = True)

#------------------------------------------------------------------------------
## Community detection with the Leiden algorithm.
#------------------------------------------------------------------------------

import leidenalg as la
from numpy import linspace
from leidenalg import find_partition
from leidenalg import CPMVertexPartition
from leidenalg import Optimiser
from progressbar import ProgressBar

# Loop to perform leidenalg community detection at n resolutions.
pbar = ProgressBar()
resolution_range = linspace(rmin,step,rmax)
profile = list()
for res in pbar(resolution_range):
    # Perfrom La clustering.
    partition = find_partition(g, CPMVertexPartition, 
            weights='weight', resolution_parameter=res)
    # Partition optimization... this takes extra time.
    optimiser = la.Optimiser()
    diff = optimiser.optimise_partition(partition,n_iterations=-1)
    # Add optimized partition to profile list.
    profile.append(partition)

#------------------------------------------------------------------------------
## Save partition profile results.
#------------------------------------------------------------------------------

from pandas import DataFrame

# Collect partition results. 
print(f"Complete! Examined network at {len(profile)} resolutions!", file = stderr)
results = {
        'Modularity' : [partition.modularity for partition in profile],
        'Membership' : [partition.membership for partition in profile],
        'Summary'    : [partition.summary() for partition in profile],
        'Resolution' : [partition.resolution_parameter for partition in profile]}

# Save cluster membership.
myfile = os.path.join(datadir, jobID + "3_" + geno + "_partitions.csv")
DataFrame(results['Membership']).to_csv(myfile)

# Save partition profile.
df = DataFrame.from_dict(results)
myfile = os.path.join(datadir, jobID + "3_" + geno + "_profile.csv")
df.to_csv(myfile)
