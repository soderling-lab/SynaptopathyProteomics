#!/usr/bin/env python3

#------------------------------------------------------------------------------
## Multiresolution partitioning of the protein co-expression network using the 
#  Leiden algorithm.
#------------------------------------------------------------------------------

import os
from igraph import Graph
import leidenalg as la

# Load ppi and pce graphs from file.
here = os.getcwd()
os.chdir(here)

g1 = Graph.Read_Adjacency("bicor_adjm.csv", sep=",", attribute="weight") 
g2 = Graph.Read_Adjacency("ppi_adjm.csv", sep=",", attribute="weight") 

# Try modularity calculation on single graphs.
p1 = la.find_partition(g1, la.CPMVertexPartition, n_iterations = -1)
p2 = la.find_partition(g2, la.ModularityVertexPartition, n_iterations = -1)
p1.summary()
p2.summary()

# How to weight the layers????
optimiser = la.Optimiser()
diff = optimiser.optimise_partition_multiplex([p1,p2], layer_weights=None, n_iterations=-1)

# Multilayer optimization
optimiser = la.Optimiser()
membership, improvement = la.find_partition_multiplex([g1,g2], la.CPMVertexPartition)


diff = optimise_partition_multiplex([part1,part2])
