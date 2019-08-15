import leidenalg as la
import igraph as ig
import numpy as np

# Generate some graphs. simple, undirected, unweighted graphs
n = 100
g1 = ig.Graph.Lattice([n], 1)
g2 = ig.Graph.Lattice([n], 1)

# Perform multiplex community detection.
optimiser = la.Optimiser()
membership, improvement = la.find_partition_multiplex([g1,g2],la.ModularityVertexPartition)
print(f"Number of clusters: {np.unique(membership).size}")

# Calculate Modularity.
part1 = la.ModularityVertexPartition(g1, membership)
part2 = la.ModularityVertexPartition(g2, membership)

# Examine number of clusters.
part1.summary()
part2.summary()
print(f"Modularity (Q)    : {part1.modularity}")
print(f"Modularity (Q)    : {part2.modularity}")

# Optimization of a single slice.
p01 = la.find_partition(g1, la.ModularityVertexPartition, n_iterations = -1)
p01.summary()
p01.modularity

p02 = la.find_partition(g2, la.ModularityVertexPartition, n_iterations = -1)
p02.summary()
p02.modularity


