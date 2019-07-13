# In order to set up HPO, we need to:
* Define a WGCNA function that can be called by Python.
* This function should take our parameters and bounds.
* Python scikit-optimize can be used to find the "best" params.
* We can min/max some quality statistic.
* Number of clusters could be optimized by minimizing distortion... the sum 
  of the squared distances to the clusters center. This idea is borrowed from
  k-means clustering. The best parameters are identified by inspection of the
  scree plot and finding its "elbow".
* Network quality can be evaluated as:
  Modularity: overall quality of the partition
  Module Cohesiveness: quality (PVE) of the individual modules. 
* Typically there is a trade off between these two statistics.

