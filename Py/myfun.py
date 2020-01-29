#!/usr/bin/env Python3

#--------------------------------------------------------------------
# xstr
def xstr(s):
    ''' Convert NoneType to blank ('') string.'''
    if s is None:
        return ''
    else:
        return str(s)
# EOF

#--------------------------------------------------------------------
# contains
def contains(mylist,value,return_index=False):
    ''' Check if list contains a value. 
    Like list.index(value) but returns False if the provided list
    does not contain value. 
    '''
    list_as_dict = dict(zip(mylist,range(len(mylist))))
    index = list_as_dict.get(value)
    if not return_index: 
        return type(index) is int # True if in list.
    else:
        return index
# EOF

#--------------------------------------------------------------------
# filter_modules
def filter_modules(partition,min_size=5,unassigned=0):
    """ Set modules with size less than minimum size to 0. """
    import numpy as np
    membership = np.array(partition.membership)
    sizes = np.array(partition.sizes())
    remove = np.array(range(len(sizes)))[sizes < min_size]
    out = [node in remove for node in membership]
    membership[out] = unassigned
    partition.set_membership(membership)
    m = np.array(partition.membership)
    unassigned = sum(m==unassigned)/len(m)
    return partition
# EOF

#--------------------------------------------------------------------
# graph_from_adjm
def graph_from_adjm(adjm,weighted=True,signed=True):
    # Simplifing graph seems to mess things up.
    import numpy as np
    from igraph import Graph
    from pandas import DataFrame
    if not signed: adjm = abs(adjm)
    # Simplify graph by keeping only upper tri...
    # Melt upper triangle into edges df.
    #idx = np.triu(np.ones(adjm.shape)).astype(np.bool)
    #adjm = adjm.where(idx)
    edges = adjm.stack().reset_index()
    edges.columns = ['nodeA','nodeB','weight']
    edges = edges[edges.weight != 0]
    edge_tuples = list(zip(edges.nodeA,edges.nodeB,edges.weight))
    if weighted: g = Graph.TupleList(edge_tuples,weights=True)
    if not weighted: g = Graph.TupleList(edge_tuples,weights=False)
    return g

#--------------------------------------------------------------------
# add_method

def add_method(cls):
    ''' Add a method to an existing class.
    Example:
    @add_method(A)
    def foo():
        print('hello world!')
    '''
    from functools import wraps 
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator

#--------------------------------------------------------------------
# Function to apply a threshold to a graph such that it is a single component.
def apply_best_threshold(subg,start=0,stop=1,num=10):
    import numpy as np
    import igraph
    from pandas import DataFrame
    # Also needs myfun...
    # Remove multiple edges.
    subg = subg.simplify(multiple=False) # has names.
    nodes = subg.vs['name']
    i = 0
    is_connected = True
    # Loop to check if graph is single component after thresholding.
    while is_connected:
        threshold = np.linspace(start,stop,num)[i]
        adjm = np.array(subg.get_adjacency(attribute="weight").data)
        mask = adjm > threshold
        df = DataFrame(adjm * mask, index=nodes,columns=nodes)
        subg_filt = myfun.graph_from_adjm(df,weighted=True,signed=True)
        components = subg_filt.components()
        is_connected = len(set(components.membership)) == 1
        i +=1
    # Ends while loop.
    return subg_filt
# Ends function.

#--------------------------------------------------------------------
# A function to perform MCL clustering.
def clusterMCL(graph, inflation=1.2, weight='weight'):
    # FIXME:: How to return updated graph?
    import os
    import igraph
    import numpy as np
    from pandas import DataFrame
    edges = graph.get_edgelist()
    nodes = [graph.vs[edge]['name'] for edge in edges]
    weights = graph.es[weight]
    edge_list = [node + [edge] for node,edge in zip(nodes,weights)]
    DataFrame(edge_list).to_csv(".tempnet.csv",sep="\t",header=False,index=False)
    # Execute MCL. Send stderr to devnull. Pipe stdout back into python.
    cmd = ["mcl",".tempnet.csv","--abc","-I",str(inflation),"-o","-"]
    process = subprocess.Popen(cmd,stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL)
    os.remove(".tempnet.csv")
    # Parse the output.
    out = process.communicate()
    modules = list(out)[0].decode('utf-8').split("\n") # decode
    modules = [module.split("\t") for module in modules]
    modules = [module for module in modules if module != ['']] # remove empty
    # Create partition.
    k = len(modules)
    module_size = [len(module) for module in modules]
    module_names = list(range(1,k+1))
    partition = dict(zip(sum(modules, []), 
        np.repeat(module_names, module_size, axis=0)))
    clusters = graph.clusters()
    membership = [partition.get(protein) for protein in graph.vs['name']]
    # Set membership.
    clusters._membership = membership
    # How to update graph?
    return(clusters)
# Done.

#--------------------------------------------------------------------
# Wrapper around clusterMCL to find best inflation parameter.
def clusterMaxMCL(graph,inflation):
    mcl_clusters = list()
    Q = list()
    for i in inflation:
        result = clusterMCL(graph,inflation=i)
        mcl_clusters.append(result)
        Q.append(result.recalculate_modularity())
    # Done.
    idx = Q.index(max(Q))
    best_Q = Q[idx]
    best_i = inflation[idx]
    best_clusters = mcl_clusters[idx]
    print("Best Inflation : {}".format(best_i) +
            "\nBest Modularity: {}".format(best_Q))
    return(best_clusters)
# Done.
