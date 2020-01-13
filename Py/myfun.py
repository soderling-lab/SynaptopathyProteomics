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
