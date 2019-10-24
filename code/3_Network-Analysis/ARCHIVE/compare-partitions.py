#!/usr/bin/env python3

from pandas import read_csv

# Load the WT(1) and KO(2) partitions.
p1 = read_csv('WT_partitions.csv',index_col=0)
p2 = read_csv('KO_partitions.csv',index_col=0)

# Compare partitions.
from sklearn.metrics import fowlkes_mallows_score

# At every resolution, evalauate the similarity of the WT and KO partitions.
fmi = [fowlkes_mallows_score(p1.iloc[i], p2.iloc[i]) for i in range(100)]

# Which resolutions are most divergent?
fmi.index(min(fmi))

# ENDOFILE
