# Code
-------------------------------------------------------------------------------
This directory contains the R and Python scripts used to perform the data 
analysis. The scripts are named in the order in which they should be executed.

* [1_Preprocessing/]() Preprocessing of the cortex and striatum data including 
simple normalization within batches, normalization across tissue specific batches 
with IRS normalization, and regression of genetic background with empirical bayes
regression.

* [2_TAMPOR/]() Preprocessed cortex and straitum datasets were combined 
using an iterative median-polish like approach termed TAMPOR. The output of this 
script is the combined cortex and striatum data. Data are scaled relative to WT 
samples across both tissue types. Thus only relative differences to WT samples
are retained.

* [3_WPCNA/]() 
