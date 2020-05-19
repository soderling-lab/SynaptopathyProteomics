#!/usr/bin/env bash
# Run the analysis.

# Parse args.
if [ $# -eq 0 ]
then
	echo "Error: Specify 'Cortex' or 'Striatum'."
	exit 1
fi

# Analysis TYPE.
TYPE="$1"

# Project ROOT.
here="$(echo $PWD)"
ROOT="$(dirname $here)"

# Preprocess raw peptide data. Analyze
# protein differential abundance by utilizing
# a glm (background ~ genotype.treatment) 
# implemented from EdgeR.
$ROOT/analysis/0_data-preprocessing.R $TYPE

# Generate protein networks using bicor from 
# WGCNA and network enhancement from neten.
$ROOT/analysis/1_create-networks.R $TYPE

# Cluster networks with LeidenAlgorithm + 
# optimization with Surprise.
$ROOT/analysis/2_leidenalg-clustering.py $TYPE

# Enforce module self-preservation using NetRep.
$ROOT/analysis/3_module-self-preservation.R $TYPE

# Perform network analysis of modules. Uses
# ideas drawn from WGCNA's summarization 
# of modules as their eigengene.
$ROOT/analysis/4_network-analysis.R $TYPE
