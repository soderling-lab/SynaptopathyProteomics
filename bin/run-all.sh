#!/usr/bin/env bash

# Parse args.
analysis_type="$1"
here="$(echo $PWD)"
root="$(dirname $here)"

# Preprocess data.
$root/analysis/1_Data-Preprocessing/0_data-preprocessing.R $analysis_type

# Generate protein networks.
$root/analysis/2_Network-Analysis/0_create-networks.R $analysis_type

# Cluster networks.
$root/analysis/2_Network-Analysis/1_leidenalg-clustering.R $analysis_type

# Enforce module self-preservation.
$root/analysis/2_Network-Analysis/2_module-self-preservation.R $analysis_type

# Perform network analysis of modules.
$root/analysis/2_Network-Analysis/3_network-analysis.R $analysis_type
