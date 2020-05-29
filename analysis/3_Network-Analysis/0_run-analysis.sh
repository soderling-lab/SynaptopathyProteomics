#!/usr/bin/env bash

# 0_run-analysis.sh - execute this script to run the analysis.

# Check if an argument was passed.
if [ $# -eq 0 ]
then
	echo "Error: Provide tissue type for analysis: 'Cortex' or 'Striatum'."
	exit
fi

# Input: should be either 'Cortex' or 'Striatum'.
CORTEX="Cortex"
STRIATUM="Striatum"
TISSUE="$1"

# Remove existing reports.
rm -f ./0[1-7]*.txt

# STEP 1.
echo "Generating "$CORTEX" and "$STRIATUM" protein networks."
./1_Generate-networks.R &> 01_Generate-networks.txt

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Passed.
else
	echo Failed.
	exit
fi

# STEP 2.
echo "Clustering the "$TISSUE" protein co-variation network."
./2_Leiden-clustering.py "$TISSUE" &> 02_"$TISSUE"_Leiden-clustering.txt

# STEP 3.
echo "Enforcing  "$TISSUE" module self-preservation."
./3_Module-preservation.R "$TISSUE" &> 03_"$TISSUE"_Module-preservation.txt

# STEP 4.
echo "Checking "$TISSUE" modules for convergent changes."
./4_Network-analysis.R "$TISSUE" &> 04_"$TISSUE"_Network-analysis.txt

# STEP 5.
echo "Analyzing "$TISSUE" modules for GO enrichment."
./5_Module_GO-Enrichment.R "$TISSUE" &> 05_"$TISSUE"_Module_GO-Enrichment.txt

# STEP 6.
echo "Testing "$TISSUE" modules for enrichment of DBD-associated genes."
./6_Module_DBD-Enrichment.R "$TISSUE" &> 06_"$TISSUE"_Module-DBD-Enrichment.txt

# STEP 7.
#./7_Select_Module_Analysis.R "$TISSUE" &> 07_"$TISSUE"_Module-DBD-Enrichment.txt
