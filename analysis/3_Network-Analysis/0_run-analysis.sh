#!/usr/bin/env bash

# 0_run-analysis.sh - execute this script to run the analysis.
# Input: should be either 'Cortex' or 'Striatum'.



# Check if an argument was passed.
if [ $# -eq 0 ]
then
	echo "Please specify the tissue for analysis: 'Cortex' or 'Striatum'."
	exit 
fi
TISSUE="$1"
	
# Check if input was Cortex or Striatum.
if [ "$TISSUE" != "Cortex" ] && [ "$TISSUE" != "Striatum" ]
then
	echo "Input should be either Cortex or Striatum (case-sensitive)."
	exit
fi

# Remove any existing reports.
rm -f ./0[1-7]*.txt

# STEP 1.
echo "Generating Cortex and Striatum protein networks."
./1_Generate-networks.R &> 01_Generate-networks.txt

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo Step 1 passed.
else
	echo Failed at step 1.
	exit
fi

# STEP 2.
echo "Clustering the "$TISSUE" protein co-variation network."
./2_Leiden-clustering.py "$TISSUE" &> 02_"$TISSUE"_Leiden-clustering.txt

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo Step 2 passed.
else
	echo Failed at step 2.
	exit
fi

# STEP 3.
echo "Enforcing  "$TISSUE" module self-preservation."
./3_Module-preservation.R "$TISSUE" &> 03_"$TISSUE"_Module-preservation.txt

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo Step 3 passed.
else
	echo Failed at step 3.
	exit
fi

# STEP 4.
echo "Checking "$TISSUE" modules for convergent changes."
./4_Network-analysis.R "$TISSUE" &> 04_"$TISSUE"_Network-analysis.txt

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo Step 4 passed.
else
	echo Failed at step 4.
	exit
fi

# STEP 5.
echo "Analyzing "$TISSUE" modules for GO enrichment."
./5_Module_GO-Enrichment.R "$TISSUE" &> 05_"$TISSUE"_Module_GO-Enrichment.txt

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 5 passed.
else
	echo Failed at step 5.
	exit
fi

# STEP 6.
echo "Testing "$TISSUE" modules for enrichment of DBD-associated genes."
./6_Module_DBD-Enrichment.R "$TISSUE" &> 06_"$TISSUE"_Module-DBD-Enrichment.txt

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 6 passed.
else
	echo Failed at step 6.
	exit
fi

# STEP 7.
#./7_Select_Module_Analysis.R "$TISSUE" &> 07_"$TISSUE"_Module-DBD-Enrichment.txt

# Check if completed successfully?
#if [ $? -eq 0 ]
#then
#	echo Step 7 passed. Well done comrade.
#else
#	echo Failed at step 7.
#	exit
#fi

# Combine reports.
cat ./0[1-9]*.txt >> "$TISSUE"_Report.txt
