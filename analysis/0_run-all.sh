#!/usr/bin/env bash

# 0_run-all.sh - execute this script to run the entire analysis.
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
rm -f Cortex.report
rm -f Striatum.report

# STEP 1a.
echo "Processing raw Cortex data."
./1_Data-Preprocessing/data-preprocessing.R Cortex &> Cortex.report

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 1a passed.
else
	echo Failed at step 1a.
	exit 0
fi

# STEP 1b.
echo "Processing raw Striatum data."
./1_Data-Preprocessing/data-preprocessing.R Striatum &> Striatum.report

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 1b passed.
else
	echo Failed at step 1b.
	exit 0
fi

# STEP 2.
echo "Combing datasets with TAMPOR normalization."
./2_TAMPOR-Normalization/tampor-normalization.R 2>&1 | tee --append Cortex.report Striatum.report > /dev/null

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 2 passed.
else
	echo Failed at step 2.
	exit
fi

# STEP 3.
echo "Performing network analysis."
./3_Network-Analysis/0_run-analysis.sh "$TISSUE" &>> "$TISSUE.report"

# NOTE:
# To redirect stderr and stdout to file AND console:
#./1_Data-Preprocessing/data-preprocessing.R Cortex 2>&1 >>log.txt | tee --append log.txt
