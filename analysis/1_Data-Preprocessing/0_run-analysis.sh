#!/usr/bin/env bash
# 0_run-analysis.sh - execute this script to preprocess the data.

# Define a simple progres spinner.
spin() {
	# From William Pursell: https://stackoverflow.com/questions/12498304/
	pid=$! # Process ID of previously executed command.
	spin='-\|/'
	i=0
	while kill -0 $pid 2>/dev/null
	do
		i=$(( (i+1) %4 ))
		printf "\r${spin:$i:1}"
		sleep 0.1 # Can be adjusted.
	done
}

# Remove any existing reports.
rm -f Cortex.report
rm -f Striatum.report

# STEP 1a.
echo "Processing raw Cortex data."
./1_data-preprocessing.R Cortex &> Cortex.report & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Successfully processed Cortex data.
else
	echo Failed to process Cortex data.
	exit 0
fi

# STEP 1b.
echo "Processing raw Striatum data."
./1_data-preprocessing.R Striatum &> Striatum.report & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Successfully processed Striatum data.
else
	echo Failed to process Striatum data.
	exit 0
fi

# STEP 2.
echo "Combing datasets with TAMPOR."
./2_TAMPOR-normalization.R 2>&1 | tee --append Cortex.report Striatum.report > /dev/null & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Data preprocessing is complete!
else
	echo Failed to combine datasets with TAMPOR.
	exit
fi