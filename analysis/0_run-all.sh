#!/usr/bin/env bash
# 0_run-all.sh - execute this script to run the entire analysis.

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
./1_Data-Preprocessing/data-preprocessing.R Cortex &> Cortex.report & spin

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
./1_Data-Preprocessing/data-preprocessing.R Striatum &> Striatum.report & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 1b passed.
else
	echo Failed at step 1b.
	exit 0
fi

# STEP 2.
echo "Combing datasets with TAMPOR."
./2_TAMPOR-Normalization/tampor-normalization.R 2>&1 | tee --append Cortex.report Striatum.report > /dev/null & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 2 passed.
else
	echo Failed at step 2.
	exit
fi

# STEP 3a.
echo "Performing Cortex network analysis."
./3_Network-Analysis/0_run-analysis.sh Cortex &>> Cortex.report & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 3a passed.
else
	echo Failed at step 3a.
	exit
fi

# STEP 3b.
echo "Performing Striatum network analysis."
./3_Network-Analysis/0_run-analysis.sh Striatum &>> Striatum.report & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo Step 3b passed.
else
	echo Failed at step 3b.
	exit
fi
