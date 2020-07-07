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
rm -f *.report

# Output log files:
CORTEX="Cortex.report"
STRIATUM="Striatum.report"

# STEP 1a.
echo "Processing raw Cortex data."
./1_*.R Cortex &> "$CORTEX" & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo -e "\tSuccessfully processed Cortex data."
else
	echo -e "\tError: Failed to process Cortex data."
	exit 0
fi

# STEP 1b.
echo "Processing raw Striatum data."
./1_*.R Striatum &> "$STRIATUM" & spin


# STEP 2.
echo "Clustering co-expression networks."
./2_*.py Cortex &> "$CORTEX" & spin
./2_*.py Striatum &> "$STRIATUM" & spin

# STEP 3.
echo "Enforcing module self-preservation."
./3_*.R Cortex &> "$CORTEX" & spin
./3_*.R Striatum &> "$STRIATUM" & spin

# STEP 4.
echo "Analyzing modules for differential abundance in Shank2, Shank3, Syngap1 and Ube3a groups."
./4_*.R Cortex &> "$CORTEX" & spin
./4_*.R Striatum &> "$STRIATUM" & spin
