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
rm -f "$TISSUE.report"

# STEP 1.
echo "Generating "$TISSUE" protein networks."
./1_Generate-networks.R "$TISSUE" &> "$TISSUE-network.report" & spin

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo -e "\tStep 1 passed."
else
	echo -e "\tFailed at step 1."
	exit 0
fi

# STEP 2.
echo "Clustering the "$TISSUE" protein co-variation network."
./2_Leiden-clustering.py "$TISSUE" &>> "$TISSUE-network.report" & spin

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo -e "\tStep 2 passed."
else
	echo -e "\tFailed at step 2."
	exit 0
fi

# STEP 3.
echo "Enforcing  "$TISSUE" module self-preservation."
./3_Module-preservation.R "$TISSUE" &>> "$TISSUE-network.report" & spin

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo -e "\tStep 3 passed."
else
	echo -e "\tFailed at step 3."
	exit 0
fi

# STEP 4.
echo "Checking "$TISSUE" modules for convergent changes." & spin
./4_Network-analysis.R "$TISSUE" &>> "$TISSUE-network.report"

# Check, did script run successfully?
if [ $? -eq 0 ]
then
	echo -e "\tStep 4 passed."
else
	echo -e "\tFailed at step 4."
	exit 0
fi

# STEP 5.
echo "Analyzing "$TISSUE" modules for GO enrichment."
./5_Module_GO-Enrichment.R "$TISSUE" &>> "$TISSUE-network.report" & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo -e "\tStep 5 passed."
else
	echo -e "\tFailed at step 5."
	exit 0
fi

# STEP 6.
echo "Testing "$TISSUE" modules for enrichment of DBD-associated genes."
./6_Module_DBD-Enrichment.R "$TISSUE" &>> "$TISSUE-network.report" & spin

# Check if completed successfully?
if [ $? -eq 0 ]
then
	echo -e "\tNetwork analysis complete!"
else
	echo -e "\tFailed at step 6."
	exit 0
fi
