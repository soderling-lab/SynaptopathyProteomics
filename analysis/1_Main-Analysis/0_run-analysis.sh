#!/usr/bin/env bash
# 0_run-analysis.sh - execute this script to run the analysis.

# Output log files:
CORTEX="Cortex.report"
STRIATUM="Striatum.report"

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

# A function that checks if script completed succcessfully.
check()  {
	if [ $? -eq 0 ]
	then
		echo -e "\tCompleted."
	else
		echo -e "\tError: Process failed."
		exit 0
	fi
}

# Remove any existing reports.
rm -f *.report

# What time is it?
now=$(date)
echo "Starting analysis at: $now."

# STEP 1.
echo "[1/8] Processing Cortex TMT data."
./1_*.R Cortex &> "$CORTEX" & spin 
check

# STEP 2.
echo "[2/8] Processing Striatum TMT data."
./1_*.R Striatum &> "$STRIATUM" & spin
check

# STEP 3.
echo "[3/8] Clustering Cortex protein covariation network."
./2_*.py Cortex &>> "$CORTEX" & spin
check

# STEP 4.
echo "[4/8] Clustering Striatum protein covariation network."
./2_*.py Striatum &>> "$STRIATUM" & spin
check

# STEP 5.
echo "[5/8] Enforcing Cortex module self-preservation."
./3_*.R Cortex &>> "$CORTEX" & spin
check

# STEP 6.
echo "[6/8] Enforcing Striatum module self-preservation."
./3_*.R Striatum &>> "$STRIATUM" & spin
check

# STEP 7.
echo "[7/8] Analyzing Cortex modules for differential abundance."
./4_*.R Cortex &>> "$CORTEX" & spin
check

# STEP 8.
echo "[8/8] Analyzing Striatum modules for differential abundance."
./4_*.R Striatum &>> "$STRIATUM" & spin
check

now=$(date)
echo "Analysis completed at: $now."
