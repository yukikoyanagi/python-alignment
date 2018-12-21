#!/usr/bin/env bash
#
#SBATCH --account austmathjea_slim
#SBATCH --nodes 1               # number of nodes
#SBATCH --time 2:00:00            # max time (HH:MM:SS)
#SBATCH --mail-type END,FAIL
#
# File: slrm.alignlocpats.sh
#
# Description: This is the job script submitted by sbatch for
#  computing alignment scores for local patterns.
# See also alignlocpats.py.
#
# Author: Yuki Koyanagi
#

echo Running on "$(hostname)"
echo Running job: "$SLURM_JOB_NAME"
echo Available nodes: "$SLURM_NODELIST"
echo Slurm_submit_dir: "$SLURM_SUBMIT_DIR"
echo Start time: "$(date)"
start=$(date +%s)

echo Clearing SCRATCH folder
rm -f ${SCRATCH}/*

echo Writing hostnames
scontrol show hostnames $SLURM_NODELIST > /tmp/nodelist

echo Activate align environment
source activate align


echo Starting servers
srun ppserver.py -p 2048 -k 10800 -t 30 &
sleep 1 # sleep a bit to ensure that the servers have started

echo Starting Python program
alignlocpats.py patterns.pkl patrot.791.txt predictions.pkl \
		--window 10 --matrix bondlength --gap -1

end=$(date +%s)
echo End time: "$(date)"
dur=$(date -d "0 $end sec - $start sec" +%T)
echo Duration: "$dur"

sleep 35 # this should kill all remote ppservers

echo Done.
