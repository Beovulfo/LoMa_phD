#!/bin/sh
#PBS -lwalltime=30:00:00
###PBS -lnodes=compute-0-3:ppn=12+compute-0-4:ppn=12
#PBS -lnodes=10
###PBS -q newnodes
#PBS -N loma2inc
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/loma_transition2/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 10 /home/toni/SVNlast/branches/loma_transition2/loma

echo "LoMa is running... they are coming..."
