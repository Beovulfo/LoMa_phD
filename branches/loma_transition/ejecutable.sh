#!/bin/sh
#PBS -lwalltime=80:00:00
###PBS -lnodes=compute-0-3:ppn=12+compute-0-4:ppn=12
#PBS -lnodes=36
###PBS -q newnodes
#PBS -N incRESx2KH02
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/loma_transition/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 36 /home/toni/SVNlast/branches/loma_transition/loma

echo "LoMa is running... they are coming..."
