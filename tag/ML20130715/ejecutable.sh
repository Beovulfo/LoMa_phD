#!/bin/sh
#PBS -lwalltime=020:30:00
#PBS -lnodes=18
#PBS -q newnodes
#PBS -N mlBCNeuman1
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVN/LowMach/branches/mixing/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for liso"
mpirun -np 18 /home/toni/SVN/LowMach/branches/mixing/liso

echo "Liso running"
