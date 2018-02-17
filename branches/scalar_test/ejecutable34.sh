#!/bin/sh
#PBS -lwalltime=100:30:00
#PBS -lnodes=34
###PBS -q oldnodes
#PBS -N mlturb34
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVN/LowMach/branches/mixing_exp/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 34 /home/toni/SVN/LowMach/branches/mixing_exp/loma

echo "LoMa running"
