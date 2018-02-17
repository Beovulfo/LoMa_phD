#!/bin/sh
#PBS -lwalltime=100:00:00
#PBS -lnodes=34
#PBS -q newnodes
#PBS -N mlRM01
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

echo "LoMa finished"
