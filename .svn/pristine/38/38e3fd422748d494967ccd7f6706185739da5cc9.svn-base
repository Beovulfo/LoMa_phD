#!/bin/sh
#PBS -lwalltime=100:00:00
#PBS -lnodes=36
#PBS -q newnodes
#PBS -N hiroll01
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVN/LowMach/branches/mixing_layer/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 36 /home/toni/SVN/LowMach/branches/mixing_layer/loma

echo "LoMa finished"
