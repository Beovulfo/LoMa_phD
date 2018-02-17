#!/bin/sh
#PBS -lwalltime=100:00:00
#PBS -lnodes=12
##PBS -q newnodes
#PBS -N testicaro2
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/mixing_layer/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 12 /home/toni/SVNlast/branches/mixing_layer/loma

echo "LoMa finished"
