#!/bin/sh
#PBS -lwalltime=100:00:00
#PBS -lnodes=34
####PBS -q newnodes
#PBS -q oldnodes
#PBS -N inc01test
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVN/LowMach/branches/scalar_test/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 34 /home/toni/SVN/LowMach/branches/scalar_test/loma

echo "LoMa finished"
