#!/bin/sh
#PBS -lwalltime=50:05:00
#PBS -lnodes=24
####PBS -q oldnodes
#PBS -N mixing_expTURB
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/mixing_exp/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for mixing_exp"
mpirun -np 24 /home/toni/SVNlast/branches/mixing_exp/loma

echo "mixing_exp finished"
