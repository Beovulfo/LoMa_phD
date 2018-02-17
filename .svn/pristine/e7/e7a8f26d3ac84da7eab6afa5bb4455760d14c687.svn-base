#!/bin/sh
#PBS -lwalltime=60:00:00
#PBS -lnodes=12
###PBS -q oldnodes
#PBS -N 3DTINCnewloma
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/newloma/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 12 /home/toni/SVNlast/branches/newloma/loma

echo "newLoMa is running... they are coming..."
