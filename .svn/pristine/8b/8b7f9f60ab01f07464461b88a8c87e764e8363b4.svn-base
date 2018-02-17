#!/bin/sh
#PBS -lwalltime=00:10:00
#PBS -lnodes=18
##PBS -q testing
#PBS -N 3Dinc
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/lomacte/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMa"
mpirun -np 18 /home/toni/SVNlast/branches/lomacte/loma

echo "LoMaCte is running... they are coming..."
