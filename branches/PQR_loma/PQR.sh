#!/bin/sh
#PBS -q mlx
#PBS -l walltime=10:00:00
#PBS -N PQR
#PBS -l nodes=1:ppn=1
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/PQR_loma/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling PQR"
./PQR_LOMA
echo "Tofis is finitto!"
