#!/bin/sh
#PBS -q default
#PBS -l walltime=12:00:00
#PBS -N PQRHZB
#PBS -l nodes=1:ppn=1
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/PQR_lomaHZB/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling PQR_lomaHZB"
./tofis
echo "PQR Ã finitto!"
