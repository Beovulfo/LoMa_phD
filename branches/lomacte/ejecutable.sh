#!/bin/sh
##PBS -q vipnew
#PBS -l walltime=24:00:00
###PBS -N Re160s20c_02
#PBS -N s80lomacte
#PBS -l nodes=6:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/lomacte/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMacte"
mpirun -np 72 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomacte/loma
echo "LoMaCte finished they are coming..."
