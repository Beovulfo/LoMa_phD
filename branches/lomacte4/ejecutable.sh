#!/bin/sh
##PBS -q prodA
#PBS -l walltime=4:00:00
#PBS -N loma4s801D
#PBS -l nodes=6:ppn=12
##PBS -l nodes=3:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/lomacte4/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMacte"
mpirun -np 72 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomacte4/loma
#mpirun -np 36 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomacte4/loma
echo "LoMaCte finished they are coming..."
