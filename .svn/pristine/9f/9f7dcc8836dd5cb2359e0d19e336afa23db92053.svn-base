#!/bin/sh
##PBS -q prodA
#PBS -q mlx
#PBS -l walltime=24:00:00
#PBS -N VDMLBM1Ds40cont
#PBS -l nodes=1:ppn=12
##PBS -l nodes=3:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/lomacte3/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMacte"
mpirun -np 12 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomacte3/loma
#mpirun -np 36 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomacte3/loma
echo "LoMaCte finished they are coming..."
