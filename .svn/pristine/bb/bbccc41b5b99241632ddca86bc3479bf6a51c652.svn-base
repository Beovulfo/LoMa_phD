#!/bin/sh
##PBS -q mlx
#PBS -q qib
#PBS -l walltime=40:59:59
##PBS -l walltime=26:00:00
#PBS -N  VORTEX
##PBS -l nodes=3:ppn=12
#PBS -l nodes=6:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/vortexHZB/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for vortexHZB"
#mpiexec -np  36 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/vortexHZB/vortexhzb
mpiexec -np  72 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/vortexHZB/vortexhzb
#mpiexec -np  20 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/vortexHZB/vortexhzb
#mpiexec -np  36 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/vortexHZB/vortexhzb
echo "VorteXHZB is running... they are coming..."
