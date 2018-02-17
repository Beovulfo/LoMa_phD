#!/bin/sh
##PBS -q newnodes
##PBS -q oldnodes
#PBS -l walltime=24:00:00
#PBS -N PantBLf03
#PBS -l nodes=11:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/lomaHZ/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMaHZ"
#mpiexec -hostfile $PBS_NODEFILE -n 18 /home/toni/SVNlast/branches/lomaHZ/lomahz
mpirun -np 132 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZ/lomahz
echo "LoMaHZ is running... they are coming..."
