#!/bin/sh
##PBS -q newnodes
##PBS -q newnodes
#PBS -l walltime=00:45:00
#PBS -N Interp
#PBS -l nodes=1:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/interp_reactml/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for Interp-ReactML"
#mpiexec -hostfile $PBS_NODEFILE -n 18 /home/toni/SVNlast/branches/lomaHZ/lomahz
mpirun -np 12 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/interp_reactml/interp
