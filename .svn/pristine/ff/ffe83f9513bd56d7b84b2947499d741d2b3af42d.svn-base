#!/bin/sh
##PBS -q newnodes
##PBS -q newnodes
#PBS -l walltime=23:59:59
#PBS -N Interp
#PBS -l nodes=1:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/interp_loma/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for Interp-LOMA"
#mpiexec -hostfile $PBS_NODEFILE -n 18 /home/toni/SVNlast/branches/lomaHZ/lomahz
mpirun -np 12 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/interp_loma/interp
