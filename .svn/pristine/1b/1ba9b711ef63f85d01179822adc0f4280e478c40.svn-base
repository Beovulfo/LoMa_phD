#!/bin/sh
##PBS -q newnodes
##PBS -q vipnew
#PBS -l walltime=23:59:59
#PBS -N bm_nonreactml4
#PBS -l nodes=6:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/nonreactml4/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for NoNReactML-loMaHZ Single precision"
#mpiexec -hostfile $PBS_NODEFILE -n 18 /home/toni/SVNlast/branches/lomaHZ/lomahz
mpirun -np 72 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/nonreactml4/nonreactml
echo "NoN-ReactML-LoMaHZ is running... they are coming..."
