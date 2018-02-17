#!/bin/sh
#PBS -q qib
##PBS -l walltime=23:59:59
#PBS -l walltime=47:59:59
##PBS -l walltime=26:00:00
#PBS -N  LF10S4Re8000_02
###PBS -N  LF03S15Re10000_01
##PBS -l nodes=compute-0-2:ppn=12+compute-0-4:ppn=12+compute-0-9:ppn=12+compute-0-10:ppn=12+compute-0-11:ppn=12+compute-0-12:ppn=12
##PBS -l nodes=compute-0-1:ppn=12+compute-0-2:ppn=12+compute-0-4:ppn=12
#PBS -l nodes=6:ppn=12
##PBS -l nodes=11:ppn=12
##PBS -l nodes=6:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
WKD=/home/toni/SVNlast/branches/lomaHZB/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for loMaHZB"
#mpiexec -np 68 /home/toni/SVNlast/branches/lomaHZB/lomahzb
#mpiexec -np 60 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZB/lomahzb
#mpiexec -np 36 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZB/lomahzb
#mpiexec -np 108 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZB/lomahzb
mpiexec -np  72 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZB/lomahzb
#mpiexec -np  132 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZB/lomahzb
#mpiexec -np 72 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZB/lomahzb
#mpiexec -np 132 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/lomaHZB/lomahzb
echo "LoMaHZ is running... they are coming..."
