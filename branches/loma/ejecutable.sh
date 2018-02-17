#!/bin/sh
###PBS -q default
##PBS -q newnodes
#PBS -l walltime=24:00:00
#PBS -l nodes=6:ppn=12
#PBS -N s20_sigma07_2
###PBS -N s40Big1_02
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/loma/
cd $WKD
echo "WKD is..."
echo $WKD
##cp hres20a.dat hre.dat
EXE=/home/toni/SVNlast/branches/loma/loma
echo "Calling mpirun for loMa"
mpirun -np 68 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/loma/loma
echo "LoMa is running... they are coming..."
