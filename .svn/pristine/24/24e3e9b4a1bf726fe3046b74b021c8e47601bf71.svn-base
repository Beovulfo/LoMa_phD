#!/bin/sh
###PBS -q default
#PBS -q prodA
#PBS -l walltime=99:00:00
#PBS -l nodes=11:ppn=12
#PBS -N s80aRe380_sig07
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/loma3/
cd $WKD
echo "WKD is..."
echo $WKD
##cp hres20a.dat hre.dat
EXE=/home/toni/SVNlast/branches/loma3/loma
echo "Calling mpirun for loMa3"
mpirun -np 132 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/loma3/loma
echo "LoMa is running... they are coming..."
