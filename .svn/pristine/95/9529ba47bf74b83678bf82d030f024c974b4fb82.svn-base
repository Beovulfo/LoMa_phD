#!/bin/sh
##PBS -q newnodes
##PBS -q newnodes
#PBS -l walltime=03:59:59
##PBS -N bm_reactml_higA
#PBS -N DEBUG_REACTML
#PBS -l nodes=3:ppn=12
#PBS -k oe
#PBS -j oe
# Parameters:
#  
##WKD=/home/toni/SVNlast/branches/newloma/
WKD=/home/toni/SVNlast/branches/reactml/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling mpirun for ReactML-loMaHZ"
#mpiexec -hostfile $PBS_NODEFILE -n 18 /home/toni/SVNlast/branches/lomaHZ/lomahz
mpirun -np 36 -hostfile $PBS_NODEFILE /home/toni/SVNlast/branches/reactml/reactml
echo "ReactML-LoMaHZ is running... they are coming..."
