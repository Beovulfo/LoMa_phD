#!/bin/sh
#PBS -lwalltime=240:30:00
#PBS -lnodes=18
#PBS -N mixlayernoise
#PBS -q oldnodes
#PBS -k oe
#PBS -j oe
# Parameters:
#  
# HAY QUE CAMBIAR EL WORKING DIRECTORY CADA VEZ!!!!
WKD=/home/toni/SVN2/LowMach/branches/mixing/
cd $WKD
echo "WKD is..."
echo $WKD
echo "Calling pc_run module with user .conf file"
mpirun -np 18 /home/toni/SVN2/LowMach/branches/mixing/liso
#$PENCIL_HOME/bin/pc_run run run -f $PENCIL_HOME/config/hosts/toni/host-icaro.uc3m.es.conf 

echo "Liso running"
