#!/bin/sh
#----------Start job description--------------------------------
#---arch = power
#@ arch = ppc64 
#@ initialdir =/gpfs/project/uc3m43/PROJECT/lomacte/ 
#@ total_tasks = 132
#@ error = /gpfs/project/uc3m43/SCRATCH/errRe160s40e_01-%j.log 
#@ output = /gpfs/project/uc3m43/SCRATCH/Re160s40my1101e_01-%j.log 
#@ wall_clock_limit = 72:00:00
#-----------End job description----------------------------------
#  
export GFORTRAN_CONVERT_UNIT=little_endian
echo "Calling srun for loMa1101"
srun -n 132 /gpfs/project/uc3m43/PROJECT/lomacte/loma1101
echo "LoMa1101 is finished!"
