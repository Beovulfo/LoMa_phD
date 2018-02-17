#!/usr/local/bin/python2.7
import sys
from readfiles import *
from fortran2vtk import *

        
#Python script to convert planes yx to VTK format --> PARAVIEW
#BaseName = '/share/drive/toni/reactml/q02D/Re1000'
#varpln = 'o3yz'
#nmin=1;nmax=10;stride = 1;
#binary = 1
#Example of use:
#python2.7 convert2vtk.py /share/drive/toni/reactml/q02D/Re1000 o3yx 1 10 1
#Read Basename
BaseName = str(sys.argv[1])
#print BaseName
#Read variable and plane
varpln = str(sys.argv[2])
pln = varpln[-2:]
#First, last and striide files
nmin=int(sys.argv[3]);nmax=int(sys.argv[4]);stride = int(sys.argv[5]);
if len(sys.argv)>6:
	binary = sys.argv[6]
else:
	binary = 1 #Default using binary

for i in range(nmin,nmax+1,stride):
    	filename = BaseName+'_'+'%3.3i'%i+'.'+varpln
	print 'Opening: ' + filename
    	[field,y,x,z] = readfieldij(filename);
	print 'Fis-Field read.'
        print 'Field shape:',np.shape(field)
	pln2vtk(x,y,z,field,varpln,filename,binary)
    	print 'VTK created at \n'
    	print filename+'.vtk'

