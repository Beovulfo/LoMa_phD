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
#python2.7 create_vtk_vector.py /share/drive/toni/reactml/q02D/Re1000 o3yx 1 10 1
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

field = [1 for j in range(3)]
if varpln[0:2]=='up':
	for i in range(nmin,nmax+1,stride):
		pl = varpln[-2:]
 		print 'pl = ' + str(pl)
		filename = BaseName+'_'+'%3.3i'%i+'.'
		print 'Opening: ' + filename
    		[field[0],y,x,z] = readfieldij(filename+'up'+pl);
    		[field[1],y,x,z] = readfieldij(filename+'vp'+pl);
    		[field[2],y,x,z] = readfieldij(filename+'wp'+pl);
        	print 'Field shape:',np.shape(field)
		create_vtk_vector(x,y,z,field,'vel'+pl,filename+'vel'+pl,binary)
    		print 'VTK created for VELOCITY vector at'
    		print filename+'vel'+pl+'.vtk'
elif varpln[0:2]=='o1':
	for i in range(nmin,nmax+1,stride):
		pl = varpln[-2:]
 		print 'pl = ' + str(pl)
		filename = BaseName+'_'+'%3.3i'%i+'.'
		print 'Opening: ' + filename
    		[field[0],y,x,z] = readfieldij(filename+'o1'+pl);
    		[field[1],y,x,z] = readfieldij(filename+'o2'+pl);
    		[field[2],y,x,z] = readfieldij(filename+'o3'+pl);
        	print 'Field shape:',np.shape(field)
		create_vtk_vector(x,y,z,field,'vor'+pl,filename+'vor'+pl,binary)
    		print 'VTK created for VORTICITY vector at'
    		print filename+'vor'+pl+'.vtk'
else:
	print 'WHaT???'

