
#!/usr/local/bin/python2.7
import sys
import re
import os.path
from readfiles import *
from fortran2vtk import *
from subprocess import call

def writehre(pathIN,pathOUT):
	hrefile = open ('hre.dat','w')
	hrefile.write('1   128 \n')
        hrefile.write('/home/toni/swap\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()


        
#Python script to convert planes yx to VTK format --> PARAVIEW
#BaseName = '/share/drive/toni/reactml/q02D/Re1000'
#varpln = 'o3yz'
#nmin=1;nmax=10;stride = 1;
#binary = 1
#Example of use:
#python2.7 make2vtk.py /share/drive/toni/Re160s20/case3/ s20c Tfyx 8 8 1 0
#Read Basename
initPath= str(sys.argv[1])
jobname = str(sys.argv[2])
BaseName = initPath+jobname 
#HomeDir = "/share/drive/toni/VDML/s80/test/04/"
#HomeDir = "/share/drive/toni/VDML/s80/test/02/"
#HomeDir = "/share/drive/toni/VDML/s80/test/00/"
#HomeDir = "/home/toni/fields/"
HomeDir = "/home/toni/fields/VDML_JFM/"
#Take the mean or not
#print BaseName
#Read variable and plane
varpln = str(sys.argv[3])
pln = varpln[-2:]
#First, last and striide files
nmin=int(sys.argv[4]);
nmax=int(sys.argv[5]);
stride = int(sys.argv[6]);
pert = int(sys.argv[7]);
#if len(sys.argv)>7:
#	binary = sys.argv[7]
#else:
binary = 1 #Default using binary

#Lets print
print "my HomeDir is: %s" %HomeDir
print "my input Dir is: %s" %BaseName
print "my varpln is: %s" %varpln
print "my pln is: %s" %pln
print "starting from: %s" %nmin
print "ending at: %s" %nmax

for i in range(nmin,nmax+1,stride):
# 
    	fileIMA = BaseName+'.'+'%3.3i'%i
	fileINPUT = HomeDir + jobname + '.'+'%3.3i'%i
	fileOUTPUT = HomeDir + jobname + '_'+'%3.3i'%i
    	filename = BaseName+'_'+'%3.3i'%i+'.'+varpln
    	filename2 = fileOUTPUT +'.'+varpln
	filecheck = fileOUTPUT + '.upxz'
	newctes3D = True
        #Check if tofis already done,filecheck is just in case we put
 	#the varpln wrong!
	#if os.path.isfile(filename2) and newctes3D==False:
	if os.path.isfile(filename2) or os.path.isfile(filecheck):
		print "Field already transformed to FIS :)"
	else:
		if os.path.isfile(fileINPUT):
			print "Image already found at %s"%fileINPUT
		else:
			print "Bringing Image home!"
			call(["rsync", "-avz",fileIMA,fileINPUT])
			print "File %s has been copied as %s" %(fileIMA,fileINPUT)
       		writehre(fileINPUT,fileOUTPUT)
        	call("./tofis")
        #call(["rm",fileINPUT])
        #Now convert 2 vtk
        	
    	[field,y,x,z] = readfieldij(filename2,pert);
	if pert==1:
		print "Taking the mean"
		filename2=filename2+'p'
 	print 'Reading %s'%(filename2)
        print 'Field shape:',np.shape(field)
	pln2vtk(x,y,z,field,varpln,filename2,binary)
    	print 'VTK created at \n'
    	print filename2+'.vtk'
#
#
#vpathIN = ['/home/toni/fields/s80a.119']
#vpathOUT = ['/home/toni/fields/s80a_119']
##vpathIN = ['/home/toni/fields/Big1.035']
##vpathOUT = ['/home/toni/fields/Big1_035']
##call(["rm", "kk2"])
##for ifile in range(len(vpathIN)): 
##	fileIN = vpathIN[ifile]
##        fileOUT = vpathOUT[ifile]
##        writehre(fileIN,fileOUT)
##        call("./tofis")
#        
        

 
