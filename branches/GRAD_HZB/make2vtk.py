
#!/usr/local/bin/python2.7
import sys
import re
import os.path
from readfiles import *
from fortran2vtk import *
from subprocess import call
import subprocess
import time
runlocal=True#Choose if run local or qsub

def writehre(pathIN,pathOUT,S,Lf, gamma,betha):
        Pr = 0.7
	sigma=0.7
	hrefile = open ('hre.dat','w')
        hrefile.write( str(sigma) +'   ' + str(Pr) + '   ' + str(S) +'   '+ str(Lf) + '   '+ str(gamma) + '  '+ str(betha) + '   1.0\n')  
        #hrefile.write('0.7   0.7   ' + str(S) +'   '+ str(Lf) + '   '+ str(gamma) + '  '+ str(betha) + '   1.0\n')  
	hrefile.write('1  32 \n')
        hrefile.write('/home/toni/swapCom\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()
#INPUT DATA
#S=1.0; Lf=1.0; gamma=6.0;betha=100.0;
#S=4.0; Lf=1.0; gamma=4.0;betha=100.0;
#H2
#S=34.5; Lf=0.3; gamma=4.0;betha=300.0;
#RML/BM/01 (Lf=1)
#S=4.0; Lf=1.0; gamma=4.0;betha=100.0;
#RML/BM/02 (Lf=2)
#S=4.0; Lf=2.0; gamma=4.0;betha=100.0;
#RML/BM/03 (Lf=0.3)
#S=4.0; Lf=0.3; gamma=4.0;betha=100.0;
#RML/BM/setA/LF1 
#S=1.0; Lf=1; gamma=4.0;betha=300.0;
#HomeDir = "/home/toni/fields/RML/BM/setA/LF1/"

#RML/BM/setA/LF03 
#S=1.0; Lf=0.3; gamma=4.0;betha=300.0;
#HomeDir = "/home/toni/fields/RML/BM/setA/LF03/"

#RML/BM/setA/LF2 
#S=1.0; Lf=2; gamma=4.0;betha=300.0;
#HomeDir = "/home/toni/fields/RML/BM/setA/LF2/"

#RML/LF03S005/ 
#S=0.5; Lf=0.3; gamma=4.0;betha=100000;
#HomeDir = "/home/toni/fields/RML/LF03S005/"

S=1.0; Lf=0.3; gamma=4.0;betha=100000;
HomeDir = "/home/toni/fields/RML/LF03S1/"#037
#Python script to convert planes yx to VTK format --> PARAVIEW
#BaseName = '/share/drive/toni/reactml/q02D/Re1000'
#varpln = 'o3yz'
#nmin=1;nmax=10;stride = 1;
#binary = 1
#Example of use:
#python2.7 make2vtk.py /share/drive/toni/Re160s10/case3/ s10c Tfyz 8 8 0
#Read Basename
initPath= str(sys.argv[1])
jobname = str(sys.argv[2])
BaseName = initPath+jobname 
#HomeDir = "/home/toni/fields/RML/BM/2D/LF1/"
#HomeDir = "/home/toni/fields/RML/BM/01/"
#HomeDir = "/home/toni/fields/RML/BM/02/"
#HomeDir = "/home/toni/fields/RML/BM/setA/LF1/"
#HomeDir = "/home/toni/fields/RML/BM/H2/"
#HomeDir  "/home/toni/fields/RML/setA/LF1/D/"
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
	filecheck = fileOUTPUT + '.qqyx'
        #Check if tofis already done,filecheck is just in case we put
 	#the varpln wrong!
	if os.path.isfile(filename2) or os.path.isfile(filecheck):
		print "Field already transformed to FIS :)"
	else:
		if os.path.isfile(fileINPUT):
			print "Image already found at %s"%fileINPUT
		else:
			print "Bringing Image home!"
			call(["rsync", "-avz",fileIMA,fileINPUT])
			print "File %s has been copied as %s" %(fileIMA,fileINPUT)
       		writehre(fileINPUT,fileOUTPUT,S,Lf,gamma,betha)
		if runlocal == True:
        		call("./GRADHZB")
        		#call("./tofis851")
                #Instead we will send this to queue!
		else:
        		#call("qsub tofis.sh")
			qsub_call = "qsub %s"
			call(qsub_call % "GRAD.sh", shell=True)
		#Now wait until is finished
			fqstat=True
			word = "GRADHZB";#Name of job
			while (fqstat==True):
				proc = subprocess.Popen('qstat', stdout=subprocess.PIPE)
				output = proc.stdout.read()
				if word in output:
					print "Still %s is queued... Waiting..."%(word)
					time.sleep(60)
					fqstat=True
				else:
					fqstat=False
        call(["rm",fileINPUT])
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
