from writefiles import create_mesh1D 

#Define mesh for writing
Ly=172/2.
my = 1053
gamma = 0.5
fmesh = '../data/mesh.txt'

#Create mesh..
print "writing mesh at %s" % fmesh

create_mesh1D(Ly,my,'tanh',gamma,fmesh)

"""
def writehre(pathIN,pathOUT):
	hrefile = open ('hre.dat','w')
	hrefile.write('1   16 \n')
        hrefile.write('/data2/toni/swap\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()

fileFirst = 41
fileLast = 41

#path1 = '/data2/toni/turb/3DincRe160'
path1 = '/data2/toni/turb/3Dinc01'
path1IN = path1 + '.'
path1OUT = path1 + '_'
from subprocess import call
for ifile in range(fileFirst,fileLast+1): 
	fileIN = path1IN + '%03d' % ifile
        fileOUT = path1OUT + '%03d' % ifile
        writehre(fileIN,fileOUT)
        call("./tofis")
        
         

""" 
