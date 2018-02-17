def writehre(pathIN,pathOUT):
	hrefile = open ('hre.dat','w')
	hrefile.write('1   16 \n')
        hrefile.write('/data2/toni/swap\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()

fileFirst = 14
fileLast = 45

path1 = '/data2/toni/turb/3DincRe160'
path1IN = path1 + '.'
path1OUT = path1 + '_'
from subprocess import call
for ifile in range(fileFirst+1,fileLast+1): 
	fileIN = path1IN + '%03d' % ifile
        fileOUT = path1OUT + '%03d' % ifile
        writehre(fileIN,fileOUT)
        call("./tofis")
        
         

 
