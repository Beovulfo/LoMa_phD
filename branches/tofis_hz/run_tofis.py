import os
import time
def writehre(pathIN,pathOUT):
	hrefile = open ('hre.dat','w')
	hrefile.write('1   16 \n')
        hrefile.write('/share/drive/toni/swap\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()

fileFirst = 1
fileLast =  20

#path1 = '/data2/toni/turb/3DincRe160'
path1 = '/share/drive/toni/Re160s80/case1/SS/s80aSS'
#path1 = '/data2/toni/turb/s80my1101resx2'
#path1 = '/data2/toni/turb/3DTmax105b'
#path1 = '/data2/toni/turb/3DTmax125bis'
path1IN = path1 + '.'
#path1OUT = '/share/drive/toni/s80aSS_'
path1OUT = path1 + '_'
#Wait for the file
flag = False
counter = 0
from subprocess import call
for ifile in range(fileFirst,fileLast+1): 
	fileIN = path1IN + '%03d' % ifile
	#fileIN = path1IN + '%03d' % ifile
 	while not flag:
		if ifile>fileFirst:
			time.sleep(120)
		flag = os.path.isfile(fileIN)
		counter = counter + 1	
		if counter > 10000:
			exit()
        fileOUT = path1OUT + '%03d' % ifile
        writehre(fileIN,fileOUT)
        call("./tofis")
        
         

 
