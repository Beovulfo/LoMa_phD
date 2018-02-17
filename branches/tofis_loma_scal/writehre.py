def writehre(pathIN,pathOUT):
	hrefile = open ('hre.dat','w')
	hrefile.write('1   16 \n')
        hrefile.write('/share/drive/toni/swap\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()

fileFirst = 7
fileLast =  9

#path1 = '/data2/toni/turb/3DincRe160'
#path1 = '/share/drive/toni/Re160s40/case3/s40c2'
#path1 = '/share/drive/toni/Re160s10/case3/s10c'
#path1 = '/share/drive/toni/Re160s20/case3/s20c2'
#path1 = '/share/drive/toni/Re160s10/case3/s10c'
path1 = '/share/drive/toni/Re160s10/case3/s10c'
#path1 = '/data2/toni/turb/s40SSmov'
#path1 = '/data2/toni/turb/s80my1101resx2'
#path1 = '/data2/toni/turb/3DTmax105b'
#path1 = '/data2/toni/turb/3DTmax125bis'
path1IN = path1 + '.'
path1OUT = path1 + '_'
from subprocess import call
for ifile in range(fileFirst,fileLast+1): 
	fileIN = path1IN + '%03d' % ifile
        fileOUT = path1OUT + '%03d' % ifile
        writehre(fileIN,fileOUT)
        call("./tofis")
        
         

 
