import sys

def writehre(pathIN,pathOUT,S,Lf, gamma,betha):
        Pr = 0.7
	sigma=0.7
	hrefile = open ('hre.dat','w')
        hrefile.write( str(sigma) +'   ' + str(Pr) + '   ' + str(S) +'   '+ str(Lf) + '   '+ str(gamma) + '  '+ str(betha) + '   1.0\n')  
        #hrefile.write('0.7   0.7   ' + str(S) +'   '+ str(Lf) + '   '+ str(gamma) + '  '+ str(betha) + '   1.0\n')  
	hrefile.write('1  64 \n')
        hrefile.write('/home/toni/swap\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()

fileFirst = 30
fileLast =  30


#path1 = '/data2/toni/turb/3DincRe160'
#path1 = '/share/drive/toni/Re160s80/case1/SS/s80aSS'
#path1 = '/share/drive/toni/reactml/q02D/Re1000'
#path0 = '/share/drive/toni/reactml/q02D'
#path0 = '/share/drive/toni/reactml/S40LF03g65_1D'
path0 = '/home/toni/fields'
#path0 = '/share/drive/toni/reactml/PantB/Re640'
#path0 = '/share/drive/toni/reactml/S04LF10g40'
#path0 = '/share/drive/toni/comet/S01LF10g6'
#path0 = '/share/drive/toni/reactml/S01LF10g6'
#path0 = '/share/drive/toni/reactml/S10LF032Dc'
#path0 = '/share/drive/toni/reactml/S10LF102Dc'
cname = 'PantB'
path1 = path0 + '/' + cname
#path1 = '/share/drive/toni/reactml/S10LF102D/Re1000'
#path1 = '/share/drive/toni/reactml/S10LF032D/Re1000'
#path1 = '/data2/toni/turb/s40SSmov'
#path1 = '/data2/toni/turb/s80my1101resx2'
#path1 = '/data2/toni/turb/3DTmax105b'
#path1 = '/data2/toni/turb/3DTmax125bis'
path1IN = path1 + '.'
path1OUT = path1 + '_'
#S = 0.0
#S = 1.0
S = 1.0
#S = 10.0
#Lf = 2.0
Lf = 1.0
#Lf = 1.0
#gamma = 0.0
#gamma = 3.73
gamma = 0.0
betha = 50
#fIC = 0
fIC = 0
from subprocess import call
#START WITH IC if wanted
if fIC == 1:
	fileIN = path0 + '/IC.dat'
	#fileIN = path0 + '/IC.dat'
        fileOUT = path1OUT + '%03d' % 0
        writehre(fileIN,fileOUT,S,Lf,gamma,betha)
        call("./tofis")
        raise SystemExit
	
	
for ifile in range(fileFirst,fileLast+1): 
	fileIN = path1IN + '%03d' % ifile
        fileOUT = path1OUT + '%03d' % ifile
        writehre(fileIN,fileOUT,S,Lf,gamma,betha)
        call("./tofis")
raise SystemExit
        
 
