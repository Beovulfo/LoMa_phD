def writehre(pathIN,pathOUT):
	hrefile = open ('hre.dat','w')
	hrefile.write('1   128 \n')
        hrefile.write('/home/toni/swap\n')
	hrefile.write(pathIN + '\n')
	hrefile.write(pathOUT + '\n')
        hrefile.close()

vpathIN = ['/home/toni/fields/VDML_JFM/s20c.006']
#vpathIN = ['/home/toni/fields/s80a.119']
vpathOUT = ['/home/toni/fields/VDML_JFM/s20c_006']
#vpathIN = ['/home/toni/fields/Big1.035']
#vpathOUT = ['/home/toni/fields/Big1_035']
from subprocess import call
for ifile in range(len(vpathIN)): 
	fileIN = vpathIN[ifile]
        fileOUT = vpathOUT[ifile]
        writehre(fileIN,fileOUT)
        call("./tofis")
        
         

 
