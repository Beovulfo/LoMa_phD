def writehre(pathIN):
	hrefile = open ('hre.dat','w')
	hrefile.write(pathIN + '\n')
        hrefile.close()

import time
#
#path1 = '/share/drive/toni/Re160s80/case1/s80a.'
#path1 = '/share/drive/toni/coefT/Re160s40/case1/s40a.'
#path1 = '/share/drive/toni/Re160s40/case1/s40a.'
#path1 = '/share/drive/toni/Re160s40/case2/s40b.'
#path1 = '/share/drive/toni/Re160s40/case3/s40c.'
#path1 = '/share/drive/toni/Re160s10/case1/s10a.'
#path1 = '/share/drive/toni/Re160s10/case2/s10b.'
#path1 = '/share/drive/toni/Re160s20/case1/s20a.'
#path1 = '/share/drive/toni/Re160s20/case3/s20c.'
#path1 = '/share/drive/toni/Re160s80/case1/y2/s80a.'
#fileFirst = 18
#fileLast =  69

#path1 = '/share/drive/toni/Re160s40/case1/y2/s40a.'
#fileFirst = 14
#fileLast =  60

path1 = '/share/drive/toni/VDML/s10/01/s10a.'
#path1 = '/share/drive/toni/Re160s80/case2/y2/s80b.'
fileFirst = 13
fileLast =  21
#path1 = '/share/drive/toni/Re160s80/case1/SS/s80aSS.'
#path1 = '/share/drive/toni/coefT/Re160s80/case1/s80aDx03.'
#path1 = '/share/drive/toni/Re160s80/case1/SS/s80aSS.'
#path1 = '/share/drive/toni/coefT/Re160s80/case1/s80aDx03.'
l_files = [path1 + '%03d' %ifile for ifile in range(fileFirst,fileLast)]
l_files_tot = l_files
#
#path1 = '/share/drive/toni/Re160s10/case2/s10b2.'
#path1 = '/share/drive/toni/Re160s80/case3/s80c2.'
#path1 = '/share/drive/toni/Re160s40/case1/s40a.'
#path1 = '/share/drive/toni/Re160s40/case3/s40c2.'
#path1 = '/share/drive/toni/Re160s80/case2/s80b2.'
#fileFirst = 101
#fileLast =  114
#l_files = [path1 + '%03d' %ifile for ifile in range(fileFirst,fileLast)]
#l_files_tot += l_files
##
##
#path1 = '/share/drive/toni/Re160s80/case2/s80b2.'
#fileFirst = 81
#fileLast =  88
#l_files = [path1 + '%03d' %ifile for ifile in range(fileFirst,fileLast)]
#l_files_tot += l_files
##
##
#
#path1 = '/share/drive/toni/Re160s80/case3/s80c.'
#fileFirst = 1
#fileLast =  8
#l_files = [path1 + '%03d' %ifile for ifile in range(fileFirst,fileLast)]
#l_files_tot += l_files
#
#path1 = '/share/drive/toni/Re160s80/case3/s80c2.'
#fileFirst = 81
#fileLast =  88
#l_files = [path1 + '%03d' %ifile for ifile in range(fileFirst,fileLast)]
#l_files_tot += l_files
#
#import subprocess
from subprocess import call
#from subprocess import call
count = 0
for ifile in l_files_tot: 
        writehre(ifile)
        print "sending new job for =%s"%(ifile)
	print "Sending ejecutable.sh"
	call(["qsub", "ejecutable.sh"])
        #call(["qsub", "ejecutableT.sh"])
	time.sleep(315)
	#subprocess.check_call(['qsub ejecutable.sh'],shell=True)
        count+=1
