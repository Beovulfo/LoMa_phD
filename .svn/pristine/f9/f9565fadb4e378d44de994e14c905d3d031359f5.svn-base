def writehre(path,iout,pout,pin,pstats,hrename):
    #hrename='hre.dat'
    Re=160.0;a0=0.020466;b0=0.04096;a=0.0
    nstep=1601;nimag=400;nhist=100;
    CFL = 0.25; ncfl=10;
    stats=1;nacum=10;
    #path='/share/drive/toni/coef/Re160s20/'
    #fileout=path+'case1/s20a'
    fileout = path+pout
    filein  = path+pin
    filestats= path+pstats
    #filein=path+'case1/IC.dat'
    #filestats=path+'case1/s20a_01'
    hrefile = open (hrename,'w')     
    #First line with Re+alpha+beta+a
    hrefile.write("CC Re"+'\t'+"alpha0"+'\t'+"beta0"+'\t'+"a0"+'\n')    
    hrefile.write(str(Re)+'\t'+str(a0)+'\t'+str(b0)+'\t'+str(a)+'\n')    
    #Second line with 
    #nstep    nimag     nhist 
    hrefile.write("CC nsteps"+'\t'+"nimag"+'\t'+"nhist"+ '\n')
    hrefile.write(str(nstep)+'\t'+str(nimag)+'\t'+str(nhist)+ '\n')
    # CFL steps cfl change        
    hrefile.write("CC CFL"+'\t'+"ncfl"+'\n')
    hrefile.write(str(CFL)+'\t'+str(ncfl)+'\n')
    # iout stats nacum
    hrefile.write("CC iout"+'\t'+"stats"+'\t'+"nacum"+ '\n')
    hrefile.write(str(iout)+'\t'+str(stats)+'\t'+str(nacum)+ '\n')
    hrefile.write("CC fileout"+'\n')                                                  
    hrefile.write(fileout+'\n')                                                  
    hrefile.write("CC filein" + '\n')                                                               
    hrefile.write(filein + '\n')                                                               
    hrefile.write("CC filestats" + '\n')   
    hrefile.write(filestats + '\n')   
    hrefile.close() 


path='/share/drive/toni/coefT/Re160s20/'
pout=[  'case1/s20a','case1/s20a',\
	'case2/s20b','case2/s20b',\
	'case3/s20c','case3/s20c']
#pout='case1/s20a'
pin=[   'case1/IC.dat','case1/s20a.004',\
	'case2/IC.dat','case2/s20b.004',\
	'case3/IC.dat','case3/s20c.004'];
#pin=['case1/IC.dat','case1/s20a.004',];
pstats=['case1/s20a_01','case1/s20a_02',\
	'case2/s20b_01','case2/s20b_02',\
	'case3/s20c_01','case3/s20c_02' ];
#pstats=['case1/s20a_01','case1/s20a_02'];
hrenames=['hres20a.dat','hres20acont.dat',\
	  'hres20b.dat','hres20bcont.dat',\
	  'hres20c.dat','hres20ccont.dat'];
#hrenames=['hres20a.dat','hres20acont.dat'];
iout=[  '1','5',\
	'1','5',\
	'1','5']

for i in range(len(hrenames)): 
	writehre(path,iout[i],pout[i],pin[i],pstats[i],hrenames[i])
        
         

 
