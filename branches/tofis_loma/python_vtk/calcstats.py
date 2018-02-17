import pylab
from scipy import *
from scipy.integrate import cumtrapz
from scipy.integrate import simps
from scipy.integrate import romb
from math  import *
#from pylab import *
from numpy import *
import numpy as np
from fortran2vtk import *
import matplotlib.pyplot as plt

def create_filenames(varname,nsteps,rksteps=3):
	""" filenames  = create_filenames(varname,nsteps,*rksteps) 
	    list       = create_filenames(str, int)
	    returns a list of filenames live varname.+001001,001002,...,nsteps/003
	Ex:
	print create_filenames("v00.", 2,0)
        >>> ['v00.001','v00.002']	
	"""
	filenames=[]
	for istep in range(nsteps):
		if rksteps != 0:
    			for rkstep in range(rksteps):
        			filename=  varname +'%03d' % (istep+1) + '%03d' % (rkstep+1)
				filenames.append(filename)
		else:
        		filename=  varname +'%03d' % (istep+1)
			filenames.append(filename)
	return filenames

def create_stanames(varname,nstepi,nstepf):
	""" filenames  = create_stanames(varname,nstepi,nstepf) 
	    list       = create_stanames(str, int,int)
	    returns a list of filenames live varname_045.sta'
	Ex:
	print create_stanames("v00.",20,21)
        >>> ['v00_020.sta','v00_021.sta']	
	"""
	filenames=[]
	for istep in range(nstepi,nstepf):
       		filename=  varname +'_'+'%03d' % (istep+1)+'.sta'
		filenames.append(filename)
	return filenames

#----------------------------------------------------
def create_spenames(varname,run,nstepi,nstepf):
	""" filenames = create_stanames(varname,run,nstepi,nstepf) 
	    list      = create_stanames(str,int, int,int)
	    returns a list of filenames like varname_01_045.spe'
	Ex:
	print create_spenames("v00.",20,21)
        >>> ['v00_01_020.spe','v00_01_021.spe']	
	"""
	filenames=[]
	for istep in range(nstepi,nstepf):
       		filename=  varname +'_' + \
		'%02d' % run +      '_' + \
		'%03d' % (istep+1)+'.spe'
		filenames.append(filename)
	return filenames

#-----------------------------------------------------------------------------------------------------#
def readASCII(filename):
	""" open a text file with "ncolumns"
	 and reads data, returns a list of lists
	 Important: data must be numbers only and same columns for all lines! 
	[list] = readASCII([str])
	"""
	import numpy as np
	f1 = open(filename,'r')
	file = f1.readlines()
	f1.close()
	#Take one line for info
	ncolumns= len(file[0].split()) 
	nlines = len(file)
	print "Reading file of %s lines and %s columns" %(nlines,ncolumns)
	x = [[] for i in range(ncolumns)]	
	for line in file:
		row = line.split()
		for col in range(len(row)):
			x[col].append(float(row[col].replace('D','E')))
	#return x
	return np.array(x)
			
			
def makeplot(x,col1=0,col2=1):
	import pylab as pl
	fig = pl.figure()
	pl.plot(x[col1,:],x[col2,:],'go')
	#fig.show()
	pl.show()
	print 'Continue computation'






#3	mom2 = rum[0]*um[0]*Ly/2.
#3	integ1 = trapz(rum*um,y)
#3	return (mom1 + mom2 - integ1)/(r0*DU**2)


def calcdmcomp(y,um,rum,rhom):
	""" Pantano definition
	"""
	from numpy import trapz
	DU = um.max()-um.min()
	if np.max(rhom)==np.mean(rhom): #density ratio = 1.0
		r0 = 1.0
		favum = um
		rhom = um/um #ones
		result = 1./(r0*DU**2.0)*trapz(rhom*(0.5*DU-favum)*(0.5*DU+favum),y)
	else: #constant density case
		r0 = (rhom[0]+rhom[-1])/2.0
		favum = rum/rhom
		result = 1./(r0*DU**2.0)*trapz(rhom*(0.5*DU-favum)*(0.5*DU+favum),y)
	return result

def calcdmcomp2(y,um,rum,rhom):
	""" Higuera definition
	"""
	from numpy import trapz
	DU = um.max()-um.min()
	r0 = 1.0
	favum = um
	result = trapz(1.0-rhom*um**2,y)
	#result = 1./(r0*DU**2.0)*trapz(rhom*(0.5*DU-favum)*(0.5*DU+favum),y)
	return result



def calcdw(y,um):
	DU = um.max()-um.min()
	m = (abs(der1(y,um))).max()
	if m ==0:
		return 0.0
	else:
		return DU/m 

def calcdwcomp(y,um,rum,rhom):
	DU = um.max()-um.min()
	if np.max(rhom)==np.mean(rhom): #density ratio = 1.0
		m = (abs(der1(y,um))).max()
	else:
		m = (abs(der1(y,rum/rhom))).max()
	if m ==0:
		return 0.0
	else:
		return DU/m 

def calcdmcomp12(y,um,rum,rhom):
	""" Pantano definition
	"""
	from numpy import trapz
	DU = um.max()-um.min()
	if np.max(rhom)==np.mean(rhom): #density ratio = 1.0
		r0 = 1.0
		favum = um
		rhom = um/um #ones
		dm1inf = 1./(r0*DU**2.0)*trapz(rhom*(0.5*DU-favum)*(0.5*DU+favum),y)/2.0
		dm2inf = 1./(r0*DU**2.0)*trapz(rhom*(0.5*DU-favum)*(0.5*DU+favum),y)/2.0
	else: #Variable density definition
		pos1=where(rhom>=1.0) #pos y where rho=1.0
		pos2=where(rhom<1.0) #pos y where rho=1.0
		r0 = (rhom[0]+rhom[-1])/2.0
		favum = rum/rhom
		dm1inf = 1./(r0*DU**2.0)*trapz(rhom[pos1]*(0.5*DU-favum[pos1])*(0.5*DU+favum[pos1]),y[pos1])
		dm2inf = 1./(r0*DU**2.0)*trapz(rhom[pos2]*(0.5*DU-favum[pos2])*(0.5*DU+favum[pos2]),y[pos2])
	return dm1inf,dm2inf




def der1(y,u):
	"""computes first derivative using Central Finite Differences
        """
	result = array(u)
        for j in range(len(y)):
 		if j==0:
			result[j]=(u[j+1]-u[j])/(y[j+1]-y[j])
		elif j==(len(y)-1):
			result[j]=(u[j]-u[j-1])/(y[j]-y[j-1])
		else:
			result[j]=(u[j+1]-u[j-1])/(y[j+1]-y[j-1])
	
	return result

def der1_3D(y,u,dim):
	"""computes first derivative using Central Finite Differences
        """
	result = u
        for j in range(len(y)):
 		if j==0:
			result[j]=(u[j+1]-u[j])/(y[j+1]-y[j])
		elif j==(len(y)-1):
			result[j]=(u[j]-u[j-1])/(y[j]-y[j-1])
		else:
			result[j]=(u[j+1]-u[j-1])/(y[j+1]-y[j-1])
	
	return result




def der2(y,u):
	"""Computes second derivative using central differences
	"""
	result = array(u)
	for j in range(len(y)):
		if j==0:
			result[j] = (u[j+2]-2*u[j+1]+u[j])/((y[j+2]-y[j])/2.)**2
		elif j==(len(y)-1):
			result[j] = (u[j-2]-2*u[j-1]+u[j])/((y[j-2]-y[j])/2.)**2
		else:
			result[j]=(u[j+1]-2*u[j]+u[j-1])/((y[j+1]-y[j-1])/2.)**2
        return result




#===============================================================#

def readspe(filename,nvar=12):
	""" This function reads the spe file generated by LOMA
#(yjsp,kx,kz,SPEC) = readspe(filename,nvar)
	(y(jsp),kx,kz,[nz/2,nspec,nvar,nx/2] = readspe (str,int)
	"""
	import scipy as sc
	import numpy as np
	import pylab
	f = open(filename,'rb')
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
                ('nx','uint32'),('ny','uint32'), \
                ('nz','uint32'),('np','uint32'), \
		('nacum','uint32'), \
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)
	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for RECORD1...!"
	#else:
    	#	print "Fiel read correctly :)"
        #Save "my"
	my = RECORD1['ny']
	nx = RECORD1['nx'][0]/2
	nz = (RECORD1['nz'][0]+1)/2
	kx = range(nx)*RECORD1['alp']
	kz = range(nz)*RECORD1['bet']
	nacum = RECORD1['nacum']
	#print kx


        rec2 = np.dtype([('dummy1','uint32'), \
                ('jsp','uint32',RECORD1['np']),\
		('yfmap',yfmap,my), \
		('dummy2','uint32')])
	RECORD2=np.fromfile(f,rec2,1)
	nspec = len(RECORD2['jsp'][0])
	#Save y array
	y = RECORD2['yfmap']['y'][0,]
	jsp = RECORD2['jsp'][0,]

	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good for RECORD2...!"
		print RECORD2['dummy1']
	#else:
    		#print "Fiel read correctly :)"
        #Need to create stats dtype

	ntotr = nz*nspec*nvar
	print "ntotr=%s,nx=%s,nz=%s,nspec=%s,nvar=%s" \
		 %(ntotr,nx,nz,nspec,nvar)

	temp1 = np.zeros(ntotr)
	temp2 = np.ndarray(shape=(nz,nspec,nvar),\
		dtype = float, order = 'F')
	SPEC = np.ndarray(shape=(nz,nspec,nvar,nx),\
		dtype = float, order = 'F')
	recspec = np.dtype([('dummy1','uint32'), \
		('data','float32',ntotr), \
		('dummy2','uint32')])
	for i in range(nx):
		readrec=np.fromfile(f,recspec,1)
		if readrec['dummy1'] != readrec['dummy2']:
			print "File read NOT good for plane %s...!" %(i+1) 			
		temp1=readrec['data'][:]
		temp2 = np.reshape(temp1,(nz,nspec,nvar),order='F')
		SPEC[:,:,:,i] = temp2[:,:,:]
	f.close()
	SPEC = SPEC/nacum
	
	#Return all
	return y[jsp-1],kx,kz,SPEC 


def worksta(filename,opt=1):
	""" This function reads the sta file generated by LOMA and so some scaling
	scale = 'dw' or 'dm'
        um,vm,wm,urms,vrms,wrms,uv = readsta(..)
        "..p" variables mean square of the variables (energy)
	(Nvar x array[my]) = readsta (str)
	opt=0 ==> no rhom yet on sta. opt=1 ==> rhom on sta
	"""
	import scipy as sc
	import numpy as np
	import pylab
	f = open(filename,'rb')
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('my','uint32'),
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
		('a0','float32'),('nacum','uint32'), \
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)
	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for RECORD1...!"
	#else:
    	#	print "Fiel read correctly :)"
        #Save "my"
	my=RECORD1['my']

        rec2 = np.dtype([('dummy1','uint32'), \
		('yfmap',yfmap,my), \
		('dummy2','uint32')])
	RECORD2=np.fromfile(f,rec2,1)

	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good for RECORD2...!"
	#else:
    		#print "Fiel read correctly :)"
	#Save y array
	y = RECORD2['yfmap']['y'][0,]
	#print y
        #Need to create stats dtype
	if opt==0: #Last version without rhom
		stats = np.dtype([('um','float64'), \
			('vm' ,'float64'), \
			('wm' ,'float64'), \
	                ('rum','float64'), \
	                ('rvm','float64'), \
	                ('rwm','float64'), \
	                ('Tm' ,'float64'), \
	                ('w1m','float64'), \
	                ('w2m','float64'), \
	                ('w3m','float64'), \
	                ('up' ,'float64'), \
	                ('vp' ,'float64'), \
	                ('wp' ,'float64'), \
	                ('uvr','float64'), \
	                ('uwr','float64'), \
	                ('vwr','float64'), \
	                ('w1p','float64'), \
	                ('w2p','float64'), \
	                ('w3p','float64'), \
	                ('Tp','float64'), \
	                ('ep','float64'), \
	                ('ruu','float64'), \
	                ('ruv','float64'), \
	                ('ruw','float64'), \
	                ('rvv','float64'), \
	                ('rvw','float64'), \
	                ('rww','float64')])

	elif opt==1:
		stats = np.dtype([('um','float64'), \
			('vm' ,'float64'), \
			('wm' ,'float64'), \
	                ('rum','float64'), \
	                ('rvm','float64'), \
	                ('rwm','float64'), \
	                ('Tm' ,'float64'), \
	                ('w1m','float64'), \
	                ('w2m','float64'), \
	                ('w3m','float64'), \
	                ('up' ,'float64'), \
	                ('vp' ,'float64'), \
	                ('wp' ,'float64'), \
	                ('uvr','float64'), \
	                ('uwr','float64'), \
	                ('vwr','float64'), \
	                ('w1p','float64'), \
	                ('w2p','float64'), \
	                ('w3p','float64'), \
	                ('Tp','float64'), \
	                ('ep','float64'), \
	                ('ruu','float64'), \
	                ('ruv','float64'), \
	                ('ruw','float64'), \
	                ('rvv','float64'), \
	                ('rvw','float64'), \
	                ('rww','float64'), \
	                ('rhom','float64')])
	elif opt==2:
		stats = np.dtype([('um','float64'), \
			('vm' ,'float64'), \
			('wm' ,'float64'), \
	                ('rum','float64'), \
	                ('rvm','float64'), \
	                ('rwm','float64'), \
	                ('Tm' ,'float64'), \
	                ('w1m','float64'), \
	                ('w2m','float64'), \
	                ('w3m','float64'), \
	                ('them','float64'), \
	                ('mum','float64'), \
	                ('up' ,'float64'), \
	                ('vp' ,'float64'), \
	                ('wp' ,'float64'), \
	                ('uvr','float64'), \
	                ('uwr','float64'), \
	                ('vwr','float64'), \
	                ('w1p','float64'), \
	                ('w2p','float64'), \
	                ('w3p','float64'), \
	                ('Tp','float64'), \
	                ('ep','float64'), \
	                ('ruu','float64'), \
	                ('ruv','float64'), \
	                ('ruw','float64'), \
	                ('rvv','float64'), \
	                ('rvw','float64'), \
	                ('rww','float64'), \
	                ('rhom','float64'), \
	                ('thep','float64'), \
	                ('theup','float64')])



	nacum = RECORD1['nacum']
	#Create type "recstats"
	recstats = np.dtype([('dummy1','uint32'), \
	 	('data',stats,my), \
		('dummy2','uint32')])

	DATA=np.fromfile(f,recstats,1)
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for STATS...!"
	#else:
    	#	print "Fiel read correctly for STATS :)"
	f.close()

	#SAVE VARIABLES
	if opt > 0:
        	ruu = DATA['data']['ruu'][0,:]/nacum
        	ruv = DATA['data']['ruv'][0,:]/nacum
        	ruw = DATA['data']['ruw'][0,:]/nacum
        	rvv = DATA['data']['rvv'][0,:]/nacum
        	rvw = DATA['data']['rvw'][0,:]/nacum
        	rww = DATA['data']['rww'][0,:]/nacum
		rhom = DATA['data']['rhom'][0,:]/nacum
	else:
		ruu = DATA['data']['up'][0,:]/nacum
		rvv = DATA['data']['vp'][0,:]/nacum
		rww = DATA['data']['wp'][0,:]/nacum
		ruv = DATA['data']['uvr'][0,:]/nacum
		ruw = DATA['data']['uwr'][0,:]/nacum
		rvw = DATA['data']['vwr'][0,:]/nacum
		rhom = DATA['data']['Tm'][0,:]/nacum #only s=1.0 has some old files with no rhom
        um  = DATA['data']['um'][0,:]/nacum 
        vm  = DATA['data']['vm'][0,:]/nacum
        wm  = DATA['data']['wm'][0,:]/nacum
	uu  = DATA['data']['up'][0,:]/nacum
	vv  = DATA['data']['vp'][0,:]/nacum
	ww  = DATA['data']['wp'][0,:]/nacum
        rum = DATA['data']['rum'][0,:]/nacum
        rvm = DATA['data']['rvm'][0,:]/nacum
        rwm = DATA['data']['rwm'][0,:]/nacum
        uv  = DATA['data']['uvr'][0,:]/nacum
        uw  = DATA['data']['uwr'][0,:]/nacum
        vw  = DATA['data']['vwr'][0,:]/nacum
	w3p = DATA['data']['w3p'][0,:]/nacum
	w3m = DATA['data']['w3m'][0,:]/nacum
	w2p = DATA['data']['w2p'][0,:]/nacum
	w2m = DATA['data']['w2m'][0,:]/nacum
	w1p = DATA['data']['w1p'][0,:]/nacum
	w1m = DATA['data']['w1m'][0,:]/nacum
        Tm = DATA['data']['Tm'][0,:]/nacum
	TT = DATA['data']['Tp'][0,:]/nacum
        ep = DATA['data']['ep'][0,:]/nacum - w3m**2 #take the mean flow dissipation
	Re = RECORD1['Re'] 
	if opt==2:
        	them  = DATA['data']['them'] [0,:]/nacum 
        	mum   = (DATA['data']['mum']  [0,:]/nacum)/Re
		thep  = DATA['data']['thep'] [0,:]/nacum 
		theup = DATA['data']['theup'][0,:]/nacum 
	else:
		mum = np.ones([len(y)])/Re #create vector

	#mu = 1.0/Re #If rho0 = 1.0 and mu constant
	my = len(y);
	nu = mum/rhom
	dyeta = 0.0*rhom #initializing
	#CALCULATIONS AND INTERPOLATIONS:
	dm = calcdmcomp(y,um,rum,rhom) #Momentum thickness
	dw = calcdwcomp(y,um,rum,rhom) #Vorticity thickness
        urms = (abs(uu-um**2.0))**0.5
        vrms = (abs(vv-vm**2.0))**0.5
        wrms = (abs(ww-wm**2.0))**0.5
	Trms = (abs(TT-Tm**2.0))**0.5
        w1rms = (abs(w1p-w1m**2.0))**0.5 #vorticity rms
        w2rms = (abs(w2p-w2m**2.0))**0.5
        w3rms = (abs(w3p-w3m**2.0))**0.5
        uv   = uv - um*vm #perturbations
        uw   = uw - um*wm
        vw   = vw - vm*wm
	#note that mean(rhou*u) = mean(rhou)*mean(u)+mean(rhou'u')
	#ruu  = ruu - rum*um #removing mean fluxes
	#rvv  = rvv - rvm*vm
	#rww  = rww - rwm*wm
	#ruv  = ruv - rum*vm
	#find ml effective region
	#pos_ml = where(abs(y)/dm<3.0)#Only valid for points within ML
	#dyeta = np.zeros([len(pos_ml)]) #initializing
	pos_ml = where(abs(y)/dw<1.0)#Only valid for points within ML
	#eta = np.zeros([len(pos_ml)]) #initializing
	#eta = 0.0*rhom[pos_ml] #initializing
	#Definitions from SARKAR et all 1989:the analysis and modelling...
	#epD = dilatational dissipation
	epD = uv*0.0 #before STATS opt 2, no epD was defined
	XI  = uv*0.0
	#epD = 2.0*mu*(der2(y,vrms**2.0))/rhom #no homegeinity condition
	#epS = ep/rhom #solenoidal part !pseudodissipation,
	epS = mum*ep/rhom #solenoidal part !pseudodissipation,
	epR = epS #Solenoidal definition is enough as seen
	if opt==2:
		thep = thep - them**2 
		theup = der1(y,theup-them*vm)
		epD  = mum/rhom*(4.0/3.0*thep-4*theup)
		epR = mum*ep/rhom + epD
		#Compressible fraction of the dissipation rate: epD/ep (without no-homo part)
		XI[pos_ml]  = 4.0/3.0*thep[pos_ml]/(4.0/3.0*thep[pos_ml]+ep[pos_ml])
	#epR = epD+epS #Reynolds definition turbulent dissipation
	EPS = np.trapz(epR,y) #integral of energy dissipation rate
	#kolmogorov scale at mid
	eta = (nu**3/epR)**0.25 #eta=(nu^3/eps)^(1/4) 
	for j in range(0,my-1):
		if eta[j]==0.0:
			dyeta[j]=0.0
		else:
			if j==0:
				dyeta[j] =1.0/(2.0*eta[j]/(y[j+1]-y[j]))
			elif j==len(y)-1:
				dyeta[j] =1.0/(2.0*eta[j]/(y[j]-y[j-1]))
			else:
				dyeta[j] =1.0/(2.0*eta[j]/(y[j+1]-y[j-1]))
	#Mean momentum perturbations
	#ru = rum - rhom*um
	#rv = rvm - rhom*vm
	#rw = rwm - rhom*wm
	#Turbulent kinetic energy according to favre aver.
	#kfav = 0.5/rhom*(ruu+rvv+rww-2*((ru**2.0+rv**2.0+rw**2.0)/rhom))

	#Turbulent reynolds stresses as Favre
	R11 = (ruu-rum*rum/rhom)/rhom
	R22 = (rvv-rvm*rvm/rhom)/rhom
	R33 = (rww-rwm*rwm/rhom)/rhom
	R12 = (ruv-rum*rvm/rhom)/rhom
	rho0 = 0.5*(rhom[0]+rhom[-1])
	dmpoint = -2.0/rho0*np.trapz(rhom*R12*der1(y,rum/rhom),y)
	if np.max(rhom)==np.mean(rhom): #density ratio = 1.0
		R11 = uu-um**2.0
		R22 = vv-vm**2.0
		R33 = ww-wm**2.0
		R12 = uv-um*vm 
		dmpoint = -2.0*np.trapz(R12*der1(y,um),y)
	#Turbulent kinetic energy and Relambda from Pantano&Sarkar.
	k = 0.5*(R11+R22+R33) 
	Relambda = 0.0*k
	llambda = 0.0*k
	llambda[pos_ml] = (2.0/3.0*k[pos_ml]*5.0*nu[pos_ml]/epR[pos_ml])**0.5
	Relambda[pos_ml] = 2*k[pos_ml]*(5.0/(epR[pos_ml]*nu[pos_ml]))**0.5
	#Relambda = 2.0*kfav*(5.0*Re/epR)**0.5
	if len(eta) == 0:
		etamin=0
	else:
		etamin=np.min(eta)

	


	#Return all
	return  {'time':RECORD1['time'][0], 'y':y,'epS':epS,'epD':epD,'epR':epR, 'EPS':EPS,'Relambda':Relambda.max(),  \
	'Re':Re,'alp':RECORD1['alp'],'bet':RECORD1['bet'], 'Relambda_v':Relambda ,'dyeta':dyeta, \
	'R11':R11,'R22':R22,'R33':R33,'R12':R12,'dm':dm,'dw':dw,'dmpoint':dmpoint,'etamin':etamin,'eta':eta, \
	'um':um,'vm':vm,'wm':wm,'urms':urms,'vrms':vrms,'wrms':wrms,'uv':uv,'Trms':Trms,'TT':TT, \
	'w1rms':w1rms,'w2rms':w2rms,'w3rms':w3rms,'w1m':w1m,'w2m':w2m,'w3m':w3m,'mum':mum,'Relambday':Relambda, \
	'Tm':Tm,'rum':rum,'rvm':rvm,'rwm':rwm,'rhom':rhom,'XI':XI,'nu':nu,'nacum':nacum,'k':k,'llambda':llambda}
	


#===============================================================#

def workspe(filename,plane=-1,nvar=12):
	""" This function reads the spe file generated by LOMA
#(yjsp,kx,kz,SPEC) = readspe(filename,nvar)
	(y(jsp),kx,kz,[nz/2,nspec,nvar,nx/2] = readspe (str,int)
	"""
	import scipy as sc
	import numpy as np
	import pylab
	f = open(filename,'rb')
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
                ('nx','uint32'),('ny','uint32'), \
                ('nz','uint32'),('np','uint32'), \
		('nacum','uint32'), \
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)
	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for RECORD1...!"
	#else:
    	#	print "Fiel read correctly :)"
        #Save "my"
	my = RECORD1['ny']
	nx = RECORD1['nx'][0]/2
	nz = (RECORD1['nz'][0]+1)/2
	kx = range(nx)*RECORD1['alp']
	kz = range(nz)*RECORD1['bet']
	nacum = RECORD1['nacum']
	time = RECORD1['time']
	#print kx


        rec2 = np.dtype([('dummy1','uint32'), \
                ('jsp','uint32',RECORD1['np']),\
		('yfmap',yfmap,my), \
		('dummy2','uint32')])
	RECORD2=np.fromfile(f,rec2,1)
	nspec = len(RECORD2['jsp'][0])
	#Save y array
	y = RECORD2['yfmap']['y'][0,]
	jsp = RECORD2['jsp'][0,]

	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good for RECORD2...!"
		print RECORD2['dummy1']
	#else:
    		#print "Fiel read correctly :)"
        #Need to create stats dtype

	ntotr = nz*nspec*nvar
	print "ntotr=%s,nx=%s,nz=%s,nspec=%s,nvar=%s" \
		 %(ntotr,nx,nz,nspec,nvar)

	temp1 = np.zeros(ntotr)
	temp2 = np.ndarray(shape=(nz,nspec,nvar),\
		dtype = float, order = 'F')
	SPEC = np.ndarray(shape=(nz,nspec,nvar,nx),\
		dtype = float, order = 'F')
	recspec = np.dtype([('dummy1','uint32'), \
		('data','float32',ntotr), \
		('dummy2','uint32')])
	for i in range(nx):
		readrec=np.fromfile(f,recspec,1)
		if readrec['dummy1'] != readrec['dummy2']:
			print "File read NOT good for plane %s...!" %(i+1) 			
		temp1=readrec['data'][:]
		temp2 = np.reshape(temp1,(nz,nspec,nvar),order='F')
		SPEC[:,:,:,i] = temp2[:,:,:]
	f.close()
	SPECplane = np.squeeze(SPEC[:,plane,:,:])/nacum
	
	#Return all
	return {'time':time,'yjsp':y[plane],'kx':kx,'kz':kz,'SPEC':SPECplane,'jsp':jsp,'y':y} 



#def hdfconvertxz(filename):
#	""" 
#	This function reads the XZ field generated by TOFIS fortran program and 
#	converts it to same filename + .h5, using gzip compression. Furthermore
#	this new file includes: xvec,zvec,y,jspecy,Re,time and the field itself.
#	{y,xvec,zvec,time} = readfieldxz(filename)
#	"""
#	import scipy as sc
#	import numpy as np
#  	import h5py
#	import pylab
#	import os.path
#	
#	if os.path.isfile(filename+'.h5')==True:
#		print "This file already exists man!Nothing done."
#		return
#	f = open(filename,'rb')
#	#Create dtypes for proper reading from Fortran unformatted
#	# binary file
#	#Declaring types
#	yfmap = np.dtype([('y','float64'),('fmap','float64')])
#	#uw00 = np.dtype([('u00','float32'),('w00','float32')])
#	rec1 = np.dtype([('dummy1','uint32'), \
#                ('time','float32'),('Re','float32'), \
#		('alp','float32'),('bet','float32'), \
#		('mgalx','uint32'),('my','uint32'),('mgalz','uint32'),\
#		('nspec','uint32'),('nacum','uint32'),\
#		('dummy2','uint32')])
#	#Read first record
#	RECORD1=np.fromfile(f,rec1,1)
#
#	#Check if first record is ok...
#	if RECORD1['dummy1'] != RECORD1['dummy2']:
#		print "File read not good...!"
#	else:
#    		print "Fiel read correctly :)"
#
#	mgalx=RECORD1['mgalx'][0]
#	my=RECORD1['my'][0]
#	mgalz=RECORD1['mgalz'][0]
#	nspec=RECORD1['nspec'][0]
#
#
#	print "nspec= %s" % nspec
#
#	rec2 = np.dtype([('dummy1','uint32'), \
# 		('jspecy','uint32',nspec),  \
#		('yfmap',yfmap,my),	\
#		('dummy2','uint32')])
#	#READ record2
#	RECORD2=np.fromfile(f,rec2,1)
#
#	#Check if record2 is ok...
#	if RECORD2['dummy1'] != RECORD2['dummy2']:
#		print "File read not good...!"
#	else:
#    		print "Fiel read correctly :)"
#
#	#Save y vector amd jspecy
#	y = RECORD2['yfmap']['y'][0,]
#        jspecy = RECORD2['jspecy'][0,]
#
#
#
#	#Create type "recplane"
#	recplaney = np.dtype([('dummy1','uint32'), \
#                 ('data','float32',mgalx*mgalz), \
#                 ('dummy2','uint32')])
#
#	#Read all planes Y info
#	FIELD1=np.ndarray(shape=(nspec,mgalx*mgalz),\
#	dtype=float, order='F')
#	for j in range(nspec):
#    		readrec = np.fromfile(f,recplaney,1)
#    		planeydata = readrec['data']
#    		#planeydata.shape=(mgalz,mgalx)
#		FIELD1[j,:] = planeydata[:]
#		#FIELD1[j,:,:] = planeydata[:,:]
#	f.close()
#        #Create vector X and Z
#	Lx = 2*3.1415/RECORD1['alp']
#	Lz = 2*3.1415/RECORD1['bet']
#	x = np.linspace(-Lx/2.,Lx/2.,mgalx)
#	z = np.linspace(-Lz/2.,Lz/2.,mgalz)
#
#
#        FIELD1.shape=(nspec,mgalz,mgalx)
#	FIELD1=FIELD1.transpose(2,0,1)
#	hf5 = h5py.File(filename+'.h5','w')
#	hf5.create_dataset('field',data=FIELD1,compression='gzip')
#	hf5.create_dataset('xvec',data=x)
#	hf5.create_dataset('zvec',data=z)
#	hf5.create_dataset('y',data=y)
#	hf5.create_dataset('time',data=RECORD1['time'])
#	hf5.create_dataset('Re',data=RECORD1['Re'])
#	hf5.create_dataset('jspecy',data=jspecy)
#	hf5.close()
#	
#        del FIELD1
#        print 'Data from time = %s' % RECORD1['time']
#        print 'mgalx = %s, my = %s, mgalz = %s' % (RECORD1['mgalx'], \
#		RECORD1['my'],RECORD1['mgalz'])
#	#Return y, FIELD
#	return {'xvec':x,'zvec':z,'y':y,'time':RECORD1['time']}  
#
##===============================================================#
#

#===============================================================#

def workspeall(filename,nvar=12):
	""" This function reads the spe file generated by LOMA 
            returning all planes
#(yjsp,kx,kz,SPEC) = readspe(filename,nvar)
	(y(jsp),kx,kz,[nz/2,nspec,nvar,nx/2] = readspe (str,int)
	"""
	import scipy as sc
	import numpy as np
	import pylab
	f = open(filename,'rb')
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
                ('nx','uint32'),('ny','uint32'), \
                ('nz','uint32'),('np','uint32'), \
		('nacum','uint32'), \
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)
	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for RECORD1...!"
	#else:
    	#	print "Fiel read correctly :)"
        #Save "my"
	my = RECORD1['ny']
	nx = RECORD1['nx'][0]/2
	nz = (RECORD1['nz'][0]+1)/2
	kx = range(nx)*RECORD1['alp']
	kz = range(nz)*RECORD1['bet']
	nacum = RECORD1['nacum']
	time = RECORD1['time']
	#print kx


        rec2 = np.dtype([('dummy1','uint32'), \
                ('jsp','uint32',RECORD1['np']),\
		('yfmap',yfmap,my), \
		('dummy2','uint32')])
	RECORD2=np.fromfile(f,rec2,1)
	nspec = len(RECORD2['jsp'][0])
	#Save y array
	y = RECORD2['yfmap']['y'][0,]
	jsp = RECORD2['jsp'][0,]

	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good for RECORD2...!"
		print RECORD2['dummy1']
	#else:
    		#print "Fiel read correctly :)"
        #Need to create stats dtype

	ntotr = nz*nspec*nvar
	print "ntotr=%s,nx=%s,nz=%s,nspec=%s,nvar=%s" \
		 %(ntotr,nx,nz,nspec,nvar)

	temp1 = np.zeros(ntotr)
	temp2 = np.ndarray(shape=(nz,nspec,nvar),\
		dtype = float, order = 'F')
	SPEC = np.ndarray(shape=(nz,nspec,nvar,nx),\
		dtype = float, order = 'F')
	recspec = np.dtype([('dummy1','uint32'), \
		('data','float32',ntotr), \
		('dummy2','uint32')])
	for i in range(nx):
		readrec=np.fromfile(f,recspec,1)
		if readrec['dummy1'] != readrec['dummy2']:
			print "File read NOT good for plane %s...!" %(i+1) 			
		temp1=readrec['data'][:]
		temp2 = np.reshape(temp1,(nz,nspec,nvar),order='F')
		SPEC[:,:,:,i] = temp2[:,:,:]
	f.close()
	SPECplane = np.squeeze(SPEC[:,:,:,:])/nacum
	
	#Return all
	return {'time':time,'yjsp':y[jsp],'kx':kx,'kz':kz,'SPEC':SPEC,'jsp':jsp,'y':y} 



#Read sta for lomaHZ

def workstaHZ(filename):
	""" This function reads the sta file generated by LOMA and so some scaling
	scale = 'dw' or 'dm'
        um,vm,wm,urms,vrms,wrms,uv = readsta(..)
        "..p" variables mean square of the variables (energy)
	(Nvar x array[my]) = readsta (str)
	opt=0 ==> no rhom yet on sta. opt=1 ==> rhom on sta
	"""
	import scipy as sc
	import numpy as np
	import pylab
	f = open(filename,'rb')
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('my','uint32'),
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
		('a0','float32'),('nacum','uint32'), \
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)
	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for RECORD1...!"
	#else:
    	#	print "Fiel read correctly :)"
        #Save "my"
	my=RECORD1['my']

        rec2 = np.dtype([('dummy1','uint32'), \
		('yfmap',yfmap,my), \
		('dummy2','uint32')])
	RECORD2=np.fromfile(f,rec2,1)

	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good for RECORD2...!"
	#else:
    		#print "Fiel read correctly :)"
	#Save y array
	y = RECORD2['yfmap']['y'][0,]
	#print y
        #Need to create stats dtype
	stats = np.dtype([('um','float64'), \
		('vm' ,'float64'), \
		('wm' ,'float64'), \
                ('rum','float64'), \
                ('rvm','float64'), \
                ('rwm','float64'), \
                ('w1m','float64'), \
                ('w2m','float64'), \
                ('w3m','float64'), \
                ('them','float64'), \
                ('Tm','float64'), \
                ('Hm','float64'), \
                ('Zm','float64'), \
                ('up' ,'float64'), \
                ('vp' ,'float64'), \
                ('wp' ,'float64'), \
                ('uvr','float64'), \
                ('uwr','float64'), \
                ('vwr','float64'), \
                ('w1p','float64'), \
                ('w2p','float64'), \
                ('w3p','float64'), \
                ('Tp','float64'), \
                ('ep','float64'), \
                ('ruu','float64'), \
                ('ruv','float64'), \
                ('ruw','float64'), \
                ('rvv','float64'), \
                ('rvw','float64'), \
                ('rww','float64'), \
                ('rhom','float64'), \
                ('thep','float64'), \
 		('theup','float64'),\
		('Hp','float64'),\
		('Zp','float64'),\
		('mum','float64')])



	nacum = RECORD1['nacum']
	#Create type "recstats"
	recstats = np.dtype([('dummy1','uint32'), \
	 	('data',stats,my), \
		('dummy2','uint32')])

	DATA=np.fromfile(f,recstats,1)
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for STATS...!"
	#else:
    	#	print "Fiel read correctly for STATS :)"
	f.close()

	#SAVE VARIABLES
        ruu = DATA['data']['ruu'][0,:]/nacum
       	ruv = DATA['data']['ruv'][0,:]/nacum
       	ruw = DATA['data']['ruw'][0,:]/nacum
       	rvv = DATA['data']['rvv'][0,:]/nacum
       	rvw = DATA['data']['rvw'][0,:]/nacum
       	rww = DATA['data']['rww'][0,:]/nacum
	rhom = DATA['data']['rhom'][0,:]/nacum

        um  = DATA['data']['um'][0,:]/nacum 
        vm  = DATA['data']['vm'][0,:]/nacum
        wm  = DATA['data']['wm'][0,:]/nacum
	uu  = DATA['data']['up'][0,:]/nacum
	vv  = DATA['data']['vp'][0,:]/nacum
	ww  = DATA['data']['wp'][0,:]/nacum
        rum = DATA['data']['rum'][0,:]/nacum
        rvm = DATA['data']['rvm'][0,:]/nacum
        rwm = DATA['data']['rwm'][0,:]/nacum
        uv  = DATA['data']['uvr'][0,:]/nacum
        uw  = DATA['data']['uwr'][0,:]/nacum
        vw  = DATA['data']['vwr'][0,:]/nacum
	w3p = DATA['data']['w3p'][0,:]/nacum
	w3m = DATA['data']['w3m'][0,:]/nacum
	w2p = DATA['data']['w2p'][0,:]/nacum
	w2m = DATA['data']['w2m'][0,:]/nacum
	w1p = DATA['data']['w1p'][0,:]/nacum
	w1m = DATA['data']['w1m'][0,:]/nacum
        Tm = DATA['data']['Tm'][0,:]/nacum
        Hm = DATA['data']['Hm'][0,:]/nacum
        Zm = DATA['data']['Zm'][0,:]/nacum
	TT = DATA['data']['Tp'][0,:]/nacum
	Hp = DATA['data']['Hp'][0,:]/nacum
	Zp = DATA['data']['Zp'][0,:]/nacum
        ep = DATA['data']['ep'][0,:]/nacum - w3m**2 #take the mean flow dissipation
	Re = RECORD1['Re'] 
        them  = DATA['data']['them'] [0,:]/nacum 
        mum   = (DATA['data']['mum'] [0,:]/nacum)/Re
	thep  = DATA['data']['thep'] [0,:]/nacum 
	theup = DATA['data']['theup'][0,:]/nacum 

	#mu = 1.0/Re #If rho0 = 1.0 and mu constant
	my = len(y);
	nu = mum/rhom
	dyeta = 0.0*rhom #initializing
	#CALCULATIONS AND INTERPOLATIONS:
	dm = calcdmcomp2(y,um,rum,rhom) #Momentum thickness
	dw = calcdwcomp(y,um,rum,rhom) #Vorticity thickness
        urms = (abs(uu-um**2.0))**0.5
        vrms = (abs(vv-vm**2.0))**0.5
        wrms = (abs(ww-wm**2.0))**0.5
	Trms = (abs(TT-Tm**2.0))**0.5
	Hp = Hp
	Zp = Zp
        w1rms = (abs(w1p-w1m**2.0))**0.5 #vorticity rms
        w2rms = (abs(w2p-w2m**2.0))**0.5
        w3rms = (abs(w3p-w3m**2.0))**0.5
        uv   = uv - um*vm #perturbations
        uw   = uw - um*wm
        vw   = vw - vm*wm
	#note that mean(rhou*u) = mean(rhou)*mean(u)+mean(rhou'u')
	#ruu  = ruu - rum*um #removing mean fluxes
	#rvv  = rvv - rvm*vm
	#rww  = rww - rwm*wm
	#ruv  = ruv - rum*vm
	#find ml effective region
	pos_ml = where(abs(y)/dw<1.0)#Only valid for points within ML
	#Definitions from SARKAR et all 1989:the analysis and modelling...
	#epD = dilatational dissipation
	epD = uv*0.0 #before STATS opt 2, no epD was defined
	XI  = uv*0.0
	#epD = 2.0*mu*(der2(y,vrms**2.0))/rhom #no homegeinity condition
	#epS = ep/rhom #solenoidal part !pseudodissipation,
	epS = mum*ep/rhom #solenoidal part !pseudodissipation,
	epR = epS #Solenoidal definition is enough as seen
	thep = thep - them**2 
	theup = der1(y,theup-them*vm)
	epD  = mum/rhom*(4.0/3.0*thep-4*theup)
	epR = mum*ep/rhom + epD
		#Compressible fraction of the dissipation rate: epD/ep (without no-homo part)
	XI[pos_ml]  = 4.0/3.0*thep[pos_ml]/(4.0/3.0*thep[pos_ml]+ep[pos_ml])
	#epR = epD+epS #Reynolds definition turbulent dissipation
	EPS = np.trapz(epR,y) #integral of energy dissipation rate
	#kolmogorov scale at mid
	eta = (nu**3/epR)**0.25 #eta=(nu^3/eps)^(1/4) 
	for j in range(0,len(y)-1):
		if eta[j]==0.0:
			dyeta[j]=0.0
		else:
			if j==0:
				dyeta[j] =1.0/(2.0*eta[j]/(y[j+1]-y[j]))
			elif j==len(y)-1:
				dyeta[j] =1.0/(2.0*eta[j]/(y[j]-y[j-1]))
			else:
				dyeta[j] =1.0/(2.0*eta[j]/(y[j+1]-y[j-1]))
	#Mean momentum perturbations
	#ru = rum - rhom*um
	#rv = rvm - rhom*vm
	#rw = rwm - rhom*wm
	#Turbulent kinetic energy according to favre aver.
	#kfav = 0.5/rhom*(ruu+rvv+rww-2*((ru**2.0+rv**2.0+rw**2.0)/rhom))

	#Turbulent reynolds stresses as Favre
	R11 = (ruu-rum*rum/rhom)/rhom
	R22 = (rvv-rvm*rvm/rhom)/rhom
	R33 = (rww-rwm*rwm/rhom)/rhom
	R12 = (ruv-rum*rvm/rhom)/rhom
	rho0 = 0.5*(rhom[0]+rhom[-1])
	dmpoint = -2.0/rho0*np.trapz(rhom*R12*der1(y,rum/rhom),y)
	#Turbulent kinetic energy and Relambda from Pantano&Sarkar.
	k = 0.5*(R11+R22+R33) 
	Relambda = 0.0*k
	Relambda[pos_ml] = 2*k[pos_ml]*(5.0/(epR[pos_ml]*nu[pos_ml]))**0.5
	#Relambda = 2.0*kfav*(5.0*Re/epR)**0.5
	


	#Return all
	return  {'time':RECORD1['time'][0], 'y':y,'epS':epS,'epD':epD,'epR':epR, 'EPS':EPS,'Relambda':Relambda.max(),  \
	'Re':Re,'alp':RECORD1['alp'],'bet':RECORD1['bet'], 'Relambda_v':Relambda ,'dyeta':dyeta, \
	'R11':R11,'R22':R22,'R33':R33,'R12':R12,'dm':dm,'dw':dw,'dmpoint':dmpoint,'etamin':eta.min(),'eta':eta, \
	'um':um,'vm':vm,'wm':wm,'urms':urms,'vrms':vrms,'wrms':wrms,'uv':uv,'Trms':Trms,'Hp':Hp,'Zp':Zp, \
	'w1rms':w1rms,'w2rms':w2rms,'w3rms':w3rms,'w1m':w1m,'w2m':w2m,'w3m':w3m,'mum':mum,'Relambday':Relambda, \
	'Tm':Tm,'Hm':Hm,'Zm':Zm,'rum':rum,'rvm':rvm,'rwm':rwm,'rhom':rhom,'XI':XI,'nu':nu,'k':k,'nacum':nacum}
	

#===============================================================#
def readHZfield(filename,my=513):
	""" This function reads the typical unformatted output generated by LOMACTE
	This case works for output files (classic vor/phi)
	(y,vor,phi,psi,scal,u00,v00,w00) = read4field(filename,my)
	(array,array[nxplanes,mz,my,2],array[nxplanes,mz,my,2],array(my,2),array(my,2) = read4field (str,int)
	"""
	import scipy as sc
	import numpy as np
	import pylab
	f = open(filename,'rb')
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	uvw00 = np.dtype([('u00','float32'),('v00','float32'),('w00','float32')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
		('a0','float32'),\
		('mx','uint32'),('my','uint32'),('mz','uint32'),\
		('yfmap',yfmap,my), \
		('uvw00',uvw00,my),
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)

	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good for RECORD1...!"
	else:
    		print "Fiel read correctly :)"
                print "time = %s " % RECORD1['time']
                print "Re = %s " % RECORD1['Re']
                print "alp = %s " % RECORD1['alp']
                print "bet = %s " % RECORD1['bet']
                print "mx = %s " % RECORD1['mx']
                print "mz = %s " % RECORD1['mz']

	mx=RECORD1['mx']
	my=RECORD1['my']
	mz=RECORD1['mz']
	ntotr  = 5*2*my*mz
        nfield = 2*my*mz
        #plane = np.ndarray(shape=(4,2*my,mz),dtype=float,order='F')
	#Create type "recplane"
	recplane = np.dtype([('dummy1','uint32'), \
                 ( 'data','float32',ntotr), \
                 ('dummy2','uint32')])

	#Read all planes info
	nplanes=mx//2

	temp1 = np.zeros(5*2*my*mz)
	temp2 = np.ndarray(shape=(5,2*my,mz),\
		 dtype=float, order='F')
	VOR = np.ndarray(shape=(2*my,mz,nplanes),\
		 dtype=float, order='F')
	PHI = np.ndarray(shape=(2*my,mz,nplanes),\
		 dtype=float, order='F')
	PSI = np.ndarray(shape=(2*my,mz,nplanes),\
		 dtype=float, order='F')
	SCAL = np.ndarray(shape=(2*my,mz,nplanes),\
		 dtype=float, order='F')
	MFZ = np.ndarray(shape=(2*my,mz,nplanes),\
		 dtype=float, order='F')

	for i in range(nplanes):
		readrec = np.fromfile(f,recplane,1)
		if recplane['dummy1'] != recplane['dummy2']:
			print "File read not good for plane %s ...!" % i
		#else:
    			#print "Fiel read correctly for plane %s :)" % i
		temp1 = readrec['data'][:]
		print "reshaping plane %s" %i
		temp2 = np.reshape(temp1,(5,2*my,mz),order='F')
		VOR[:,:,i] = temp2[0,:,:]
		PHI[:,:,i] = temp2[1,:,:]
		PSI[:,:,i] = temp2[2,:,:]
		SCAL[:,:,i] = temp2[3,:,:]
		MFZ[:,:,i] = temp2[4,:,:]

	f.close()
        print 'Data from time = %s' % RECORD1['time']
        print 'mx = %s, my = %s, mz = %s' % (RECORD1['mx'],RECORD1['my'],RECORD1['mz'])
	#Return y, FIELD
	return {'y':RECORD1['yfmap']['y'][0,:],'vor':VOR, 'phi':PHI, 'psi': PSI, \
	'scal':SCAL,'u00':RECORD1['uvw00']['u00'],'v00':RECORD1['uvw00']['v00'],\
	'w00':RECORD1['uvw00']['w00'],'MFZ':MFZ}

#COMPUTE STATS
def compstatsij(filename):
	""" 
	This function reads the field generated by TOFIS fortran program.
	Computing some stats
	() = compstatsij(filename)
	"""
	import scipy as sc
	import numpy as np
	import pylab
	f = open(filename,'rb')
        planeij = filename[-2:]
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	#uw00 = np.dtype([('u00','float32'),('w00','float32')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
		('mgalx','uint32'),('my','uint32'),('mgalz','uint32'),\
		('nspec','uint32'),('plane','uint32'),\
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)

	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good...!"
	else:
    		print "Fiel read correctly :)"

	mgalx=RECORD1['mgalx'][0]
	my=RECORD1['my'][0]
	mgalz=RECORD1['mgalz'][0]
	nspec=RECORD1['nspec'][0]
	time = RECORD1['time'][0]
	print 'time = %s' % time
                
        plane = RECORD1['plane']

	print "nspec= %s" % nspec

	rec2 = np.dtype([('dummy1','uint32'), \
 		('jspecy','uint32',nspec),  \
		('yfmap',yfmap,my),	\
		('dummy2','uint32')])
	#READ record2
	RECORD2=np.fromfile(f,rec2,1)

	#Check if record2 is ok...
	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good...!"
	else:
    		print "Fiel read correctly :)"

	#Save y vector amd jspecy
	y = RECORD2['yfmap']['y'][0,]
        jspecy = RECORD2['jspecy'][0,]


	#READING Differently depnding of plane typ
        if planeij =='yz':
		print 'Reading PLANE YZ...\n'
		rectZ = np.dtype([('dummy1','uint32'), \
                	 ('data','float32',mgalz), \
                 	('dummy2','uint32')])
		#Read all planes Y info
		FIELD1=np.ndarray(shape=(mgalz,my),\
			 dtype=float, order='F')
		for j in range(my):
    			readrec = np.fromfile(f,rectZ,1)
    			rectZdata = readrec['data']
			FIELD1[:,j] = rectZdata[0,:]
	elif planeij =='yx':
		print 'Reading PLANE YX...\n'
		rectx = np.dtype([('dummy1','uint32'), \
                 ('data','float32',mgalx), \
                 ('dummy2','uint32')])
		#Read all planes Y info
		FIELD1=np.ndarray(shape=(mgalx,my),\
			 dtype=float, order='F')
		for j in range(my):
    			readrec = np.fromfile(f,rectx,1)
    			rectxdata = readrec['data']
			FIELD1[:,j] = rectxdata[0,:]
	elif planeij =='xz':
		print 'Reading VOLUME XZ...\n'
	#Create type "recplane"
		recplaney = np.dtype([('dummy1','uint32'), \
                	 ('data','float32',mgalx*mgalz), \
                 	('dummy2','uint32')])

		#Read all planes Y info
		FIELD1=np.ndarray(shape=(mgalx,nspec,mgalz),\
		dtype=float, order='F')
		#FIELD1=np.ndarray(shape=(nspec,mgalx*mgalz),\
		#dtype=float, order='F')
		print "Number of y planes:%s" % nspec
		Npoints=mgalx*mgalz
		print "Npoints=%s" %Npoints
		vxmean=  0.0*np.arange(0,nspec)
		vstddev= 0.0*np.arange(0,nspec)
		vskew=   0.0*np.arange(0,nspec)
		vkurt=   0.0*np.arange(0,nspec)
		for j in range(nspec):
    			readrec = np.fromfile(f,recplaney,1)
    			planeydata = readrec['data']
			planey=planeydata.ravel()
    			#planeydata.shape=(mgalz,mgalx) 
			#NOW WE START COMPUTING STATS
			#1) Compute average
			xmean=np.mean(planey)
			vxmean[j]=xmean
			order1=0.0
			order2=0.0
			order3=0.0
			order4=0.0
			print "plane at y = %s" % y[jspecy[j]-1]
			#Accumulate stats
			for ipoint in range(0,Npoints):
				order1=order1+planey[ipoint]-xmean	
				order2=order2+(planey[ipoint]-xmean)**2	
				order3=order3+(planey[ipoint]-xmean)**3	
				order4=order4+(planey[ipoint]-xmean)**4	
			stddev=np.sqrt(order2/Npoints)
			vstddev[j]=stddev
			if stddev == 0.0:
				vskew[j]=0.0
				vkurt[j]=0.0
			else:
				vskew[j] = (order3/stddev**3)/Npoints
				vkurt[j] = (order4/stddev**4)/Npoints
		y =y[jspecy-1]
		statsname = filename + '.stats.csv'
		stddevMAX = np.max(vstddev)
		posMAX = np.argmax(vstddev)
		print 'max rms is %s, located  at y=%s' %(stddevMAX,posMAX)
		np.savetxt(statsname,(y,vxmean,vstddev,vskew,vkurt),delimiter=",",header="y,mean,stddev,skew,kurt",fmt='%1.6e')
	f.close()
        #Create vector X and Z
	Lx = 2*3.1415/RECORD1['alp']
	Lz = 2*3.1415/RECORD1['bet']
	x = np.linspace(-Lx/2.,Lx/2.,mgalx)
	z = np.linspace(-Lz/2.,Lz/2.,mgalz)

        print 'Data from time = %s' % RECORD1['time']
        print 'mgalx = %s, my = %s, mgalz = %s' % (RECORD1['mgalx'], \
		RECORD1['my'],RECORD1['mgalz'])
#	#Return y, FIELD
#	if planeij =='yz':
#		return FIELD1,y,x[plane-1],z
#	elif planeij =='yx':
#		return FIELD1,y,x,z[plane-1]
#	elif planeij =='xz':
#		print np.shape(FIELD1)
#                #print len(y[jspecy-1])
#		if nspec==1:
#			y = [y[jspecy-1]]
#			print 'Plane at: %s' % y
#			return FIELD1,y,x,z
#		else:
#			return FIELD1,y[jspecy-1],x,z
#      #Return all
#        return  {'time':RECORD1['time'][0], 'y':y,'epS':epS,'epD':epD,'epR':epR, 'EPS':EPS,'Relambda':Relambda.max(),  \
#        'Re':Re,'alp':RECORD1['alp'],'bet':RECORD1['bet'], 'Relambda_v':Relambda ,'dyeta':dyeta, \
#        'R11':R11,'R22':R22,'R33':R33,'R12':R12,'dm':dm,'dw':dw,'dmpoint':dmpoint,'etamin':etamin,'eta':eta, \
#        'um':um,'vm':vm,'wm':wm,'urms':urms,'vrms':vrms,'wrms':wrms,'uv':uv,'Trms':Trms,'TT':TT, \
#        'w1rms':w1rms,'w2rms':w2rms,'w3rms':w3rms,'w1m':w1m,'w2m':w2m,'w3m':w3m,'mum':mum,'Relambday':Relambda, \
#        'Tm':Tm,'rum':rum,'rvm':rvm,'rwm':rwm,'rhom':rhom,'XI':XI,'nu':nu,'nacum':nacum,'k':k,'llambda':llambda}
#


def statsfromXZ(filename):
	""" 
	This function reads the field generated by TOFIS fortran program.
	Computing some stats
	() = statsfromXZ(filename)
	"""
	import scipy as sc
	import numpy as np
	import pylab
	import fortran2vtk
	f = open(filename,'rb')
        planeij = filename[-2:]
	#Create dtypes for proper reading from Fortran unformatted
	# binary file
	#Declaring types
	yfmap = np.dtype([('y','float64'),('fmap','float64')])
	#uw00 = np.dtype([('u00','float32'),('w00','float32')])
	rec1 = np.dtype([('dummy1','uint32'), \
                ('time','float32'),('Re','float32'), \
		('alp','float32'),('bet','float32'), \
		('mgalx','uint32'),('my','uint32'),('mgalz','uint32'),\
		('nspec','uint32'),('plane','uint32'),\
		('dummy2','uint32')])
	#Read first record
	RECORD1=np.fromfile(f,rec1,1)

	#Check if first record is ok...
	if RECORD1['dummy1'] != RECORD1['dummy2']:
		print "File read not good...!"
	else:
    		print "Fiel read correctly :)"

	mgalx=RECORD1['mgalx'][0]
	my=RECORD1['my'][0]
	mgalz=RECORD1['mgalz'][0]
	nspec=RECORD1['nspec'][0]
                
        plane = RECORD1['plane']

	print "nspec= %s" % nspec

	rec2 = np.dtype([('dummy1','uint32'), \
 		('jspecy','uint32',nspec),  \
		('yfmap',yfmap,my),	\
		('dummy2','uint32')])
	#READ record2
	RECORD2=np.fromfile(f,rec2,1)

	#Check if record2 is ok...
	if RECORD2['dummy1'] != RECORD2['dummy2']:
		print "File read not good...!"
	else:
    		print "Fiel read correctly :)"

	#Save y vector amd jspecy
	y = RECORD2['yfmap']['y'][0,]
        jspecy = RECORD2['jspecy'][0,]

	Lx = 2*3.1415/RECORD1['alp']
	Lz = 2*3.1415/RECORD1['bet']
	x = np.linspace(-Lx/2.,Lx/2.,mgalx)
	z = np.linspace(-Lz/2.,Lz/2.,mgalz)
	#These depends on TOFIS
	#READING Differently depnding of plane typ
	if ((planeij =='yz') or (planeij=='yx')):
		print 'STATS are taken from XZ fields \n'
	else:
		print 'Reading VOLUME XZ...\n'
		#Create type "recplane"
		recplaney = np.dtype([('dummy1','uint32'), \
                	 ('data','float32',mgalx*mgalz), \
                 	('dummy2','uint32')])

		#Read all planes Y info
		print "Number of y planes:%s" % nspec
		Npoints=mgalx*mgalz
		print "Npoints=%s" %Npoints
		vxmean=  0.0*np.arange(0,nspec)
		vstddev= 0.0*np.arange(0,nspec)
		vskew=   0.0*np.arange(0,nspec)
		vkurt=   0.0*np.arange(0,nspec)
		for j in range(nspec):
    			readrec = np.fromfile(f,recplaney,1)
    			planeydata = readrec['data']
			planey = planeydata.ravel()
    			#planeydata.shape=(mgalz,mgalx) 
			#NOW WE START COMPUTING STATS
			#1) Compute average
			xmean=np.mean(planey)
			fmin = np.min(planey)
			fmax = np.max(planey)
			vxmean[j]=xmean
			order1=0.0
			order2=0.0
			order3=0.0
			order4=0.0
			#print "plane at y = %s" % y[jspecy[j]-1]
			#Accumulate stats
			for ipoint in range(0,Npoints):
				order1=order1+planey[ipoint]-xmean	
				order2=order2+(planey[ipoint]-xmean)**2	
				order3=order3+(planey[ipoint]-xmean)**3	
				order4=order4+(planey[ipoint]-xmean)**4	
			stddev=np.sqrt(order2/Npoints)
			vstddev[j]=stddev
			if stddev == 0.0:
				vskew[j]=0.0
				vkurt[j]=0.0
			else:
				vskew[j] = (order3/stddev**3)/Npoints
				vkurt[j] = (order4/stddev**4)/Npoints
		y =y[jspecy-1]
		statsname = filename + '.stats.csv'
		np.savetxt(statsname,(y,vxmean,vstddev,vskew,vkurt),delimiter=",",header="y,mean,stddev,skew,kurt",fmt='%1.6e')
	f.close()
	print "file %s created!" % statsname
        print 'Data from time = %s' % RECORD1['time']
        print 'mgalx = %s, my = %s, mgalz = %s' % (RECORD1['mgalx'], \
		RECORD1['my'],RECORD1['mgalz'])

