from scipy import *
from scipy.integrate import cumtrapz
from scipy.integrate import simps
from scipy.integrate import romb
from math  import *
from pylab import *
from numpy import *
from numba import autojit

def isNaN(num):
    return num != num

def myinteg(y,u):
	result = np.zeros(len(y))
	for j in range(0,len(y)-2):
		result[j+1]=result[j]+0.5*(y[j+1]-y[j])*(u[j+1]+u[j])
	#result = result - 0.5* result[-1]
	return result


def der1(y,u):
	"""computes first derivative using Central Finite Differences
        """
	result = zeros(len(y))
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
	result = zeros(len(y))
	for j in range(len(y)):
		if j==0:
			result[j] = (u[j+2]-2*u[j+1]+u[j])/((y[j+2]-y[j])/2.)**2
		elif j==(len(y)-1):
			result[j] = (u[j-2]-2*u[j-1]+u[j])/((y[j-2]-y[j])/2.)**2
		else:
			result[j]=(u[j+1]-2*u[j]+u[j-1])/((y[j+1]-y[j-1])/2.)**2
        return result

def integ(y,u):
        result = cumtrapz(u,y,initial=0.0)
	#return append(result,result[-1])
	return result

def geninit(my,L,deltaT,option=None):
	""" Generate from my and DTemperature the initial arrays
		(y,ru00,rv00,t00) = geninit(my,L,deltaT,option=None)
	"""
	y = linspace(-L,L,my)
	u00 = zeros(my)
	T00 = zeros(my) 
	rhou00 = zeros(my)
	rhov00 = zeros(my)
	for j in range(my):
		u00[j] = erf(pi**0.5*y[j])
		T00[j] = 1.+deltaT/2.*erf(pi**0.5*y[j])
		rhou00[j] = 1./T00[j]*u00[j]
	return y,rhou00,rhov00,T00

def geninit2(my,L,Tmean,deltaT,option=None):
	""" Generate from my and DTemperature the initial arrays
		(y,ru00,rv00,t00) = geninit(my,L,Tmean,deltaT,option=None)
	"""
	y = linspace(-L,L,my)
	u00 = zeros(my)
	T00 = zeros(my) 
	rhou00 = zeros(my)
	rhov00 = zeros(my)
	u0 = 0.5
	for j in range(my):
		#u00[j] = erf(pi**0.5*y[j])
		u00[j] = u0*tanh(-y[j]/2.)
		T00[j] = Tmean+deltaT/2.*tanh(-0.5*y[j])
		rhou00[j] = 1./T00[j]*u00[j]
	return y,rhou00,rhov00,T00



def dtvisc(y,re):
	""" Estiamtes Dtvisc and Dtconv for umax=1
	"""
	#Estimate Dt
	Dtvisc = re/2.*(diff(y)[0])**2
	Dtconv = diff(y)[0]/1.0
	return Dtvisc, Dtconv

def calcdm(y,u00):
	integ1 = trapz(u00**2.0,y)/(u00[-1]-u00[0])**2.0
	return 0.25*(y[-1]-y[0])-integ1	


#def RK3evolve(my,CFL,nstep,Tmean,deltaT,RKopt, param=0.0,time0=None, y=[], ru00=[], rv00=[], t00=[]):
def RK3evolve(my,CFL,nstep,Tmean,deltaT,RKopt, param=0.0,param0=0.0,verb = False,L=86.0,time0=None, y=[], ru00=[], rv00=[], t00=[]):
	""" This functions evolve rhou00,T00 and updates rhov00 using drho,
	    reproduces 00 mode evolution of loma
		(time,y,rhou00,rhov00,T00,timev,dmv) = evolve(my,CFL,nstep,deltaT,RKopt,time0,y,ru00,rv00,t00)
	"""
	if verb == True:
		print "Running Runge-Kutta evolution..."
	#RKNEW
	if RKopt == 1: #RK 3 
		xi    = [0   ,-17.0/60.0,  -5.0/12.0]
		gamma = [8./15.,    5./12.,  3./4.    ]
		beta  = [gamma[0]-param0,    param,  (param*(gamma[0]*gamma[2]+gamma[1]*gamma[2]+xi[1]*gamma[2]+gamma[0]*xi[2])-gamma[0]*gamma[1]*gamma[2]+gamma[0]*gamma[1]*xi[2])/(param-gamma[0]*gamma[1])]
		alpha = [param0,   gamma[1]+xi[1]-param, gamma[2]+xi[2]-beta[2]   ]
		#beta  = [8./15.,    17/15.,  1/6.    ]
		if verb == True:
			print "using RK3 with Low Storage"

	elif RKopt == 2:
		alpha = [0,    -1.,  1/6.    ]
		gamma = [8./15.,    5./12.,  3./4.    ]
		beta  = [8./15.,    param,  (1-gamma[0]+param*gamma[0])/(param*gamma[0])]
		#beta  = [8./15.,    17/15.,  1/6.    ]
		xi    = [0.0   ,-17.0/60.0,  -5.0/12.0]
		#ibeta = [15./4.,       15.,  6.0      ]
		if verb == True:
			print "Using RK3 new option (not Spalart)"
		
	else:
		#NOTE THAT RKopt==3 gives EXPLICIT RK3
		#RKspalart
		gamma = [8./15,      5./12.,  3./4.    ]
		xi    = [0.0  ,   -17./60.0,  -5.0/12.0]
		alpha = [29./96.,   -3./40., 1./6.   ]
		beta  = [37./160.,   5./24.,  1./6.    ]
		ibeta = [160./37.,  24./5.0,  6.0      ]
		if verb == True:
			print "Using RK3 Spalart, EXP-IMP"
	
	#------------Physics -------------------------#
	ire = 1./160. #inverse of Reynolds
	ipe = 1./160. #inverse of Peclet
	#---------------------------------------------#
	#Initialization of variables (arrays)

	if RKopt == 3 and verb==True:
		print "using RK3 Spalart..."
	time=0.0
	dm=zeros(nstep)
	timev = zeros(nstep)
	T00last = zeros(my) 
	T00last1 = zeros(my) 
	drho = zeros(my) 
	drhodt = zeros(my) 
	drholast = zeros(my) 
	intdrho = zeros(my) 
	DT00 = zeros(my) 
	if time0 is None:
		(y,rhou00,rhov00,T00) = geninit2(my,L,Tmean,deltaT)
	else: #restarting
		print "Restarting from previous run..."
		rhou00 = ru00
		rhov00 = rv00
		T00    = t00
		time   = time0
	#density ratio
	sratio = T00[0]/T00[-1]
	if verb == True:
		print "density ratio = %s" % sratio
	#Start RHS
	rhs1=zeros([my,3]) #RHS rhou00
	rhs4=zeros([my,3]) #RHS T00

	Dtvisc = min(diff(y))**2./(2.*ire)
	#Start istep 
	if RKopt==2: #Update using N(rho) directly (NOT low-storage)
		for istep in range(nstep):
		#Obtain timestep
			Dt = max(1./Dtvisc,max(rhov00*T00)/min(diff(y)))
			Dt = CFL/Dt
			if istep==0 and verb == True:
		#		Dt = 0.01*Dtvisc #first step is important
				print "Initial Dt = %s" % Dt
			for rkstep in range(3):
			#For each rkstep
				#Define the RHS's
				rhs1[:,rkstep] = -der1(y,rhou00*rhov00*T00) +ire * der2(y,rhou00*T00)
				rhs4[:,rkstep] = -rhov00*T00*der1(y,T00) + ipe*T00 * der2(y,T00)   
				#Define Dt*RKparameters
				dtxi = Dt* xi[rkstep]
				dtgamma = Dt* gamma[rkstep]
				dtbeta = Dt* beta[rkstep]
				T00last = T00 #Save T00 before evolution
			#Explicit Evolution:
				if rkstep == 0: #First rkstep
					rhou00  = rhou00 + dtgamma*rhs1[:,rkstep] #Evolve rhou00 
					DT00 = dtgamma*rhs4[:,rkstep] #Save the increment of T00
					T00  = T00 + dtgamma*rhs4[:,rkstep] 
	                        	drhodt = -rhs4[:,rkstep]/T00last**2.0  #First substep
				else:
					rhou00  = rhou00[:] + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
					DT00   = dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
					T00    = T00 + dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
	                        	drhodt=(1.0-xi[rkstep]/beta[rkstep])*(-rhs4[:,rkstep]/T00**2.)+xi[rkstep]/beta[rkstep]*(-rhs4[:,rkstep-1]/T00last**2.) 
				#Option 1, not sure about value of v00 at boundaries:
				rhov00 = trapz(drhodt,y)/(sratio+1.)-cumtrapz(drhodt,y,initial=0.0)
				#Option 2, we know this should be zero
				#rhov00 = -cumtrapz(drhodt,y,initial=0.0) #v00(-Ly) = 0.0
			time += Dt
			dm[istep]=calcdm(y,rhou00*T00)
			timev[istep]=time
			if isNaN(dm[istep]):
				print 'NaN found in this simulation!'
				return time, y, rhou00,rhov00, T00,timev,dm
			#print "istep = %s, rkstep = %s" % (istep, rkstep)
	elif RKopt==1: #Low Storage RK3
		for istep in range(nstep):
		#Obtain timestep
			Dt = max(1./Dtvisc,max(rhov00*T00)/min(diff(y)))
			Dt = CFL/Dt
			if istep==0 and verb == True:
		#		Dt = 0.01*Dtvisc #first step is important
				print "Initial Dt = %s" % Dt
			for rkstep in range(3):
			#For each rkstep
				#Define the RHS's
				rhs1[:,rkstep] = -der1(y,rhou00*rhov00*T00) +ire * der2(y,rhou00*T00)
				rhs4[:,rkstep] = -rhov00*T00*der1(y,T00) + ipe*T00 * der2(y,T00)   
				#Define Dt*RKparameters
				dtxi = Dt* xi[rkstep]
				dtgamma = Dt* gamma[rkstep]
				dtbeta = Dt* beta[rkstep]
				T00last = T00
				if rkstep == 0: #First rkstep
					rhou00  = rhou00 + dtgamma*rhs1[:,rkstep] #Evolve rhou00 
					DT00 = dtgamma*rhs4[:,rkstep] #Save the increment of T00
					T00 = T00 + DT00	
					drhodt = -DT00/(beta[rkstep]*Dt*T00last*T00)-alpha[rkstep]/beta[rkstep]*drhodt
					#drhodt = -rhs4[:,rkstep]/T00last**2.
				else:
					rhou00  = rhou00[:] + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
					DT00   = dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
					T00 = T00 + DT00	
					drhodt = -DT00/(beta[rkstep]*Dt*T00last*T00)-alpha[rkstep]/beta[rkstep]*drhodt
				rhov00 = trapz(drhodt,y)/(sratio+1.)-cumtrapz(drhodt,y,initial=0.0)
			#istep loop
			time += Dt
			dm[istep]=calcdm(y,rhou00*T00)
			timev[istep]=time
			if isNaN(dm[istep]):
				print 'NaN found in this simulation!'
				return time, y, rhou00,rhov00, T00,timev,dm
	elif RKopt==3: #Update using N(rho) directly (NOT low-storage)
		#FIRST STEP is different
		timev[0] = 0.0
		dm[0] = calcdm(y,rhou00*T00)
		Dt = max(1./Dtvisc,max(rhov00*T00)/min(diff(y)))
		Dt = CFL/Dt
		print "Initial Dt = %s" % Dt
		T00last = T00 #Save T00 at step n
		Dtkp1 = 0.0
		for rkstep in range(3):
			Dtkp1 += Dt*(gamma[rkstep]+xi[rkstep])
		#For each rkstep
			#Define the RHS's
			rhs1[:,rkstep] = -der1(y,rhou00*rhov00*T00) +ire * der2(y,rhou00*T00)
			rhs4[:,rkstep] = -rhov00*T00*der1(y,T00) + ipe*T00 * der2(y,T00)   
			#Define Dt*RKparameters
			dtxi = Dt* xi[rkstep]
			dtgamma = Dt* gamma[rkstep]
			dtbeta = Dt* beta[rkstep]
			#Explicit Evolution:
			if rkstep == 0: #First rkstep
				rhou00  = rhou00 + dtgamma*rhs1[:,rkstep] #Evolve rhou00 
				DT00 = dtgamma*rhs4[:,rkstep] #Save the increment of T00
				T00  = T00 + dtgamma*rhs4[:,rkstep] 
			else:
				rhou00  = rhou00[:] + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
				DT00   = dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
				T00    = T00 + dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
			drhodt = ((1.0/T00)-(1.00/T00last))/Dtkp1
			#Option 1, not sure about value of v00 at boundaries:
			rhov00 = trapz(drhodt,y)/(sratio+1.)-cumtrapz(drhodt,y,initial=0.0)
				#Option 2, we know this should be zero
				#rhov00 = -cumtrapz(drhodt,y,initial=0.0) #v00(-Ly) = 0.0
		dm[1]=calcdm(y,rhou00*T00)
		time += Dt
		timev[1]=time
		#End OF FIRST STEP
		#--------------------------------------------------------------#
		for istep in range(1,nstep-1):
			T00last1 = T00last
			T00last = T00 #Save T00 before evolution
			Dtkp1 = 0.0
			Dt1 = Dt #PREVIOUS DT
		#Obtain timestep
			Dt = max(1./Dtvisc,max(rhov00*T00)/min(diff(y)))
			Dt = CFL/Dt
			for rkstep in range(3):
				Dtkp1 += Dt*(gamma[rkstep]+xi[rkstep])
				ALP =  Dt1/Dtkp1
			#For each rkstep
				#Define the RHS's
				rhs1[:,rkstep] = -der1(y,rhou00*rhov00*T00) +ire * der2(y,rhou00*T00)
				rhs4[:,rkstep] = -rhov00*T00*der1(y,T00) + ipe*T00 * der2(y,T00)   
				#Define Dt*RKparameters
				dtxi = Dt* xi[rkstep]
				dtgamma = Dt* gamma[rkstep]
				dtbeta = Dt* beta[rkstep]
			#Explicit Evolution:
				if rkstep == 0: #First rkstep
					rhou00  = rhou00 + dtgamma*rhs1[:,rkstep] #Evolve rhou00 
					T00  = T00 + dtgamma*rhs4[:,rkstep] 
				else:
					rhou00  = rhou00[:] + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
					T00    = T00 + dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
				drhodt = (((1.+ALP)**2.0-1.0)*(1.0/T00)-(1.0+ALP)**2.0*(1.0/T00last)+1.0/T00last1)/((1.+ALP)*Dt1)
				#Option 1, not sure about value of v00 at boundaries:
				rhov00 = trapz(drhodt,y)/(sratio+1.)-cumtrapz(drhodt,y,initial=0.0)
				#Option 2, we know this should be zero
				#rhov00 = -cumtrapz(drhodt,y,initial=0.0) #v00(-Ly) = 0.0
			dm[istep+1]=calcdm(y,rhou00*T00)
			time += Dt
			timev[istep+1]=time
			if isNaN(dm[istep]):
				print 'NaN found in this simulation!'
				return time, y, rhou00,rhov00, T00,timev,dm
			#print "istep = %s, rkstep = %s" % (istep, rkstep)

	if verb == True:
		print "last Dt = %s" % Dt
		print "final time = %s, dm = %s " % (timev[-1],dm[-1])
	return time, y, rhou00,rhov00, T00,timev,dm



		#=======================END OF evolve function ==================================#

def RKCNevolve(my,CFL,nstep,Tmean,deltaT,RKopt, time0=None, y=[], ru00=[], rv00=[], t00=[]):
	""" This functions evolve rhou00,T00 and updates rhov00 using drho,
	    reproduces 00 mode evolution of loma
		(time,y,rhou00,rhov00,T00) = evolve(my,CFL,nstep,deltaT,RKopt,time0,y,ru00,rv00,t00)
	"""
	print "Running Runge-Kutta with Crank-Nicholson evolution..."
	aa = 0.5
	bb = 0.5
	#RKNEW
	if RKopt == 1: #RK order 2
                alpha = [1./5. ,   0.4,  1./6.  ]
                beta  = [-0.1  , 1./2.,  1./6.  ]
                gamma = [1./10.,    5.,  3./4.  ]
                xi    = [0.0   ,  -4.1, -0.4167 ]
                ibeta = [-10.  ,    2.,  6.0    ]
		print "Using RK order 2"

	elif RKopt == 2:
		alpha = [4./15.,    1./15.,  1./6.    ]
		beta  = [4./15.,    1./15.,  1./6.    ]
		gamma = [8./15.,    5./12.,  3./4.    ]
		xi    = [0.0   ,-17.0/60.0,  -5.0/12.0]
		ibeta = [15./4.,       15.,  6.0      ]
		print "Using new RK3 (not Spalart)"
		
	else:
	#RKspalart
		gamma = [8./15,      5./12.,  3./4.    ]
		xi    = [0.0  ,   -17./60.0,  -5.0/12.0]
		alpha = [29./96.,   -3./40., 1./6.   ]
		beta  = [37./160.,   5./24.,  1./6.    ]
		ibeta = [160./37.,  24./5.0,  6.0      ]
		print "Using Spalart RK3 EXP-IMP"
	
	#------------Physics -------------------------#
	ire = 1./160. #inverse of Reynolds
	ipe = 1./160. #inverse of Peclet
	L = 10.
	#---------------------------------------------#
	#Initialization of variables (arrays)
	time=0.0
	dm=zeros(nstep)
	timev = zeros(nstep)
	T00last = zeros(my) 
	drho = zeros(my) 
	intdrho = zeros(my) 
	DT00 = zeros(my) 
	if time0 is None:
		(y,rhou00,rhov00,T00) = geninit2(my,L,Tmean,deltaT)
	else: #restarting
		print "Restarting from previous run..."
		rhou00 = ru00
		rhov00 = rv00
		T00    = t00
		time   = time0
	#Start RHS
	rhs1=zeros([my,3]) #RHS rhou00
	rhs4=zeros([my,3]) #RHS T00
	figure()

	Dtvisc = min(diff(y))**2/(2*ire)
	#Start istep 
	for istep in range(nstep):
		#Obtain timestep
		Dt = max(1./Dtvisc,max(rhov00*T00)/min(diff(y)))
		Dt = CFL/Dt
		DT00 = zeros(my)
		for rkstep in range(3):
		#For each rkstep
			#RHS
			rhs1[:,rkstep] = -der1(y,rhou00*rhov00*T00) +ire * der2(y,rhou00*T00)
			rhs4[:,rkstep] = -rhov00*T00*der1(y,T00) + ipe * der2(y,T00)   
			#Dt * RKparameters
	                dtxi = Dt* xi[rkstep]
			dtgamma = Dt* gamma[rkstep]
			dtbeta = Dt* beta[rkstep]
			#Explicit Evolution:
			if rkstep == 0: #First rkstep
				T00last = T00
				rhou00  = rhou00 + dtgamma*rhs1[:,rkstep] 
				T00     = T00 + dtgamma*rhs4[:,rkstep] 
				DT00    = DT00 + dtgamma*rhs4[:,rkstep]
			else:
				rhou00  = rhou00 + dtgamma*rhs1[:,rkstep] + dtxi*rhs1[:,rkstep-1]
				T00     = T00 + dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]
				DT00    = DT00 + dtgamma*rhs4[:,rkstep] + dtxi*rhs4[:,rkstep-1]		
				#We need the total increment of T00 por each step (full 3 substeps)
		#end of rkstep loop
		#Out of RKstep loop
		drho = -DT00/(T00last*T00)
	        intdrho = integ(y,2*drho/Dt) -0.5*trapz(2*drho/Dt,y) #integrates drho
		# use Crank-Nicholson
		rhov00 = -rhov00-intdrho
		#istep loop
		time += Dt

		plot(y,rhov00)
		suptitle('rhov00 evolution', fontsize=14)
		xlabel('y (vertical axis)')
		ylabel('rhov00')

		dm[istep]=calcdm(y,rhou00*T00)
		timev[istep]=time
	print "final time = %s, dm = %s " % (timev[-1],dm[-1])
	return time, y, rhou00,rhov00, T00,timev,dm
	

def ABCNevolve(my,CFL,nstep,Tmean,deltaT,RKopt=None, time0=None, y=[], ru00=[], rv00=[], t00=[]):
	""" This functions evolve rhou00,T00 and updates rhov00 using drho,
	    reproduces 00 mode evolution of loma
		(time,y,rhou00,rhov00,T00,timev,dmv) = evolve(my,CFL,nstep,deltaT,RKopt,time0,y,ru00,rv00,t00)
	"""
	print "Running Adam-Basfort with Crank-Nicholson evolution..."
	aa = 0.5
	bb = 0.5
	#------------Physics -------------------------#
	ire = 1./160. #inverse of Reynolds
	ipe = 1./160. #inverse of Peclet
	L = 84.
	#---------------------------------------------#
	#Initialization of variables (arrays)
	time=0.0
	dm=zeros(nstep)
	timev = zeros(nstep)
	T00last = zeros(my) 
	drho = zeros(my) 
	DT00 = zeros(my) 
	intdrho = zeros(my) 
	figure()
	if time0 is None:
		(y,rhou00,rhov00,T00) = geninit2(my,L,Tmean,deltaT)
	else: #restarting
		print "Restarting from previous run..."
		rhou00 = ru00
		rhov00 = rv00
		T00    = t00
		time   = time0
	#Density ratio
	sratio = T00[0]/T00[-1]
	#Start RHS
	rhs1=zeros(my) #RHS rhou00
	rhs1last=zeros(my) #RHS rhou00
	rhs4=zeros(my) #RHS T00
	rhs4last=zeros(my) #RHS T00
	#Dt visc definition
	Dtvisc = min(diff(y))**2/(2*ire)
	#Start istep 
	for istep in range(nstep):
		#Obtain timestep
		Dt = max(1./Dtvisc,max(rhov00*T00)/min(diff(y)))
		Dt = CFL/Dt
		T00last = T00
		rhs1last = rhs1
		rhs4last = rhs4   
		rhs1 = -der1(y,rhou00*rhov00*T00) +ire * der2(y,rhou00*T00)
		rhs4 = -rhov00*T00*der1(y,T00) + ipe *T00* der2(y,T00)   
		if istep == 0: #first step solved using Euler
			rhou00 = rhou00 + rhs1*Dt	
			T00    = T00 + rhs4*Dt
			drho    = -rhs4/T00**2.
		else:
			#A-B Explicit Evolution:
			rhou00  = rhou00 + 3./2.*Dt*rhs1 - 1./2.*Dt*rhs1last
			T00     =    T00 + 3./2.*Dt*rhs4 - 1./2.*Dt*rhs4last 
			drho    =  -2.0*rhs4/T00**2.+rhs4last/T00last**2.0
		rhov00 = trapz(drho,y)/(sratio+1.)-cumtrapz(drho,y,initial=0.0)
		time += Dt
		dm[istep]=calcdm(y,rhou00*T00)
		timev[istep]=time
	print "final time = %s, dm = %s " % (timev[-1],dm[-1])
	return time, y, rhou00,rhov00, T00,timev,dm

#Manual cumtrapz
def mycumtrapz(y,f):
	"""make cumtrapz same way used on Loma Code
	  returns cumintegral
	"""
	cumintegral = zeros(len(f))
	cumintegral[0] = 0.0
	for j in range(len(y)-1):
    		cumintegral[j+1]=cumintegral[j]+0.5*(y[j+1]-y[j])*(f[j+1]+f[j])
	return cumintegral
	

#TEST FUNCTIONS
def test():
	#TEST derivatives
	my = 1025
	y = linspace(-2*pi,2*pi,my)
	f = sin(y)
	df_anal = cos(y)
	d2f_anal = -sin(y)
	df = der1(y,f)
	d2f = der2(y,f)
	print "Testing derivatives with f = sin(y)"
	figure()
	plot(y,df_anal,'r',y,df,'b--')
	suptitle('First derivative comparition')
	xlabel('y')
	figure()
	plot(y,d2f_anal,'r',y,d2f,'b--')
	suptitle('Second derivative comparition')
	xlabel('y')
	e1 = sum(abs(df-df_anal))/my
	e2 = sum(abs(d2f-d2f_anal))/my
	assert  e1 < 1e-3, \
		"First derivative not good enough. Error = %20.15e" % e1
	assert  e2 < 1e-3, \
		"Second derivative not good enough.Error = %20.15e" % e2

	#TEST integ1
	intf_anal = -cos(y)
	intf = integ(y,f) - 1.0 #-1 is the primitive at j=0	
	intfcumtrapz = cumtrapz(f,y,initial = 0)-1
	e3 = sum(abs(intf-intf_anal))/my
	figure()
	plot(y,intf_anal,'r',y,intf,'b--',y,intfcumtrapz,'r--')
	suptitle('Integration')
	xlabel('y')
	assert e3 < 1e-3, \
		"Itegration looks like SHIT. Error = %20.15e" %e3
	


