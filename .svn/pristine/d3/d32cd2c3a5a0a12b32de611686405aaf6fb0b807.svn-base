def create_mesh1D(Ly,my,type,alpha,filename='mesh.txt',beta=0.8,gamma=0.05):
	"""
	This function will create 1D mesh according to given shape and type
 	y = create_mesh1D(Ly,my,type,alpha,filename)
		Ex. y = create_mesh1D(172./2,1025,'tanh',0.5,'mesh.txt')
                will create a mesh with more 1/0.5 more resolution close to
 		boundaries than in the rest.
	"""
	import matplotlib.pyplot as plt
	import numpy as np
	from scipy.integrate import cumtrapz
	
	if type == 'tanh':
		y=np.linspace(-1,1,my)
		Dy = np.array(np.ones(my))
		Dy=(1+alpha)+alpha*np.tanh(-(y-0.7)/0.05)
		for j in range(0,my/2):
			Dy[j] =Dy[-1-j] 
		y2 = cumtrapz(Dy,y,initial=0)
		ynew = 2*y2*Ly/y2[-1] - Ly
		np.savetxt(filename, ynew)
	        return ynew 

	elif type =='test':
		y=np.linspace(-1,1,my)
		Dy = np.array(np.ones(my))
		Dy=(1+alpha)+alpha*np.tanh((y-beta)/gamma)
		for j in range(0,my/2):
			Dy[j] =Dy[-1-j] 
		y2 = cumtrapz(Dy,y,initial=0)
		ynew = 2*y2*Ly/y2[-1] - Ly
		np.savetxt(filename, ynew)
	        return ynew 

	elif type =='line':
		ynew=np.linspace(-Ly,Ly,my)
		np.savetxt(filename, ynew)
	        return ynew 




		
	else:
		print "That MESH TYPE is not ready yet MAN..."
		return np.ones(my)
	#plt.plot(ynew[:-1],np.diff(ynew),'b.')	
	#plt.show()


def write4field(filename,time,Re,alp,bet,a0,nx,my,mz,y,VOR,PHI,PSI,SCAL,\
	rhou00,	rhov00,rhow00):
	""" This function writes the typical unformatted input file for LOMA
	none = write4field(filename,input_variables)
	"""
	import numpy as np
	import pylab
	import matplotlib.pyplot as plt
	import sys
        #fmap not used by code
	fmap = np.ones(my)
        nplanes = nx/2+1
	#RECORD 1
	rec1_size = np.dtype('int32')
	rec1_real4 = 5+3*my
	rec1_real8 = 2*my
	rec1_integer = 3
	rec1_size  = rec1_real4*4+rec1_real8*8+rec1_integer*4
        #Y-FMAP and ZERO modes structs
	yfmap = np.zeros((my,),dtype=('f8,f8'))
	uvw00 = np.zeros((my,),dtype=('f8,f8,f8'))
        #FILL y and fmap, and rhou,rhov,rhow	
	for j in range(0,my):
	    yfmap[j]=(y[j],fmap[j])

	for j in range(0,my):
	    uvw00[j]=(rhou00[j],rhov00[j],rhow00[j])
        #Create RECORD1 as is...	
	record1 = np.array([(rec1_size, \
	        time,                 \
	        Re,                   \
	        alp,                  \
	        bet,                  \
	        a0,                   \
	        nx,                   \
	        my,                   \
	        mz,                   \
	        yfmap,                \
	        uvw00,                \
	        rec1_size)],             \
	        dtype=[('dummy1','i4'),  \
	        ('time','f4'),\
	        ('Re','f4'),  \
	        ('alp','f4'), \
	        ('bet','f4'), \
	        ('a0','f4'), \
	        ('mx','i4'),\
	        ('my','i4'),\
	        ('mz','i4'),\
	        ('yfmap','f8,f8',my),\
	        ('uvw00','f4,f4,f4',my),\
	        ('dummy2','i4')])
	          
	#OPEN FILE
	output_file = open(filename,'wb')
	#WRITE RECORD 1 TO FILE
	record1.tofile(output_file)
	#CREATE RECORD 2 FOR EACH PLANE
	ntotr=8*my*mz #size of each plane  
	rec2_size = ntotr*4 
	#..................................................#
	#           Writing plane by plane                 #
	#Prepare buffers to avoid ORDER errors with ARRAYS
	temp1 = np.zeros(ntotr,dtype='f4')
	temp2 = np.reshape(temp1,(4,2*my,mz),order='F')
	ireal = range(0,2*my,2)
	iimag = range(1,2*my,2)
	for i in range(nplanes):
		for k in range(mz):
			temp2[0,ireal,k] = VOR[:,k,i].real
	        	temp2[0,iimag,k] = VOR[:,k,i].imag
	        	temp2[1,ireal,k] = PHI[:,k,i].real
	        	temp2[1,iimag,k] = PHI[:,k,i].imag
	        	temp2[2,ireal,k] = PSI[:,k,i].real
	        	temp2[2,iimag,k] = PSI[:,k,i].imag
	        	temp2[3,ireal,k] = SCAL[:,k,i].real
	        	temp2[3,iimag,k] = SCAL[:,k,i].imag
	        	temp1 = np.reshape(temp2,8*my*mz,order='F')
                	#struct for PLANE RECORD
	        	record2 = np.array([(rec2_size, \
	                        	temp1[:],\
	                             	rec2_size)], \
					dtype=[('dummy1','i4'),\
	                              ('recplane','f4',8*my*mz),\
	                              ('dummy2','i4')])
		#Writing plane...
		record2.tofile(output_file)
        #Close FILE	
	output_file.close()

	

