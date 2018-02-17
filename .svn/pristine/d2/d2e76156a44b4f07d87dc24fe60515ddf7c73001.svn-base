import numpy as np
import mayavi.mlab as mlab
from readfiles import *
filename1='turb/3DincRe160_025.o1xz'
filename2='turb/3DincRe160_025.o2xz'
filename3='turb/3DincRe160_025.o3xz'
(y,jspecy,field1,x,z) = readfieldxz(filename1)
(y,jspecy,field2,x,z) = readfieldxz(filename2)
(y,jspecy,field3,x,z) = readfieldxz(filename3)
[xx1,yy1,zz1] = np.meshgrid(x[0:200],y,z[0:200],indexing='ij')
s1 = field1[0:200,:,0:200]
s2 = field2[0:200,:,0:200]
s3 = field3[0:200,:,0:200]
s = 0.5*(s1**2+s2**2+s3**2)
mlab.figure()
mlab.contour3d(xx1,yy1,zz1,s,contours=10,opacity=0.3)
mlab.show()


