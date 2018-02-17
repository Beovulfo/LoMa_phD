#!/usr/local/bin/python2.7
import numpy as np
from pyvtk import *
def pln2vtk(x,y,z,field,fieldname,fname,ASCII=1):
    """Converts a plane from lomaHZ to .vtk ASCII"""
    mgalx = len(x);my = len(y); mgalz=len(z);
    if mgalz==1:
        npoints = my*mgalx
        mydata=np.float32(np.reshape(np.transpose(field),[npoints]))
    elif mgalx==1:
        npoints =my*mgalz
        mydata=np.float32(np.reshape(np.transpose(field),[npoints]))
    else:
        npoints = my*mgalx*mgalz
        mydata=np.float32(np.reshape(np.transpose(field),[npoints]))
#mydata = np.array(mydata)
#print np.shape(mydata)
    y2 = np.float32(y)
    x2 = np.float32(x)
    z2 = np.float32(z)
    vtk = VtkData(RectilinearGrid(list(x2),list(y2),list(z2)),
             PointData(Scalars(mydata,name=fieldname)))
    if ASCII==1:
        vtk.tofile(fname)
    else:
        vtk.tofile(fname,'binary')
 


def create_vtk_vector(x,y,z,field,fieldname,filename,ASCII=0):
    """
    Used to create vectors such velocity or vorticity
    previously the components of the vector should have been loaded as:
    varpln1='upyx';varpln2='vpyx';varpln3='wpyx'
    field = [ [] for i in range(3)];

    [field[0],y,x,z] = readfieldij(filename+'.' + varpln1);
    [field[1],y,x,z] = readfieldij(filename+'.'+varpln2);
    [field[2],y,x,z] = readfieldij(filename+'.'+varpln3);


    
    """
    fname = filename
    data = [ 1 for i in range(3)];
    mgalx = len(x);my = len(y); mgalz=len(z);
    #fieldname = 'vel'
    if mgalz==1:
        npoints = my*mgalx
        data = np.float32(np.reshape(np.transpose(field),[npoints*3]))
    elif mgalx==1:
	npoints = my*mgalz
        data = np.float32(np.reshape(np.transpose(field),[npoints*3]))
    else:
	npoints = my*mgalx*mgalz
        data = np.float32(np.reshape(np.transpose(field),[npoints*3]))

#mydata = np.array(mydata)
#print np.shape(mydata)
    y2 = np.float32(y)
    x2 = np.float32(x)
    z2 = np.float32(z)
    vtk = VtkData(RectilinearGrid(list(x2),list(y2),list(z2)),
             PointData(Vectors(data,name=fieldname)))
    if ASCII==1:
        vtk.tofile(fname)
    else:
        vtk.tofile(fname,'binary')



def Two2vtk(x,y,z,field1,field2,fieldname1,fieldname2,fname,ASCII=0):
    """Converts a plane from lomaHZ to .vtk ASCII"""
    mgalx = len(x);my = len(y); mgalz=len(z);
    if mgalz==1:
        npoints = my*mgalx
        mydata=np.float32(np.reshape(np.transpose(field1),[npoints]))
        mydata2=np.float32(np.reshape(np.transpose(field2),[npoints]))
    elif mgalx==1:
        npoints =my*mgalz
        mydata=np.float32(np.reshape(np.transpose(field1),[npoints]))
        mydata2=np.float32(np.reshape(np.transpose(field2),[npoints]))
    else:
        npoints = my*mgalx*mgalz
        mydata=np.reshape(np.transpose(field1),[npoints])
        del field1
        mydata=np.float32(mydata)
        #mydata=np.float32(np.reshape(np.transpose(field1),[npoints]))
#mydata = np.array(mydata)
#print np.shape(mydata)
    y2 = np.float32(y)
    x2 = np.float32(x)
    z2 = np.float32(z)
    vtk = VtkData(RectilinearGrid(list(x2),list(y2),list(z2)),
             PointData(Scalars(mydata,name=fieldname1)))
    del mydata
    mydata=np.reshape(np.transpose(field2),[npoints])
    mydata=np.float32(mydata)
    vtk.point_data.append(Scalars(np.float32(mydata),name=fieldname2))
    if ASCII==1:
        vtk.tofile(fname)
    else:
        vtk.tofile(fname,'binary')
 



