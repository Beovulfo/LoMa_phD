
# coding: utf-8

# In[1]:

import pylab
import matplotlib.pyplot as plt
import numpy as np
#get_ipython().magic(u'matplotlib inline')


# In[2]:

from writefiles import *


# In[3]:

#y = create_mesh1D(258./2,1025,'tanh',0.5,'mesh.txt')
y = create_mesh1D(20,257,'line',0.5,'1Dmesh.txt')


# In[ ]:

plt.plot(y[0:-1],np.diff(y),'b*')
print "Min Dy = %s; Max Dy = %s"%(np.min(np.diff(y)),np.max(np.diff(y)))
plt.show()


