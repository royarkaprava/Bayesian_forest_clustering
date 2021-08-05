#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


from sklearn import cluster, datasets, mixture
from sklearn.preprocessing import StandardScaler


# In[2]:


p =200


# In[ ]:





# In[3]:


import sklearn

np.random.seed(9)

# mu = sklearn.datasets.make_blobs(n_samples=p, n_features=2, centers=3,shuffle=False)
mu = sklearn.datasets.make_moons(n_samples=p, noise=0.05,shuffle=False)



Y = mu[0]


plt.scatter(Y[:,0],Y[:,1])


# In[4]:


import numba
from numba import jit



import importlib
import gibbs_sampling as gs


# In[5]:


importlib.reload(gs)


# In[6]:


tree = gs.SpanningTree(Y)


# In[7]:


prob= tree.computeMarginalProb()


# In[8]:


# prob[prob>0.1]=0.1


# In[9]:


np.quantile(tree.D[tree.D>0],0.1)


# In[10]:


trace = tree.runMCMC()


# In[11]:


A1 = gs.getA(trace[0][0])
A2 = gs.getA(trace[24][0])
A3 = gs.getA(trace[99][0])


# In[12]:


from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
rcParams['figure.dpi'] = 300


# In[13]:


import networkx as nx


# In[14]:


width =2


# In[15]:


M= nx.Graph(A1)


# In[16]:


def pltGraph(A, color='r', usingWeight=True):
    M= nx.Graph(A)
    edges = M.edges()
    weights = [prob[u][v] for u,v in edges]
    weights = np.log(weights)
    weights=(weights-np.min(weights))/(np.max(weights)-np.min(weights))*4
    
    if usingWeight:
        nx.draw(M,pos=Y,edge_color=color,width=weights)
    else:
        nx.draw(M,pos=Y,edge_color=color,width=4, node_size=200)


# In[17]:


f = plt.figure()
pltGraph(A1,'red',False)
f.savefig("moons1.png")


# In[18]:


f = plt.figure()
pltGraph(A3,'red',False)
f.savefig("moons2.png")


# In[19]:


prob[prob>0.6]=0.6


# In[ ]:





# In[20]:


import seaborn as sns


# In[21]:


sns.diverging_palette(240, 10, n=9)


# In[22]:


cmap = sns.diverging_palette(240, 10, n=9,as_cmap=True)


# In[23]:


f = plt.figure(figsize=[4,3])

plt.imshow(prob, vmin=0.0,vmax=1,cmap=cmap)
plt.colorbar()

f.savefig("moonsMarginal.png")

