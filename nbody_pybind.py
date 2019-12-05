#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('pylab', '')
import matplotlib.animation as animation
plt.ion()


# In[2]:


from build_release import nbody
import importlib
importlib.reload(nbody)
nbody.init_threads()


# In[13]:


npart = 500
xs = np.random.uniform(-20,20,size=(npart,2))
#vs = np.random.normal(0,5,(npart,2))
vs = np.zeros((npart+1,2))
ms = np.random.exponential(3,npart).reshape(-1,1)

xs_orig,vs_orig,ms_orig = xs.copy(),vs.copy(),ms.copy()


# In[6]:


nbody.step(xs,vs,ms,2)
plt.scatter(xs[:,0], xs[:,1], s=4*ms)


# In[4]:





# In[16]:


from IPython.display import HTML
get_ipython().run_line_magic('matplotlib', 'qt')
def animate(i):
    global vs
    nbody.step(xs,vs,ms,50)
    #if i%100 == 0:    
    #vs = -vs
    return datas,


fig = plt.figure()
ax = plt.axes()
ax.set_xlim(-7, 7)
ax.set_ylim(-7, 7)
datas = ax.scatter([], [], [])
datas.set_offsets(xs)
datas.set_sizes(ms.ravel() * 3)

ani = animation.FuncAnimation(fig, animate, frames=range(1000), interval=1000/25, blit = False)
#HTML(ani.to_jshtml())


# In[5]:


npart = 100
xs = np.random.uniform(-0.5,0.5,size=(npart,2))
#vs = np.random.normal(0,5,(npart,2))
vs = np.zeros((npart,2))
ms = np.random.exponential(3,npart).reshape(-1,1)
get_ipython().run_line_magic('timeit', 'nbody.step(xs,vs,ms,100)')


# In[8]:


xs


# In[ ]:




