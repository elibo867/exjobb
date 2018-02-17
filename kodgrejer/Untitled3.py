
# coding: utf-8

# In[46]:


import numpy as np
density_arrays=np.load('T0518.npz')
#score_arrays=np.load('T0518_scores.npz')

l=len(density_arrays.files)
x=np.zeros((l,11,120,120,120))


for i in range(l): 
    x[i, :, : , :]=density_arrays[density_arrays.files[i]]

print np.max(x)


# In[44]:


len(density_arrays.files)

