
# coding: utf-8

# In[141]:


import numpy as np

density_arrays=np.load('T0515_TEST.npz')
score_arrays=np.load('T0515_scores.npz')

l=len(density_arrays.files)
#print density_arrays['arr_0']
#print density_arrays.shape


#print density_arrays['arr_0'].item()[1].shape

x=np.zeros((11,120,120,120))



for i in range(l): 
    if i==0:   
        x=np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)
        print x.shape
    else: 
        x=np.concatenate((x,np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)))
        print x.shape


# In[3]:




