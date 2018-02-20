
# coding: utf-8

# In[2]:


import numpy as np

density_arrays=np.load('../../../proj/x_elibo/generate_dataset/casp9/casp9_results/T0515.npz')
score_arrays=np.load('../../../proj/x_elibo/generate_dataset/casp9/casp9_results/T0515_scores.npz')

l=len(density_arrays.files)


x=np.zeros((11,120,120,120))
#concatenate all protein density maps into one array with right shape 
for i in range(l): 
    if i==0:   
        x=np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)
        print x.shape
    else: 
        x=np.concatenate((x,np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)))
        print x.shape


# In[3]:




