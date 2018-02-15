
# coding: utf-8

# In[2]:


import numpy as np

all_arrays=np.load('x_train_6.npz')
x_train=all_arrays['arr_0']
print x_train.shape

