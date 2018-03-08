
# coding: utf-8

# In[21]:


import numpy as np
training_array=np.load('training_data_6prot.npz')
testing_array=np.load('test_data_3prot.npz')
y=np.concatenate((training_array['all_scores'],testing_array['all_scores']),axis=0)
print y

