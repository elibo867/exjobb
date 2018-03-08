
# coding: utf-8

# In[9]:


import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.models import load_model 
from keras.layers import Conv3D, Activation, MaxPooling3D
from keras.layers import BatchNormalization, Reshape, Dense
from keras import optimizers
from keras import losses


# In[14]:



density_arrays=np.load('../kodgrejer/T0518.npz')
score_arrays=np.load('../kodgrejer/T0518_scores.npz')
print density_arrays.files
print score_arrays.files


l=4

x=np.zeros((11,120,120,120))
#concatenate all protein density maps into one array with right shape 
#first, create array, thereafter add the rest! 
for i in range(l): 
    if i==0:   
        x=np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)
        print x.shape
    else: 
        x=np.concatenate((x,np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)))
        print x.shape
            
y=score_arrays['all_scores'][:4]
   
#split the data into training (80%) and testing (20%)  
x_train,x_test, y_train, y_test=train_test_split(
        x,y,test_size=0.2) #random_state=?? 42?? 

model=load_model('test_model.h5')
model.fit(x_train,y_train, batch_size=9 , epochs=3, shuffle=True, )

