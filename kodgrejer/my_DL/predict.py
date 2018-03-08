
# coding: utf-8

# In[9]:


import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.models import load_model 
from keras.callbacks import ModelCheckpoint

from keras.layers import Conv3D, Activation, MaxPooling3D
from keras.layers import BatchNormalization, Reshape, Dense
from keras import optimizers
from keras import losses


def main():
	print('loading x')
	T0518_array=np.load('../../disk/casp9_results/T0518_array.npz')
	x=T0518_array['x'][:9]
	print ('loading model')
	model=load_model('trained_model.h5')

	print('prediction:')
	predictions=model.predict(x, batch_size=None, verbose=1)
	print (predictions)

	print ('loading true')
	score_array=np.load('../../disk/casp9_results/T0518_scores.npz')
	y=score_array['all_scores'][:9]
	print ('true')
	print (y)
main()