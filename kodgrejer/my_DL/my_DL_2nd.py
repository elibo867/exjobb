
# coding: utf-8

# In[9]:

import sys
import os
import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.models import load_model 
from keras.callbacks import ModelCheckpoint

from keras.layers import Conv3D, Activation, MaxPooling3D
from keras.layers import BatchNormalization, Reshape, Dense
from keras import optimizers
from keras import losses

 
def load_data(target_name):
	#if there isn't a file with a full target array
	if not os.path.isfile('../../disk/casp9_results/{0}_array.npz'.format(target_name)):
		#assemble an array and save it 
		density_arrays=np.load('../../disk/casp9_results/{0}.npz'.format(target_name))
		
		l=len(density_arrays.files)


		x=np.ndarray(shape=(l,11,120,120,120), dtype=np.float32)
		#concatenate all protein density maps into one array with right shape 
		#first, create array, thereafter add the rest! 
		print ('Loading...')
		for i in range(l): 
			x[i,:,:,:,:]=density_arrays[density_arrays.files[i]].astype(np.float32)
			if (i+1)%5==0:
				print (i+1, 'of', l)

		print ('saving...')
		
		#save array into .npz file 
		np.savez_compressed('../../disk/casp9_results/{0}_array.npz'.format(target_name), x=x)
		print ('{0}-- klar'.format(target_name))
		
	else: 
		density_arrays=np.load('../../disk/casp9_results/{0}_array.npz'.format(target_name))
		x=density_arrays['x']
		print ('{0}-- klar'.format(target_name))
		
	#load GDT-TS score array 
		
	score_arrays=np.load('../../disk/casp9_results/{0}_scores.npz'.format(target_name))
	y=score_arrays['all_scores']
   

	#split the data into training (80%) and testing (20%)  
	#print ('borjar splitta')
	#x_train,x_test, y_train, y_test=train_test_split(
		#x,y,test_size=0.2, shuffle=False) #random_state=?? 42?? 
	#print ('splittat klart')
	return x, y#, x_test, y_test


def main():
	target_names=sys.argv[1:]

	print ('loading', target_names[0])    
	x_train,y_train=load_data(target_names[0])	
	
	model=load_model('trained_model.h5')

	filepath='./output/weights.{epoch:02d}-{loss:.4f}.hdf5'
	checkpoint=ModelCheckpoint(filepath, 
							  monitor='loss',
							  verbose=1,
							  period=5,
							  save_best_only=True,
							  mode='min')

	model.fit(x_train,y_train, batch_size=9 , epochs=25, callbacks=[checkpoint])

	model.save('trained_model.h5', overwrite=True)

main()