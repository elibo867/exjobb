
# coding: utf-8

# In[ ]:


import sys 
import os 
import numpy as np
#from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.callbacks import ModelCheckpoint
from keras.layers import Conv3D, Activation, MaxPooling3D
from keras.layers import BatchNormalization, Reshape, Dense
from keras import optimizers
from keras import losses


# In[ ]:


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


# In[ ]:


def main():
    target_names=sys.argv[1:]
    print ('loading', target_names[0])    
    x_train,y_train=load_data(target_names[0])

    model=Sequential()
    #the first layer in a Sequential model 
    #(and only the first, because following layers can do automatic shape inference) 
    #needs to receive information about its input shape.

    #for i in range target_names:

  
    #x_train, y_train,x_test,y_test=load_data(target_name)

    filter_size=(3,3,3) 
    
    #The network structure such as described in 3DCNN article

    #1
    model.add(Conv3D(16, filter_size, strides=(1,1,1), padding='valid', data_format='channels_first',input_shape=(11,120,120,120)))
    model.add(Activation('relu'))
    model.add(MaxPooling3D(pool_size=(3,3,3), strides=(2,2,2), data_format='channels_first'))
    #4
    model.add(Conv3D(32,filter_size,strides=1,padding='valid', data_format='channels_first'))
    model.add(BatchNormalization(axis=1, scale=False))
    model.add(Activation('relu'))
    model.add(MaxPooling3D(pool_size=(3,3,3), strides=(2,2,2),data_format='channels_first'))
    #8
    model.add(Conv3D(32,filter_size, strides=1, padding='valid',data_format='channels_first'))
    model.add(BatchNormalization(axis=1, scale=False))
    model.add(Activation('relu'))
    #11
    model.add(Conv3D(64,filter_size, strides=1, padding='valid',data_format='channels_first'))
    model.add(BatchNormalization(axis=1, scale=False))
    model.add(Activation('relu'))
    model.add(MaxPooling3D(pool_size=(3,3,3), strides=(2,2,2),data_format='channels_first'))
    
    #15
    model.add(Conv3D(128,filter_size, strides=1,padding='valid', data_format='channels_first'))
    model.add(BatchNormalization(axis=1, scale=False))
    model.add(Activation('relu'))
    model.add(Conv3D(128,filter_size, strides=1, padding='valid',data_format='channels_first'))
    model.add(BatchNormalization(axis=1, scale=False))
    model.add(Activation('relu'))
    model.add(Conv3D(256,filter_size, strides=1,padding='valid', data_format='channels_first'))
    model.add(BatchNormalization(axis=1, scale=False))
    model.add(Activation('relu'))
    model.add(Conv3D(512,filter_size, strides=1, padding='valid',data_format='channels_first'))
    model.add(BatchNormalization(axis=1, scale=False))
    model.add(Activation('relu'))
    model.add(MaxPooling3D(pool_size=(3,3,3), strides=(2,2,2),data_format='channels_first'))
    
    #28
    model.add(Reshape((512,)))
    
    #29
    model.add(Dense(256, activation='linear'))
    model.add(Activation('relu'))
    model.add(Dense(128, activation='linear'))
    model.add(Activation('relu'))
    model.add(Dense(1, activation='linear'))
    
    #print all the layers and output dimensions 
    #model.summary()
  

    #Compiler. There are more arguments that could be added. 
    adam=optimizers.Adam(lr=0.0003, decay=0.01)
    model.compile(optimizer=adam, loss='mse', metrics=['mae'])

    
    #saving checkpoints 
    filepath='./output/weights.{epoch:02d}-{loss:.4f}.hdf5'
    checkpoint=ModelCheckpoint(filepath, 
                              monitor='loss',
                              verbose=1,
                              period=5,
                              save_best_only=True,
                              mode='min')
    

    #Training. Could save it into a parameter e.g. 'history', for the posibility of plotting metrics later. 
    #Why 52 epochs? Need to avoid overfitting and to run if the model doesn't get better. TRY
    #add validation split to compare models?? 
    model.fit(x_train,y_train, batch_size=9 , epochs=25, callbacks=[checkpoint] )

    model.save('trained_model.h5', overwrite=True)
    
  
    #use a batch as large as you can afford without going out of memory --???
    #print(model.evaluate(x_test,y_test, batch_size=8, verbose=1))'''

    

main()

