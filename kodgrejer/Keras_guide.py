
# coding: utf-8

# In[43]:


import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.callbacks import ModelCheckpoint
from keras.layers import Conv3D, Activation, MaxPooling3D
from keras.layers import BatchNormalization, Reshape, Dense
from keras import optimizers
from keras import losses



# In[ ]:


def load_data():
    '''

    #THIS IS FOR THE LARGE FILES 
    density_arrays=np.load('../casp9_results/T0515.npz')
    score_arrays=np.load('../casp9_results/T0515_scores.npz')

    l=len(density_arrays.files)


    x=np.zeros((11,120,120,120))
    #concatenate all protein density maps into one array with right shape 
    #first, create array, thereafter add the rest! 
    for i in range(l): 
        if i==0:   
            x=np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)
            print ('x.shape= ',x.shape)
        else: 
            x=np.concatenate((x,np.expand_dims(density_arrays[density_arrays.files[i]], axis=0)))
            if i%5==0:
                print ('x.shape= ',x.shape)
    
            
    y=score_arrays['arr_0']
    
    x=x.astype(float32)
    y=y.astype(float32)
    
    #split the data into training (80%) and testing (20%)  
    x_train,x_test, y_train, y_test=train_test_split(
        x,y,test_size=0.2) #random_state=?? 42?? 
    
    '''
    #load x-training data from npz file.     
    training_array=np.load('training_data_6prot.npz')
    testing_array=np.load('test_data_3prot.npz')
    x=np.concatenate((training_array['all_arrays'], testing_array['all_arrays']), axis=0)
    y=np.concatenate((training_array['all_scores'],testing_array['all_scores']),axis=0)
    
    x=x.astype(np.float32)
    y=y.astype(np.float32)
    
    
    x_train,x_test, y_train, y_test=train_test_split(
    x,y,test_size=0.2) #random_state=?? 42??'''

    
    return x_train, y_train, x_test, y_test


# In[ ]:


def main():
    model=Sequential()
    #the first layer in a Sequential model 
    #(and only the first, because following layers can do automatic shape inference) 
    #needs to receive information about its input shape.
    
    x_train, y_train,x_test,y_test=load_data()

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
    filepath='./output/weights.{epoch:02d}-{loss:.2f}.hdf5'
    checkpoint=ModelCheckpoint(filepath, 
                              monitor='loss',
                              verbose=1,
                              period=5)
    

    #Training. Could save it into a parameter e.g. 'history', for the posibility of plotting metrics later. 
    #Why 52 epochs? Need to avoid overfitting and to run if the model doesn't get better. TRY
    #add validation split to compare models?? 
    model.fit(x_train,y_train, batch_size=9 , epochs=55, callbacks=[checkpoint] )
    
  
    #use a batch as large as you can afford without going out of memory --???
    print(model.evaluate(x_test,y_test))
          
    model.save('test_model.h5', overwrite=True)
    

main()

