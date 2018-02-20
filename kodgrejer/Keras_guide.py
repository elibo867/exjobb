
# coding: utf-8

# In[ ]:


import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Conv3D, Activation, MaxPooling3D
from keras.layers import BatchNormalization, Reshape, Dense
from keras import optimizers
from keras import losses



# In[ ]:


def load_data(): 
    
    density_arrays=np.load('../../../proj/x_elibo/generate_dataset/casp9/casp9_results/T0515.npz')
    score_arrays=np.load('../../../proj/x_elibo/generate_dataset/casp9/casp9_results/T0515_scores.npz')

    l=len(density_arrays.files)


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
            
    y=score_arrays['arr_0']
    
    x_train,x_test, y_train, y_test=train_test_split(
        x,y,test_size=0.33) #random_state=?? 42?? 
    
    '''
    
    #load x-training data from npz file. 
    #First argument in x_training array is the data,(X,11,120,120,120)
    x_training_array=np.load('x_train_6.npz')
    x_train=x_training_array['arr_0']
  
    #load x-test data
    x_testing_array=np.load('x_test_3.npz')
    x_test=x_testing_array['arr_0']
    
    #load y-train data
    y_training_array=np.load('y_train_6.npz')
    y_train=y_training_array['arr_0']

    #load y-test data
    y_testing_array=np.load('y_test_3.npz')
    y_test=y_testing_array['arr_0']
    
    training_array=np.load('training_data_6prot.npz')
    x_train=training_array['all_arrays']
    y_train=training_array['all_scores']
    
    testing_array=np.load('test_data_3prot.npz')
    x_test=testing_array['all_arrays']
    y_test=testing_array['all_scores']'''  
    
    
    
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
  

    #Compiler! The optimizer is correct but the loss is just for testing.
    adam=optimizers.Adam(lr=0.0003, decay=0.01)
    model.compile(optimizer=adam, loss='mean_squared_error', metrics=['accuracy'])
    

    #Training
    model.fit(x_train,y_train, epochs=6)
    
  
    
    model.evaluate(x_test,y_test)
    

main()

