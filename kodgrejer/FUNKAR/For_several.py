
# coding: utf-8

# In[ ]:


import Get_Densities
import sys
import numpy as np


# In[ ]:


def make_4Darray(filename):
    '''Extract density maps from dictionary to np.array'''
    #get the dictionary with atom densities
    density_dict=Get_Densities.main(filename)
    
    #Concatenate the maps into one single array
    for key in sorted(density_dict): 
        #The first time, need to  create the new array
        if key==1: 
            density_array=np.expand_dims((density_dict[key]),axis=0)
        else:
            density_array=np.concatenate((density_array, np.expand_dims((density_dict[key]),axis=0)), axis=0)
            
    return density_array

   
        
  


# In[ ]:


def main():
    
    filenames=sys.argv[1:]
    #for every file create 11 densitymaps. Generates dict of dicts. Kanske inte är smart?
    print ('Computing atom densities')
    counter=1
    #only making one for starters, the last file will be stored in array
    for filename in filenames: 
        dens_array=make_4Darray(filename)
    print (dens_array.shape)
   # np.save('outarray', dens_array)
    print ('ready')
    
        
main()

