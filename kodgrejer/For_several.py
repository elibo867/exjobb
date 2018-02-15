
# coding: utf-8

# In[ ]:


import Get_Densities
import Get_GDT_TS
import sys
import os
import re
import numpy as np


# In[ ]:


def make_4Darray(filename):
    '''Extract density maps from dictionary to np.array'''
    #get the dictionary with atom densities - why not do an array immediately?
    density_dict=Get_Densities.main(filename)
    
    #Concatenate the 11 maps into one single array
    for key in sorted(density_dict): 
        #The first time, need to  create the new array
        if key==1: 
            density_array=np.expand_dims((density_dict[key]),axis=0)
        else:
            density_array=np.concatenate((density_array, np.expand_dims((density_dict[key]),axis=0)), axis=0)
            
    return density_array


# In[ ]:


def get_target(filename):
    f=open(filename)
    f.close()
    f_path=os.path.realpath(f.name)
    
    target_match=re.search('(T\d\d\d\d)',f_path)
    
    if target_match:
        target_name=target_match.group(1)
    
        return target_name 
    
    


# In[ ]:


def main():
    
    filenames=sys.argv[1:]
    #for every file create 11 densitymaps. 
    print ('Computing atom densities')
    counter=1
    no_passed=0
    all_arrays=[]
    all_scores=[]
    
    
    print (counter, 'of', len(filenames))
    
    for filename in filenames: 
        
        if counter%25==0:
            print (counter,'of', len(filenames))
        
        #Try to create densitymaps and collect GDT score
        try: 
            #compute the 11 density maps
            dens_array=make_4Darray(filename)
            #create list of arrays: [all_maps_prot1, all_maps_prot2,...] 
            all_arrays.append(dens_array)
            
            #find GDT_TS score
            filename=filename+('.fixed.TM')
            GDT=Get_GDT_TS.main(filename)
            all_scores.append(GDT)
            
            counter+=1
        
        except IOError: 
            print ('cannot find', filename)
            no_passed+=1
            counter+=1
            continue
    
    print (counter-1,'of', len(filenames))
    
    target_name=get_target(filenames[0])
    
    
    

    #generates zip_file with one array shape (x,11,120,120,120) where x is number of proteins used
    np.savez_compressed(target_name, all_arrays=all_arrays, all_scores=all_scores)
    print (no_passed, 'files ignored')

        
main()

