
# coding: utf-8

# In[ ]:


import Get_Densities
import Get_GDT_TS
import sys
import os
import re
import numpy as np
import zipfile
from tempfile import TemporaryFile


# In[ ]:


def get_target(filename):
    f=open(filename)
    f.close()
    f_path=os.path.realpath(f.name)
    
    #Looking at the path to see which target that's beeing processed
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
    
    target_name=get_target(filenames[0])
    zip_name=target_name+'.npz'
    
    #better to check if the filename exists and then add something to the name . instead of overwriting
    try:
        os.remove(zip_name)
    except FileNotFoundError:
        pass
    
    with zipfile.ZipFile(zip_name, mode='a', compression=zipfile.ZIP_DEFLATED) as zf:
        for filename in filenames: 

            print (counter,'of', len(filenames), '--', filename)    

            #Try collect density maps and GDT score. If it doesn't work - 
            try: 
                #compute the 11 density maps
                dens_array=Get_Densities.main(filename)
                #create list of arrays: [all_maps_prot1, all_maps_prot2,...] 
                #all_arrays.append(dens_array)
            
                #find GDT_TS score - if it exists. Otherwise, returns 'No GDT_TS-score in TM file' 
                filename=filename+('.fixed.TM')
                GDT=Get_GDT_TS.main(filename)
            
                #If there is no GDT score in the TM file 
                if isinstance(GDT, str):
                    print (GDT)
                    print(filename)
                    no_passed+=1
                    counter+=1
                    continue
                    
                tmpfilename='{0}arr_{1}.npy'.format(target_name+'_', counter-1-no_passed)
                np.save(tmpfilename, dens_array)
                zf.write(tmpfilename)
                
                os.remove(tmpfilename)
                
                all_scores.append(GDT)
                counter+=1
        
            except Exception as e:
                print (e)
                print (filename)
                no_passed+=1
                counter+=1
                continue
    
    #if counter!=2 and (counter+1)%5!=0:
        #print (counter-1,'of', len(filenames))
    
    
    #generates zip_file with one array shape (x,11,120,120,120) where x is number of proteins used
    np.savez_compressed(target_name+'_scores', all_scores=all_scores)
    print (no_passed, 'files ignored')

        
main()

