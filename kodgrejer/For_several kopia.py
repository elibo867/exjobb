
# coding: utf-8

# In[ ]:


import Get_Densities
import sys

def input_to_4Darray(filename):
    density_maps[filename]=Get_Densities.main(filename)
        
    for key in density_maps[filename].keys(): # key= atomtyp int 1-11
        print key
    print counter, 'of', len(filenames)
    counter+=1

def main():
    density_maps={}
    filenames=sys.argv[1:]
    #for every file create 11 densitymaps. Generates dict of dicts. Kanske inte är smart?
    print 'Computing atom densities'
    counter=1
    
    for filename in filenames: 
        input_to_4Darray(filename)
        
        
    print 'ready'
        
main()

