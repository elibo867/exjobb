
# coding: utf-8

# In[ ]:


import Get_Densities
import sys


def main():
    density_maps={}
    filenames=sys.argv[1:]
    #for every file create 11 densitymaps. Generates dict of dicts. Kanske inte är smart?
    print 'Computing atom densities'
    counter=1
    for filename in filenames: 
        density_maps[filename]=Get_Densities.main(filename)
        print counter, 'of', len(filenames)
        counter+=1
        
    print 'ready'
        
main()

