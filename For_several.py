
# coding: utf-8

# In[ ]:


import input_mod
import sys

def main():
    density_maps={}
    filenames=sys.argv[1:]
    for filename in filenames: 
        density_maps[filename]=input_mod.main(filename)
        
main()

