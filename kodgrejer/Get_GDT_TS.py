
# coding: utf-8

# In[5]:


import sys 
import re
import numpy as np


# In[6]:


def extract_score(filename):
    f=open(filename,'r')
    text=f.read()
    
    #get GDT_TS_score
    GDT_TS_match=re.search('GDT-TS-score=\s(\d.\d\d\d\d)',text)
    try:
        score=float(GDT_TS_match.group(1))
        return score
    
    except Exception as e:
        print (e)
        return 'No GDT_TS-score in TM file'
        
        
     


# In[7]:


def main(argv): 

    filename=argv
    scores=extract_score(filename)
    
           
    #scores_array=np.asarray(scores) 
    
    #np.savez_compressed('y_train_3', scores)
    
    return scores
    
if __name__ == '__main__':
    main(sys.argv) 
        

