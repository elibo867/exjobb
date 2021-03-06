
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
    score=float(GDT_TS_match.group(1))
    
    return score 


# In[7]:


def main(argv): 
    #filenames=sys.argv[1:]
    filename=argv
    scores=[]
    #for filename in filenames: 
        #scores.append(extract_score(filename))
    scores=extract_score(filename)
    

    
           
    #scores_array=np.asarray(scores) 
    
    #np.savez_compressed('y_train_3', scores)
    print ('ready')
    
    return scores
    
if __name__ == '__main__':
    main(sys.argv) 
        

