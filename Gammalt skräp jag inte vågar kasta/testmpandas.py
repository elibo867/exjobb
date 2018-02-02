
# coding: utf-8

# In[108]:


from biopandas.pdb import PandasPdb
from math import exp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import spatial
from mpl_toolkits.mplot3d import Axes3D
#get_ipython().magic(u'matplotlib inline')


# In[109]:


def centre_model(atoms):
    '''Calculate the centre of gravity and "moves" the protein so the CG is located in origo. 
    The decimals of the coordinates are wrong (too many). 
    How to solve? '''
    
    CG=[atoms.x_coord.mean(),atoms.y_coord.mean(), atoms.z_coord.mean()]     
    atoms['x_coord']=atoms.x_coord.subtract(CG[0]) #subtracts CG coordinate from all coordinates
    atoms['y_coord']=atoms.y_coord.subtract(CG[1])
    atoms['z_coord']=atoms.z_coord.subtract(CG[2])



# In[110]:


def create_density_map(atoms):
    
    dens_map=np.zeros((120,120,120)) #120x120x120 matrix filled w/ zeroes

    x_s=atoms.x_coord.subtract(-60) #adding 60 (subtracting -60, how to du addition???????)
    y_s= atoms.y_coord.subtract(-60)
    z_s=atoms.z_coord.subtract(-60)
    
    coords=[x_s, y_s, z_s]
    coords=np.array(coords).T.tolist() #transpose the coordinate list
 
    for i in range(len(coords)):
        x=int(round(coords[i][0]))
        y=int(round(coords[i][1]))
        z=int(round(coords[i][2]))
        
        for  k in range(x-2,x+3): #need to take +3 to make it get to +2
            for l in range(y-2,y+3):
                for m in range(z-2,z+3): 
                    
                    if k==x and l==y and m==z: #if xyz = klm , r=0 
                        dens_map[k][l][m]+=exp(0)

                    elif ((x-1)<= k <=(x+1)) and ((y-1)<=l<=(y+1)) and ((z-1)<=m<=(z+1)):       
                        dens_map[k][l][m]+=exp(1/2)
                    else:
                        dens_map[k][l][m]+=exp(-2)
  

    return dens_map


# In[111]:


def plot_map(dens_map):

    values=[]
    for i in range(120):
        for j in range(120):
            for k in range(120):
                if dens_map[i][j][k]!=0:
                    values.append([i , j , k , dens_map[i][j][k]])

    values=np.transpose(values)
    x=values[0]
    y=values[1]
    z=values[2]
    density=values[3]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.scatter(x, y, z, c=density, cmap=cm.coolwarm)
    plt.show()
  
  


# In[112]:


def main():
    ppdb=PandasPdb()
    ppdb=ppdb.read_pdb('5eh6.pdb')#('HHpredAQ_TS1.pdb')
    atoms=ppdb.df['ATOM']
    centre_model(atoms)

    density_map=create_density_map(atoms)
    
    plot_map(density_map)
    
    

   

main()


# In[ ]:



    

