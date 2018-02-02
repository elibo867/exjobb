
# coding: utf-8

# In[122]:


from biopandas.pdb import PandasPdb
from math import exp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import spatial
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib inline  


# In[123]:


def centre_model(atoms):
    '''Calculate the centre of gravity and "moves" the protein so the CG is located in origo. 
    The decimals of the coordinates are wrong (too many). 
    How to solve? '''
    
    CG=[atoms.x_coord.mean(),atoms.y_coord.mean(), atoms.z_coord.mean()]     
    atoms['x_coord']=atoms.x_coord.subtract(CG[0]) #subtracts CG coordinate from all coordinates
    atoms['y_coord']=atoms.y_coord.subtract(CG[1])
    atoms['z_coord']=atoms.z_coord.subtract(CG[2])



# In[124]:


def atoms_to_map(atoms):
    
    
    
    atom_types_dict={('ALA', 'CB') : 'carbon_sp3',
                ('ARG', 'CB') : 'carbon_sp3',
                ('ARG', 'CG') : 'carbon_sp3',
                ('ARG', 'CD') : 'carbon_sp3',
                ('ASN', 'CB'): 'carbon_sp3' ,
                ( 'ASP', 'CB') : 'carbon_sp3',
                ('CYS', 'CB') : 'carbon_sp3',
                ('GLN', 'CB') : 'carbon_sp3',
                ('GLN', 'CG') : 'carbon_sp3',
                ('GLU', 'CB') : 'carbon_sp3',
                ('GLU', 'CG') : 'carbon_sp3' ,
                ('HIS', 'CB') : 'carbon_sp3' ,
                ('ILE' , 'CB') : 'carbon_sp3',
                ('ILE' , 'CG1') : 'carbon_sp3',
                ('ILE' , 'CG2') : 'carbon_sp3',
                ('ILE' , 'CG3') : 'carbon_sp3',
                ('ILE' , 'CD1') : 'carbon_sp3',
                ('LEU' , 'CB') : 'carbon_sp3',
                ('LEU' , 'CG') : 'carbon_sp3',
                ('LEU' , 'CD1') : 'carbon_sp3',
                ('LEU' , 'CD2') : 'carbon_sp3',
                ('LEU' , 'CD3') : 'carbon_sp3',
                ('LYS' , 'CB') : 'carbon_sp3',
                ('LYS' , 'CG') : 'carbon_sp3',
                ('LYS' , 'CD') : 'carbon_sp3',
                ('LYS' , 'CE') : 'carbon_sp3',
                ('MET' , 'CB') : 'carbon_sp3',
                ('MET' , 'CG') : 'carbon_sp3',
                ('MET' , 'CE') : 'carbon_sp3',
                ('MSE' , 'CB') : 'carbon_sp3',
                ('MSE' , 'CG') : 'carbon_sp3',
                ('MSE' , 'CE') : 'carbon_sp3',
                ('PHE' , 'CB') : 'carbon_sp3',
                ('PRO' , 'CB') : 'carbon_sp3',
                ('PRO' , 'CG') : 'carbon_sp3',
                ('PRO' , 'CD') : 'carbon_sp3',
                ('SER' , 'CB') : 'carbon_sp3',
                ('THR' , 'CB') : 'carbon_sp3',
                ('THR' , 'CG2') : 'carbon_sp3',
                ('TRP' , 'CB') : 'carbon_sp3',
                ('TYR' , 'CB') : 'carbon_sp3',
                ('VAL' , 'CB') : 'carbon_sp3', 
                ('VAL' , 'CG1') : 'carbon_sp3',
                ('VAL' , 'CG2') : 'carbon_sp3',
                ('VAL' , 'CG3') : 'carbon_sp3',
                #and backbone CA - but don't know how
                }
    indices=[]
    for index,row in atoms.iterrows(): 
        if (row['residue_name'], row['atom_name']) not in atom_types_dict:
            indices.append(index)

    #drop?????????????
            


# In[ ]:






# In[125]:


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


# In[126]:


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
  
  


# In[127]:


def main():
    ppdb=PandasPdb()
    ppdb=ppdb.read_pdb('5eh6.pdb')#('HHpredAQ_TS1.pdb')
    atoms=ppdb.df['ATOM']
    centre_model(atoms)
    
    new_atoms=atoms_to_map(atoms)
    
    
    
    
    
    density_map=create_density_map(atoms)
    
    plot_map(density_map)
    
    

   
main()


# In[ ]:



    

