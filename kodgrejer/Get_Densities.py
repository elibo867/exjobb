
# coding: utf-8

# In[23]:


import numpy as np
import matplotlib.pyplot as plt
import sys

from biopandas.pdb import PandasPdb
from math import exp
from matplotlib import cm
#from scipy import spatial
from mpl_toolkits.mplot3d import Axes3D
#remember to comment following row if running code in terminal 
#%matplotlib inline  


# In[24]:


def centre_model(atoms):
    '''Calculate the centre of gravity and "moves" the protein so the CG is located in origo.'''
    
    #Find centre of gravity and subtract from every coordinate 
    CG=[atoms.x_coord.mean(),atoms.y_coord.mean(), atoms.z_coord.mean()]     
    atoms['x_coord']=atoms.x_coord.subtract(CG[0])
    atoms['y_coord']=atoms.y_coord.subtract(CG[1])
    atoms['z_coord']=atoms.z_coord.subtract(CG[2])
    
    return atoms


# In[25]:


def atoms_to_map(atoms, atom_type):  
    '''Selects the rows of the atoms df of the given atom type, returns them in a df'''

    #Create a dictionary with all different atom types. (residue name, atom name) or (atom name) as key and atom type as value
    
    atom_types_dict={ 
                    #Type1 - Sulfur/ selenium
                    ('CYS', 'SG') : 1, ('MET', 'SD') : 1, ('MSE', 'SE') : 1,
                     #Type2 - Nitrogen (amide)
                    ('ASN', 'ND2'): 2 ,('GLN', 'NE2') : 2, ('N'):2,  #backbone N inc N-terminal
                    #Type3 - Nitrogen (aromatic)
                    ('HIS', 'ND1') : 3, ('HIS', 'NE1') : 3, ('TRP', 'NE1') : 3,
                    #Type4 - Nitrogen (guanidinium)
                    ('ARG', 'NE'): 4,('ARG', 'NH1'): 4 ,('ARG', 'NH2'): 4 ,('ARG', 'NH3'): 4 ,
                    #Type5 - Nitrogen (ammonium)
                    ('LYS', 'NZ'): 5 ,
                    #Type6 - Oxygen (carbonyl)
                    ('ASN', 'OD1') : 6 , ('GLN', 'OE1'): 6 , ('O'): 6, # backbone O (should be exept C terminal BUT HOW?)
                    #Type7 - oxygen (hydroxyl)
                    ('SER', 'OG'): 7 , ('THR', 'OG1'): 7, ('TYR', 'OH'): 7 ,
                    #Type8 - Oxygen (carboxyl)
                    ('ASP', 'OD1') : 8,('ASP', 'OD2') : 8, ('ASP', 'OD3') : 8, ('GLU', 'OE1') : 8, 
                    ('GLU', 'OE2') : 8, ('GLU', 'OE3') : 8, ('OXT'):8, #C-terminal O NEEDS TO BE ADDED
                    #Type9 - Carbon (sp2)
                    ('ARG', 'CZ') : 9, ('ASN', 'CG') : 9,('ASP', 'CG') : 9, ('GLN', 'CD') : 9, 
                    ('GLU', 'CD') : 9, ('C'):9, #backbone C
                    #Type10 - carbon (aromatic)
                    ('HIS', 'CG') : 10,('HIS', 'CD2') : 10,('HIS', 'CE1') : 10,('PHE', 'CG') : 10,
                    ('PHE', 'CD1') : 10,('PHE', 'CD2') : 10,('PHE', 'CD3') : 10,('PHE', 'CE1') : 10,
                    ('PHE', 'CE2') : 10,('PHE', 'CE3') : 10,('PHE', 'CZ') : 10,('TRP', 'CG') : 10,
                    ('TRP', 'CD1') : 10,('TRP', 'CD2') : 10, ('TRP', 'CD3') : 10, ('TRP', 'CE1') : 10,
                    ('TRP', 'CE2') : 10,('TRP', 'CE3') : 10,('TRP', 'CZ1') : 10, ('TRP', 'CZ2') : 10,
                    ('TRP', 'CZ3') : 10, ('TRP', 'CH2') : 10,('TYR', 'CG') : 10,('TYR', 'CD1') : 10,
                    ('TYR', 'CD2') : 10,('TYR', 'CD3') : 10,('TYR', 'CE1') : 10,('TYR', 'CE2') : 10,
                    ('TYR', 'CE3') : 10, ('TYR', 'CZ') : 10,
                    # type11
                    ('ALA', 'CB'): 11, ('ARG', 'CB') : 11,('ARG', 'CG') : 11, ('ARG', 'CD') : 11,
                    ('ASN', 'CB'): 11 , ( 'ASP', 'CB') : 11,('CYS', 'CB') : 11, ('GLN', 'CB') : 11, 
                    ('GLN', 'CG'): 11, ('GLU', 'CB') : 11,('GLU', 'CG') : 11 , ('HIS', 'CB') : 11, 
                    ('ILE', 'CB'): 11, ('ILE' , 'CG1') : 11, ('ILE' , 'CG2') : 11, ('ILE' , 'CG3') : 11,
                    ('ILE', 'CD1'): 11, ('LEU' , 'CB') : 11,('LEU' , 'CG') : 11,('LEU' , 'CD1') : 11,
                    ('LEU', 'CD2'): 11,('LEU' , 'CD3') : 11,('LYS' , 'CB') : 11, ('LYS' , 'CG') : 11,
                    ('LYS', 'CD') : 11, ('LYS' , 'CE') : 11,('MET' , 'CB') : 11,('MET' , 'CG') : 11,
                    ('MET', 'CE') : 11,('MSE' , 'CB') : 11,('MSE' , 'CG') : 11,('MSE' , 'CE') : 11,
                    ('PHE', 'CB') : 11, ('PRO' , 'CB') : 11, ('PRO' , 'CG') : 11, ('PRO' , 'CD') : 11,
                    ('SER', 'CB') : 11,('THR' , 'CB') : 11, ('THR' , 'CG2') : 11,('TRP' , 'CB') : 11,
                    ('TYR', 'CB') : 11, ('VAL' , 'CB') : 11,('VAL' , 'CG1') : 11, ('VAL' , 'CG2') : 11, 
                    ('VAL', 'CG3') : 11,('CA'): 11 #backbone CA 
                    }
    
    #Iterate over every atom row and delete if not of the given type. 
    #If the (res, atom_name) tuple doesn't exist as key in atom_types_dict. 
    #See if just the (atom_name) does. If not, print what atom and remove it. 
    #If the key exists, see if the atom type is right, otherwise, remove the row. 

    for index,row in atoms.iterrows():        
        if (row['residue_name'], row['atom_name']) not in atom_types_dict.keys(): 
            if (row['atom_name']) not in atom_types_dict.keys(): 
                #print ('No atom type defined for' , (row['residue_name'], row['atom_name']))
                atoms=atoms[atoms.index != index]
            
            elif atom_types_dict[(row['atom_name'])] != atom_type:
                atoms=atoms[atoms.index != index]
   
        elif atom_types_dict[(row['residue_name'], row['atom_name'])] != atom_type:
            atoms=atoms[atoms.index != index]

    #returns a dataframe with only atoms of the wanted type
    return atoms


# In[26]:


def create_density_map(atoms):
    '''Create a 120x120x120 matrix and compute the atom density for every coordinate'''

     #120x120x120 matrix filled w/ zeroes    
    dens_map=np.zeros((120,120,120))

    # Adding 60 to each coordinate - so they can be used for indexing the matrix
    #Actually subtracting -60, because I don't know how to add. 
    x_s=atoms.x_coord.subtract(-60) 
    y_s= atoms.y_coord.subtract(-60)
    z_s=atoms.z_coord.subtract(-60)

    #Controlling if the atom fits into the 120x120x120 grid
    #If the maximum value of any coordinate is more than 119, the flag is assigned a 1
    #and the function is interrupted
    #Maybe 118 is better, since since the program is based on looking at the two closest neighbours? 
    #What happens at the edges? 

    flag=0

    if x_s.max()>118 or y_s.max()>118 or z_s.max()>118: 
        flag=1
        #return dens_map, flag

    
    #Only goes this far if the protein fits in the grid. 
    #The coordinates (a list of several [x y z]'s) are transposed to a list of l=[x....][y...][z...]   
    coords=[x_s, y_s, z_s]
    coords=list(zip(*coords)) 
    
    #Looping over every coordinate in the grid
    #Since the wanted resolution of the grid is 1Angstrom, the coordinates are converted to integers.
    for i in range(len(coords)):
        x=int(coords[i][0])         #used round before turning into integer before.. Don't think it's nessecary
        y=int(coords[i][1])
        z=int(coords[i][2])
        
        if x>0 and x<118 and y>0 and y<118 and z>0 and z<118:
        

            #For all positions neighbouring (within a distance of 2 Angstrom)  
            for  k in range(x-2,x+3): #need to take +3 to make it get to +2
                for l in range(y-2,y+3):
                    for m in range(z-2,z+3):

            #If the position k,l,m  is the coordinate we are investigating, r=0. 
                        if k==x and l==y and m==z:
                            dens_map[k][l][m]+=exp(0)

        #If the position k,l,m is within 1Angstrom from x,y,z r=1
                        elif ((x-1)<= k <=(x+1)) and ((y-1)<=l<=(y+1)) and ((z-1)<=m<=(z+1)):       
                            dens_map[k][l][m]+=exp(1/2)

            #Otherwise, r=2 (since we're just looking at positions with a maximum distance of 2A)
                        else:
                            dens_map[k][l][m]+=exp(-2)

    #Return density map, and a flag that is 0 if the protein fits in the grid, else 1  
    return dens_map, flag


# In[27]:


def plot_map(dens_map):

    #convert df to list of lists to be able to plot it. 
    #Think there is a np.tolist function that might be better?

    all_data=[]
    for i in range(120):
        for j in range(120):
            for k in range(120):
                if dens_map[i][j][k]!=0:
                    all_data.append([i , j , k , dens_map[i][j][k]])

    #Transpose the list. 
    all_data=np.transpose(all_data)
    x=all_data[0]
    y=all_data[1]
    z=all_data[2]
    density=all_data[3]

    #Pick a colormap 
    cmap=cm.jet
    
    #no idea what's going on here
    den_vals=np.array(density)
    colors=cmap(density)

    density=np.array(density)
    colors[:,-1]=den_vals/den_vals.max()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.scatter(x, y, z, c=colors)
    plt.show()


# In[ ]:


def make_4Darray(density_dict):
    '''Extract density maps from dictionary to np.array'''
    #get the dictionary with atom densities - why not do an array immediately?
    #density_dict=Get_Densities.main(filename)
    
    #Concatenate the 11 maps into one single array
    for key in sorted(density_dict): 
        #The first time, need to  create the new array
        if key==1: 
            density_array=np.expand_dims((density_dict[key]),axis=0)
        else:
            density_array=np.concatenate((density_array, np.expand_dims((density_dict[key]),axis=0)), axis=0)
            
    return density_array


# In[28]:


def main(argv):

    
    ppdb=PandasPdb()
    ppdb=ppdb.read_pdb(argv) 
    atoms=ppdb.df['ATOM']
    
    #relocate centre of gravity to origo
    atoms=centre_model(atoms)
    
    #there are 11 atom types
    atom_types=list(range(1,12)) 
    
    #flag for indication if the protein is to big (having atoms outside of the 120*120*120 grid)
    hint=0
    
    all_density_maps={}
    
    for atom_type in atom_types:
        #create 120x120x120 densitymap (array only containting the atoms for choosen type)
        new_atoms=atoms_to_map(atoms, atom_type) 
        density_map,flag=create_density_map(new_atoms)
        
        #if any of the density maps are to big, the hint will be set to one
        if flag==1: 
            hint=1
        
        #if there are atom densities in the map - add to dictionary and plot 
        if not np.any(density_map):  
            all_density_maps[atom_type]=np.zeros((120,120,120))    
            
        else:
            all_density_maps[atom_type]=density_map
            #plot_map(density_map)
      
    #print if the protein was to large
    if hint==1:
        print (argv, 'had atoms outside the grid')
             
    #make array out of dictionary.. Why making a dict at all? 
    dens_array=make_4Darray(all_density_maps)
    
    return dens_array

if __name__ == '__main__':
    main(sys.argv) 
     

