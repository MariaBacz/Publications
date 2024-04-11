import MDAnalysis as mda
import MDAnalysis.analysis.distances
import MDAnalysis.analysis.rdf as rdf_calc

import numpy as np
from numpy import linalg as LA
from numpy.linalg import norm

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm.auto import tqdm

# File and folder names
u = mda.Universe("topology.tpr", "trajectory.xtc")
u1 = u.select_atoms("protein")

#Funcion for contacts calculations (code by Raquel Lopez-Rios De Castro and Chris Lorenz, modified by Maria Baczynska)
def contact_matrix(protein_1, protein_2):
    
    contacts = np.zeros((len(np.unique(protein_1.resids)), len(np.unique(protein_2.resids))))
    
    
    if np.unique(protein_1.resids)[0] != 0: #in case the first resid of the protein is not 0
        
        sub_1 = np.unique(protein_1.resids)[0]
        
    if np.unique(protein_2.resids)[0] != 0:
        
        sub_2 = np.unique(protein_2.resids)[0]
        
    else:
        
        sub_1 = 1
        sub_2 = 1
    
    

    for ts in tqdm(u.trajectory[::50]):

        for i in (np.unique(protein_1.resids)):
    
            a = protein_1.select_atoms('resid '+str(i))
        
            for k in (np.unique(protein_2.resids)): 
        
                b = protein_2.select_atoms('resid '+str(k))
        
                if MDAnalysis.analysis.distances.distance_array(a.center_of_mass(), b.center_of_mass(), box=[u.dimensions[0], u.dimensions[1], u.dimensions[2], u.dimensions[3], u.dimensions[4], u.dimensions[5]]) < 10.0:
            
                    contacts[i-sub_1,k-sub_2] += 1
    
    return contacts

#Use the contacs funcion on each selection id est helix within the trimer Lcl NTD
NTD_A_sel = u.select_atoms("segid seg_0_PROA")
print("protein helix 1:", NTD_A_sel)
  
NTD_B_sel = u.select_atoms("segid seg_1_PROB")
print("protein helix 2:", NTD_B_sel)

NTD_C_sel = u.select_atoms("segid seg_2_PROC")
print("protein helix 2:", NTD_C_sel)

contacts_CB_1 = contact_matrix(NTD_C_sel, NTD_B_sel)
contacts_AB_1 = contact_matrix(NTD_A_sel, NTD_B_sel)
contacts_AC_1 = contact_matrix(NTD_A_sel, NTD_C_sel)

contacts = (np.array(contacts_CB_1) + np.array(contacts_AB_1) + np.array(contacts_AC_1)) / 3

print("contacs:", contacts, contacts.shape)

# Plotting heatmap
string_array_columns=[]
for i in range(len(u1[0:30])):
    string_array_columns.append(str(u1.residues[i]))
    
string_array_rows=[]
for i in range(len(u1[0:30])):
    string_array_rows.append(str(u1.residues[i]))

#Create a pandas dataframe
df = pd.DataFrame(contacts, string_array_columns)

xticks= string_array_columns
xticks_pre= [x[:-1] for x in string_array_columns]
xticks= [x[+8:] for x in xticks_pre]

yticks= string_array_rows
yticks_pre= [x[:-1] for x in string_array_rows]
yticks= [x[+8:] for x in yticks_pre]

print(xticks)

# Plotting
fig, ax = plt.subplots(figsize=(11,8))
sns.heatmap(df, xticklabels=xticks, yticklabels=yticks, cmap="Greys", ax=ax, linewidths=.5, vmin=0, vmax=100)

ax.invert_yaxis() 

plt.rcParams["font.family"] = "Arial"
ax.tick_params(labelsize=18)
cax = ax.figure.axes[-1]
cax.tick_params(labelsize=18)

plt.xlabel("Residues", fontsize = 22)
plt.ylabel("Residues", fontsize = 22)
plt.title("Number of contacs between the monomers", fontsize = 27, x=0.5, y=1.02)

# Save the plot
plt.tight_layout()
plt.savefig("Contacs_within_trimer_replicated.png")
