import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals, distances, contacts
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable

# Load input files (tpr and xtc)
u = mda.Universe("topology.tpr", "trajectory.xtc")

# Select the groups of atoms
NTD_A_sel = u.select_atoms("segid seg_0_PROA and name CA")
print("protein helix 1:", len(NTD_A_sel))
  
NTD_B_sel = u.select_atoms("segid seg_1_PROB and name CA")
print("protein helix 2:", len(NTD_B_sel))

NTD_C_sel = u.select_atoms("segid seg_2_PROC and name CA")
print("protein helix 2:", len(NTD_C_sel))

# Contacts between NTD_A_sel and NTD_B_sel
def contacts_within_cutoff(u, NTD_A_sel, NTD_B_sel, radius=10):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(NTD_A_sel.positions, NTD_B_sel.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts]) 
    return np.array(timeseries)

NTD_A_B = contacts_within_cutoff(u, NTD_A_sel, NTD_B_sel, radius=10)
NTD_A_B_df = pd.DataFrame(NTD_A_B, columns=['Frame',
                                  'Contacts'])
print(NTD_A_B_df.head())

# Contacts between NTD_C_sel and NTD_B_sel
def contacts_within_cutoff(u, NTD_C_sel, NTD_B_sel, radius=10):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(NTD_C_sel.positions, NTD_B_sel.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

NTD_C_B = contacts_within_cutoff(u, NTD_C_sel, NTD_B_sel, radius=10)
NTD_C_B_df = pd.DataFrame(NTD_C_B, columns=['Frame',
                                  'Contacts'])
print(NTD_C_B_df.head())

# Contacts between NTD_C_sel and NTD_A_sel
def contacts_within_cutoff(u, NTD_C_sel, NTD_A_sel, radius=10):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(NTD_C_sel.positions, NTD_A_sel.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

NTD_C_A = contacts_within_cutoff(u, NTD_C_sel, NTD_A_sel, radius=10)
NTD_C_A_df = pd.DataFrame(NTD_C_A, columns=['Frame',
                                  'Contacts'])
print(NTD_C_A_df.head())

# Ploting
NTD_C_A_df.plot(x='Frame', color="black", figsize=(7, 3), xlim=(0, 10000), ylim=(0, 100), ylabel=('Number of Contacts'), fontsize= 13)
plt.tick_params(axis='both', labelsize=13)
plt.xlabel('Frame', fontsize=13)
plt.ylabel('Number of Contacts', fontsize=13)
plt.title('Contacts between monomer 1 and 2 within 10Å cutoff', fontsize=15)
plt.savefig("Contacts between monomer 1 and 2 within 10Å cutoff.png", dpi=1200)
plt.rcParams["font.family"] = "Arial"

NTD_C_B_df.plot(x='Frame', color="black", figsize=(7, 3), xlim=(0, 10000), ylim=(0, 100), ylabel=('Number of Contacts'), fontsize= 13)
plt.tick_params(axis='both', labelsize=13)
plt.xlabel('Frame', fontsize=13)
plt.ylabel('Number of Contacts', fontsize=13)
plt.title('Contacts between monomer 1 and 3 within 10Å cutoff', fontsize=15)
plt.savefig("Contacts between monomer 1 and 3 within 10Å cutoff.png", dpi=1200)
plt.rcParams["font.family"] = "Arial"

NTD_A_B_df.plot(x='Frame', color="black", figsize=(7, 3), xlim=(0, 10000), ylim=(0, 120), ylabel=('Number of Contacts'), fontsize= 13)
plt.tick_params(axis='both', labelsize=13)
plt.xlabel('Frame', fontsize=13)
plt.ylabel('Number of Contacts', fontsize=13)
plt.title('Contacts between monomer 2 and 3 within 10Å cutoff', fontsize=15)
plt.savefig("Contacts between monomer 2 and 3 within 10Å cutoff.png", dpi=1200)
plt.rcParams["font.family"] = "Arial"