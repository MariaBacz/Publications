# Importing the necessary libraries for plotting
import matplotlib.pyplot as plt
import numpy as np

# Loading data from the rmsf.xvg file, ignoring lines starting with "@" or "#"
res,rmsf = np.loadtxt("tri_sol_full_nojump_rmsd.xvg", comments=["@", "#"], unpack=True)
res1,rmsf1 = np.loadtxt("tri_longer_full_nojump_rmsd.xvg", comments=["@", "#"], unpack=True)

# Creating a figure and axis objects for the line plot
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)

# Plotting the line representing the RMSD values
ax.plot(res,rmsf, color="black", label='Original simulation') 
ax.plot(res1,rmsf1, color="grey", linestyle="-", label='Replica') 

# Setting labels for the x-axis (residue) and y-axis (RMSD value)
plt.rcParams['font.family'] = 'Ariel'
ax.set_xlabel("Time (µs)", fontsize=18) 
ax.set_ylabel(r"RMSD (nm)", fontsize=18) 
plt.legend(fontsize=15, loc = 'upper right')

# Saving the plot as a PNG image with higher resolution (300 dpi)
plt.title('RMSD of the Lcl NTD Trimer in Solution', fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(0, 1.2)
plt.tight_layout()
plt.savefig("RMSD of the Lcl NTD Trimer in Solution.png", format="png", dpi=300)      
