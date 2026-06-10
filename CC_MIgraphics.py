import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

Mut_inf_A = np.load("MI_Be.npy")/2
Mut_inf_B = np.load("MI_Ne2-2.npy")/2
Mut_inf_C = np.load("MI_Ne2-4.npy")/2

fig, axes = plt.subplots(1, 3, figsize=(15,7))

cmap =  mcolors.LinearSegmentedColormap.from_list("white_blue", ["white", "mediumblue"])
plt.subplots_adjust(bottom=0.18)


r=sns.heatmap(Mut_inf_A, ax=axes[0], cmap=cmap, vmin=0, vmax=0.2, square=True, cbar=False)
r.axvline(6, color='black', linewidth=2)
r.axhline(6, color='black', linewidth=2)
r.axvline(4, ymin=0, ymax=0.5, color='black', linestyle='dashed')
r.axvline(10, ymin=0.5, ymax=1,color='black', linestyle='dashed')
r.axhline(4, xmin=0, xmax=0.5,color='black', linestyle='dashed')
r.axhline(10, xmin=0.5, xmax=1,color='black', linestyle='dashed')
r.invert_yaxis()   

r2=sns.heatmap(Mut_inf_B, ax=axes[1], cmap=cmap, vmin=0, vmax=0.2, square=True, cbar=False)
r2.axvline(12, color='black', linewidth=2)
r2.axhline(12, color='black', linewidth=2)
r2.axvline(6, ymin=0, ymax=0.5, color='black', linestyle='dashed')
r2.axvline(8, ymin=0, ymax=0.5, color='black', linestyle='dashed')
r2.axvline(12, ymin=0.5, ymax=1,color='black', linestyle='dashed')
r2.axvline(18, ymin=0.5, ymax=1,color='black', linestyle='dashed')
r2.axvline(20, ymin=0.5, ymax=1,color='black', linestyle='dashed')
r2.axhline(6, xmin=0, xmax=0.5,color='black', linestyle='dashed')
r2.axhline(8, xmin=0, xmax=0.5,color='black', linestyle='dashed')
r2.axhline(12, xmin=0.5, xmax=1,color='black', linestyle='dashed')
r2.axhline(18, xmin=0.5, xmax=1,color='black', linestyle='dashed')
r2.axhline(20, xmin=0.5, xmax=1,color='black', linestyle='dashed')
r2.invert_yaxis() 

r3=sns.heatmap(Mut_inf_C, ax=axes[2], cmap=cmap, vmin=0, vmax=0.2, square=True, cbar=False)
r3.axvline(12, color='black', linewidth=2)
r3.axhline(12, color='black', linewidth=2)
r3.axvline(6, ymin=0, ymax=0.5, color='black', linestyle='dashed')
r3.axvline(8, ymin=0, ymax=0.5, color='black', linestyle='dashed')
r3.axvline(12, ymin=0.5, ymax=1,color='black', linestyle='dashed')
r3.axvline(18, ymin=0.5, ymax=1,color='black', linestyle='dashed')
r3.axvline(20, ymin=0.5, ymax=1,color='black', linestyle='dashed')
r3.axhline(6, xmin=0, xmax=0.5,color='black', linestyle='dashed')
r3.axhline(8, xmin=0, xmax=0.5,color='black', linestyle='dashed')
r3.axhline(12, xmin=0.5, xmax=1,color='black', linestyle='dashed')
r3.axhline(18, xmin=0.5, xmax=1,color='black', linestyle='dashed')
r3.axhline(20, xmin=0.5, xmax=1,color='black', linestyle='dashed')
r3.invert_yaxis() 

axes[0].tick_params(axis='y', labelrotation=0)
axes[1].tick_params(axis='y', labelrotation=0)
axes[2].tick_params(axis='y', labelrotation=0)

axes[1].set_yticks(np.arange(0, 24, 2))
axes[1].set_yticklabels(np.arange(0, 24, 2))
axes[2].set_yticks(np.arange(0, 24, 2))
axes[2].set_yticklabels(np.arange(0, 24, 2))

axes[0].set_title(r'$^{8}\mathrm{Be}$', fontsize=16)
axes[1].set_title(r'$^{20}\mathrm{Ne}$', fontsize=16)
axes[2].set_title(r'$^{22}\mathrm{Ne}$', fontsize=16)

cbar = fig.colorbar(r2.collections[0],   # agafa el mapeig del heatmap
    ax=axes, orientation='horizontal',
    fraction=0.03, pad=0.1)

#cbar.set_label(r'$S_{ij}$', fontsize=10)
cbar.ax.set_ylabel(r'$S_{ij}$', fontsize=14, rotation=0, labelpad=20)
cbar.ax.yaxis.set_label_coords(1.1, -0.4)


#plt.tight_layout()
plt.show()