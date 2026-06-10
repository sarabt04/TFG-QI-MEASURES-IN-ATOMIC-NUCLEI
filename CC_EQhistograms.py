import numpy as np
import matplotlib.pyplot as plt

EQ_A = np.load("EQ_Be.npy")
EQ_B = np.load("Eq_Ne2-2.npy")
EQ_C = np.load("Eq_Ne2-4.npy")

fig, axes = plt.subplots(1, 3, figsize=(16,2.5))

datasets = [EQ_A, EQ_B, EQ_C]
titols = [r'$^8$Be', r'$^{20}$Ne', r'$^{22}$Ne']
bins = np.arange(1, 7 + 0.2, 0.2)
x0_vals = [1.98, 3.10, 2.57]
x1_vals = [3.67, 5.08, 5.8]

for ax, data, titol, x0, x1 in zip(axes, datasets, titols, x0_vals, x1_vals):
    new_values=[v for v in [x0, x1] if v not in data]
    data_final=np.append(data, new_values)
    n, bins_edges, patches = ax.hist(data_final, bins=bins, color="orange", ec="black")
    
    idx0 = np.digitize(x0, bins_edges) - 1
    if 0 <= idx0 < len(patches):
        patches[idx0].set_facecolor('brown')
        
    idx1 = np.digitize(x1, bins_edges) - 1         
    if 0 <= idx1 < len(patches):
        patches[idx1].set_facecolor('olive')
        
    ax.set_title(titol, fontsize=20)
    ax.set_yscale("log")
    ax.set_xlabel("S")
    
axes[0].set_xlim(1,4.5)
axes[0].set_ylabel(r"$N_p$")
axes[0].set_ylim(0,5e2)
axes[1].set_xlim(1,7)
axes[2].set_xlim(1,7)

plt.tight_layout()
plt.show()