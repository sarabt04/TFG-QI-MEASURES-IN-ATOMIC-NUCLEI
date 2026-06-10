import numpy as np
import matplotlib.pyplot as plt

PS = np.load("PS_Be8.npy")

bins = np.linspace(-1, 1, 40)
plt.hist(PS, bins=bins, color="orange", ec="black", lw=2)
plt.xlim(-1,1)
plt.xlabel(r'$\langle\psi|P|\psi\rangle$', size=10)
plt.ylabel("Frequency", size=10)
plt.show()