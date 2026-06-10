#ADVERTISMENT: this code takes a long time to compile, because the input arguments to the function all_equipartition_entropy are huge
import numpy as np

#1-DADES INICIALS: MATRIU, ESTATS, DIMENSIÓ VALENCE SPACE-----------------------------------
data=np.loadtxt(r"data\hamiltonian_matrix_sd_2_4.dat")
with open(r"data\many_body_states_sd_2_4.dat") as f:
    many_body = f.read().splitlines()[1:]  # salta la primera línia
d=24

#2-ENTANGLEMENT MEASURES--------------------------------------------------------------------
from AA_Functions_ent import ground_state, occ, bitstrings, sing_orb, Mutt_inf
Vec=ground_state(data)
Occ=occ(many_body)
Bitstrings=bitstrings(Occ, d)

#2.1-Single-particle orbital entropy
Occ_prob=sing_orb(Vec, d, Occ)[0]
S_sing_orb=sing_orb(Vec, d, Occ)[1]

print("Occupation probability=", Occ_prob)
print("Single particle entropy=", S_sing_orb)

#2.2-Mutual information
MI=Mutt_inf(d, Vec, S_sing_orb, Bitstrings)
np.save("MI_Ne2-4.npy", MI)

#2.3-Equipartitions p-n, m<0-m>0
from AA_Functions_ent import single_equip_entropy, all_equip, all_equip_entropy

protons=[0,1,2,3,4,5,6,7,8,9,10,11]
Sp=single_equip_entropy(Vec, protons, Bitstrings)
print(r"S_{pn}=", Sp)

mpos=[0,1,2,6,8,9,12,13,14,18,20,21]
Sm=single_equip_entropy(Vec, mpos, Bitstrings)
print(r"S_{m}=", Sm)

#2.4-All Equipartitions
perc=0.01
T=all_equip(d, perc)
S_all_eq=all_equip_entropy(Vec, Bitstrings, d, perc)
np.save("Eq_Ne2-4.npy", S_all_eq)
