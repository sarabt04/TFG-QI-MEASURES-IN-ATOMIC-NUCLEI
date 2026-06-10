import numpy as np

#1-DADES INICIALS: MATRIU, ESTATS, DIMENSIÓ VALENCE SPACE-----------------------------------
data=np.loadtxt(r"data\hamiltonian_matrix_p_2_2.dat")
with open(r"data\many_body_states_p_2_2.dat") as f:
    many_body = f.read().splitlines()[1:]  # salta la primera línia
d=12

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
np.save("MI_Be.npy", MI)

#2.3-Equipartitions p-n, m<0-m>0
from AA_Functions_ent import single_equip_entropy, all_equip, all_equip_entropy

protons=[0,1,2,3,4,5]
Sp=single_equip_entropy(Vec, protons, Bitstrings)
print(r"S_{pn}=", Sp)

mpos=[0,1,4,6,7,10]
Sm=single_equip_entropy(Vec, mpos, Bitstrings)
print(r"S_{m}=", Sm)

#2.4-All Equipartitions
perc=1
T=all_equip(d, perc)
S_all_eq=all_equip_entropy(Vec, Bitstrings, d, perc)
np.save("Eq_Be.npy", S_all_eq)

#3-MAGIC MEASURES---------------------------------------------------------------------------
from AA_Functions_magic import Cp_list, Magic
#lo unic q canviaria es obtenir directament els non-zero ps o fer un afuncio per aixo
cp=Cp_list(Vec, Bitstrings)
Cp_nz=[]
for key, value in cp.items():
    if value!=0:
        Cp_nz.append(value)
        
Mlin=Magic(Cp_nz, d)[0]
M1=Magic(Cp_nz, d)[1]
M2=Magic(Cp_nz, d)[2]

np.save("PS_Be8.npy", Cp_nz)
print(r"M$_{lin}$=", Mlin)
print(r"M$_{1}$=", M1)
print(r"M$_{2}$=", M2)

import math
th_num_PS=4**d
Card=th_num_PS-len(Cp_nz)
print(Card)
nu=math.log2(Card/(2**d))
print(nu)