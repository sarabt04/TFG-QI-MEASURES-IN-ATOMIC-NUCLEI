import numpy as np
from scipy.sparse import coo_matrix
import math
from scipy.sparse.linalg import eigsh
from itertools import combinations
import random

#0.1-CONSTUÏR+DIAGONALITZAR LA MATRIU------------------------------------------------------
def ground_state(data_matrix):
    files=data_matrix[:, 0].astype(int)
    columnes=data_matrix[:, 1].astype(int)
    valors=data_matrix[:, 2]
    n = max(files.max(), columnes.max()) + 1
    H_sparse = coo_matrix((valors, (files, columnes)), shape=(n, n))

    E, vecs = eigsh(H_sparse, k=5, which='SA')
    vec1=vecs[:,0]   #Estat fonamental
    return vec1

#0.2-CONVERTIR AUTOESTATS A BITSTRINGS-----------------------------------------------------  
def occ(many_body_config):
    occ=[]
    for line in many_body_config:
        parts=line.split(",",1)
        orbital_str=parts[1].strip().strip("()")
        orbital_tupple=[int(x) for x in orbital_str.split(",")]
        occ.append(orbital_tupple)  
    return occ
 
def bitstrings(occ, d):
    bitstrings=[]
    for i in occ:
        buit=[0]*d
        for j in i:
            buit[j]=1
        bitstrings.append(buit)   
    return bitstrings

#1-ENTROPIA DE VON NEUMANN PER SINGLE-PARTICLE ORBITAL-------------------------------------
def sing_orb(vec1, dim_vs, occ):
    prob = np.abs(vec1)**2
    gamma_i=np.zeros(dim_vs)
    for idx,element in enumerate(occ):
        gamma_i[element]+=prob[idx]
    S_i=[-gamma*math.log2(gamma)-(1-gamma)*math.log2(1-gamma) for gamma in gamma_i]
    return gamma_i, S_i

#2.MUTUAL INFORMATION----------------------------------------------------------------------
def bits_to_int(bits):
    idx = 0
    for b in bits:
        idx = 2*idx + b
    return idx

def Von_Neumann_entr(sing_values):
    S=0
    for i in sing_values:
        eigenv=i**2
        if i>1e-5:
            S+=-eigenv*math.log2(eigenv)
    return S

def Mutt_inf(d, vec1, S_sing_orb, bitstrings):
    Mut_inf=np.zeros((d,d))
    for i in range(len(bitstrings[0])):
        for j in range(i):
            #Singular Value Decomposition
            idx_A=[]
            idx_B=[]
            element=[]
            for idx, a in enumerate(vec1):
                total=bitstrings[idx]
                partA=[x for k,x in enumerate(total) if k in (i,j)]
                partB=[x for k,x in enumerate(total) if k not in (i,j)]
                idx_A.append(bits_to_int(partA))
                idx_B.append(bits_to_int(partB))
                element.append(a)

            #Reenumerate B indices
            _, compact_idxB = np.unique(idx_B, return_inverse=True)
            C=np.zeros((max(idx_A)+1, compact_idxB.max()+1))
            np.add.at(C, (idx_A, compact_idxB), element)
                
            rhoA = C @ C.conj().T
            eigvals = np.linalg.eigvalsh(rhoA)
            eigvals = eigvals[eigvals > 1e-12]
            S_joint=-np.sum(eigvals*np.log2(eigvals))
            
            S_mut=S_sing_orb[i]+S_sing_orb[j]-S_joint
            
            Mut_inf[i,j]=S_mut
            Mut_inf[j,i]=S_mut
    return Mut_inf

#3.EQUIPARITIONS----------------------------------------------------------------------
def single_equip_entropy(Vec, p_list, bitstrings):
    idxA=[]
    idxB=[]
    element=[]
    for idx, a in enumerate(Vec):
        total=bitstrings[idx]
        partA=[x for k,x in enumerate(total) if k in p_list]
        partB=[x for k,x in enumerate(total) if k not in p_list]
        idxA.append(bits_to_int(partA))
        idxB.append(bits_to_int(partB))
        element.append(a)

    #Reenumerate A, B indices
    _, compact_idxA = np.unique(idxA, return_inverse=True)
    _, compact_idxB = np.unique(idxB, return_inverse=True)
    C_equip=np.zeros((compact_idxA.max()+1, compact_idxB.max()+1))
    np.add.at(C_equip, (compact_idxA, compact_idxB), element)
        
    singval_equip=np.linalg.svd(C_equip, compute_uv=False)
    return Von_Neumann_entr(singval_equip)

def all_equip(d, percentage):
    TOTAL=list(range(d))
    comb=list(combinations(TOTAL, d//2))
    Combinations=[x for k,x in enumerate(comb) if 0 in x]
    
    n_sample=int(len(Combinations)*percentage)
    mostra_random=random.sample(Combinations, n_sample)
    return mostra_random

def all_equip_entropy(Vec, bitstrings, d, percentage):
    S_equip=[]
    Partitions=all_equip(d,percentage)
    for part in Partitions:
        part_s=set(part)
        S_equip_i=single_equip_entropy(Vec, part_s, bitstrings)
        S_equip.append(S_equip_i)
    return S_equip