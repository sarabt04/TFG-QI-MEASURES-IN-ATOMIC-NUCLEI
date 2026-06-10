from collections import defaultdict
from itertools import product
import math

#5-MAGIC---------------------------------------------------------------------------------
def local_choices(ket, bra):
    if ket==bra: #Only operators that will give a contribution neq0 I or Z
        return [(0, 1), (3, 1 if ket == 0 else -1),]

    else: #Only operators that will give a contribution neq0 X or Y
        return [(1, 1), (2, 1j if ket == 0 else -1j),]


def Cp_list(vec, bitstring):
    Cp_list=defaultdict(float)
    
    for idx1,a in enumerate(vec):
        branomod=bitstring[idx1]
        for idx2,b in enumerate(vec):
            kettomod=bitstring[idx2]
            
            PS=[local_choices(k,b) for k,b in zip(kettomod,branomod)]  
            
            for op in product(*PS):
                phase=1
                pauli_lab=[]
                for sigma, ph in op:
                    phase*=ph
                    pauli_lab.append(sigma)
                    
                value=a*b*phase
                val=value.real
                Cp_list[tuple(pauli_lab)] += val

    return Cp_list

def Magic(cp_values,d):
    D=2**(d)

    Xi_p=[i**(2)/D for i in cp_values]
    Xi_p2=[i**2 for i in Xi_p]    
    Xi_pd=[i*D for i in Xi_p]
    Mlin=1-D*sum(Xi_p2)
    M1=-sum(p*math.log2(q) for p, q in zip(Xi_p, Xi_pd))
    M2=-math.log2(D*sum(Xi_p2))
    
    return Mlin, M1, M2


