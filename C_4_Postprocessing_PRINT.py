#_____________________________________________________________________________________________________________
# C_4 POSTPROCESSING, CALCULATING PROBABILITIES
#_____________________________________________________________________________________________________________


from C_3_ACES_Ledger_PRINT import *

#____________________________________________________________________________________________________________
# 1. Symplectic Representation, Symplectic Product

def symplectic(pauli, bool = False):
    """
        Returns tuple of sign and symplectic representation of pauli: symplectic('-XY') = (-1, [[1, 0], [1, 1]])
    INPUT:
        |Variable    | Type           | Comment
        |____________|________________|____________________________________________________
        | pauli      | str            | String of Paulis, e. g. 'XIY'
        | bool       | bool           | It True returns [True, False] instead of [0, 1]
        |            |                |    Default: False
    """
    P = []
    SymplecticP = []
    if type(pauli) == str:
        if pauli[0]=='-':
            sgn = -1
            pauli = pauli[1::]
        else:
            sgn = 1
        for i in pauli:
            P.append(qinf.Pauli(i))
    for i in P:
        SymplecticP.append([qinf.Pauli(i).x[0], qinf.Pauli(i).z[0]])
    if bool == True:
        return sgn, SymplecticP
    else:
        nonboolP = []
        for i in SymplecticP:
            subpair = []
            for j in i:
                if j == True:
                    subpair.append(1)
                else:
                    subpair.append(0)
            nonboolP.append(subpair)
        return sgn, nonboolP

def symplectic_product(pauli1, pauli2, sgn = False):
    """
    RETURNS:
        Symplectic product of two Pauli words: <a, b> = a^T (𝜰 + 𝜰^T) b mod2
        P_a given by its symplectic representation a = a_1 a_2 ... a_2n
        such that P_a ∝ X^a_1 Z^a_2 ... X^a_2n-1 Z^a_2n
        with        n
        (𝜰 + 𝜰^T) = ⊕  [[0, 1], [1, 0]]
                   k=1

    INPUT:
        |Variable    | Type           | Comment
        |____________|________________|____________________________________________________
        | pauli1,    | str            | String of Paulis, e. g. 'XIY'
        |   pauli2   | str            |
        | sgn        | bool           | Carry sign from paulis. Default: False
    """
    a = symplectic(pauli1, bool=False)
    b = symplectic(pauli2, bool=False)
    Ymatrix = array([[0, 1], [1, 0]])
    if (sgn == True) and (a[0]*b[0]==-1):
        res = 1
    else:
        res = 0
    for i in eachindex(a[1]):
        res += array(a[1])[i].dot(Ymatrix.dot(array(b[1])[i]))
    return res % 2

#____________________________________________________________________________________________________________
# 2. Probabilities Walsh-hadamard Transform

# ___2.1 Calculate probabilities

def probabilities(lambdas, paulierrors = ['I', 'X', 'Y', 'Z']):
    """
    RETURNS:
        Probabilities from circuit eigenvalues for a single-gate experiment, e.g.
        probabilities({'I': 1.0, 'X': 0.4299, 'Y': 0.7699, 'Z': 0.58}, ['I', 'X', 'Y', 'Z']) =
        {'I': 0.695, 'X': 0.0199, 'Y': 0.1899, 'Z': 0.095}

    INPUT:
        |Variable     | Type | Comment
        |_____________|______|______________________________________________________________
        | lambdas     | dict | Circuit eigenvalues dict from ACES experiments
        | paulierrors | list | List of Pauli errors we want to learn the probabilities for
    """
    P = {}
    for i in paulierrors:
        p = 0
        for j in lambdas:
            p += (-1)**symplectic_product(i, j)*lambdas[j]/len(paulierrors)
        P[i] = p
    return P

# ___2.2 Visualize Probabilities
import matplotlib.pyplot as plt

def visualize_probs(probs):
    positions = range(0, len(list(probs.values())))
    labels = list(probs.keys())
    plt.xticks(positions, labels)
    plt.plot(range(0, len(list(ps.values()))), list(ps.values()), '+')
    plt.ylim([0,1])
