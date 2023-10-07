from C_0_1_Packages_PRINT import *

# ___0.1 Iterate over whole List

def eachindex(Array):
    return list(range(len(Array)))

# ___0.2 Count Occurances in Lists

def occurance_counter(lst):
    """
    RETURNS:
        Dict with unique values and numbers of occurances in list, e.g.
        occurance_counter(lst) = {'>0': 4, '<0': 3}

    INPUT:
        |Variable   | Type          | Comment
        |___________|_______________|______________________________________________________
        | lst       | list          | Any list
    """
    keys = list(Counter(lst).keys()); values = list(Counter(lst).values())
    return {keys[i]: values[i] for i in eachindex(keys)}


# ___0.3 Shortscut Matrix Multiplication

def mm(Cl1, Cl2): # Matrix multilication
    return matmul(qinf.Operator(Cl1), qinf.Operator(Cl2))

# ___0.4 List flattening

def flatlist(x):
    """
    RETURNS:
        Flattens list of lists to by depth 1.
        flatlist([['A', 'B', 'C'], ['B', 'C', 'A'], ['C', 'A', 'B']]) =
        ['A', 'B', 'C', 'B', 'C', 'A', 'C', 'A', 'B']
    """
    res = x[0]
    for i in x[1::]:
        for j in i:
            res.append(j)
    return res

# ___0.5 Transforming outputs

def tuptostr(list):
    """
    INPUT:
        Tuple of tuples, e.g. (('I','X'),( 'I',),( 'X',)).

    RETURNS:
        List of strings, e.g. ['IX', 'I', 'X'].
    """
    res = []
    for i in eachindex(list):
        res.append(''.join(list[i]))
    return res
    del res

# ___0.5 Rotating permutations of list

def rotperms(x):
    """
    RETURNS:
        All rotation permutations of elements of list x:
        rotperms(['A', 'B', 'C']) = [['A', 'B', 'C'], ['B', 'C', 'A'], ['C', 'A', 'B']]
    """
    rotpermutations = []
    for i in eachindex(x):
        res = []
        for j in eachindex(x):
            #print((i + j) % len(x))
            res.append(x[(i + j) % len(x)])
        #print(res)
        rotpermutations.append(res)
    return rotpermutations
