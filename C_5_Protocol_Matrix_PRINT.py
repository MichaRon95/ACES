#_____________________________________________________________________________________________________________
# C_5 CONSTRUCT PROTOCOL MATRIX
#_____________________________________________________________________________________________________________


from C_4_Postprocessing_PRINT import *

#_____________________________________________________________________________________________________________
# 1. Protocol Matrix Row

def protocol_matrix_row(gates, input_pauli, noId = False, paulis = ['I', 'X', 'Y', 'Z'], sign = False):
    """
    RETURNS:
        Row of protocol matrix for a pair of a circuit give by gates and an input Pauli, e.g.
        protocol_matrix_row(randclifflist((1, 2)), 'X', paulis = ['I', 'X', 'Y']) =
        [1, -1, -1, 1, -1, -1, -1, 1, 1]

    INPUT:
        |Variable | Type           | Comment
        |_________|________________|___________________________________________________________
        | gates   | list(Clifford) | List of Clifford gates defining a cicuit.
        | pauli   | str            | String of Paulis, e. g. 'XIY'
        | noId    | bool           | If = True all Pauliwords with I are excluded
        |         |                |     Default noId = False
        | paulis  | str,           | Possible choice of Pauli-errors as str or list of str.
        |         |  list(str)     |     Default: paulis = ['I', 'X', 'Y', 'Z']
        | sgn     | bool           | Carry sign from paulis.
        |         |                |     Default: False

    """
    dims = (len(gates[0].output_dims()), len(gates)); dims
    perrors = allpauliconfigs(dims, paulis = ['X', 'Y', 'Z']) if noId == True else allpauliconfigs(dims, paulis); perrors
    if type(input_pauli) == str:
        input_pauli = qinf.Pauli(input_pauli); input_pauli
    pauli_circuit = [input_pauli]
    for i in eachindex(gates):
        pauli_circuit.append(pauli_circuit[i].evolve(gates[i]))
    pauli_circuit
    sign_vector = []
    for i in perrors:
        exponent = 0
        for j in eachindex(i.split(', ')):
            exponent += symplectic_product(pauli_circuit[j].to_label(), i.split(', ')[j], sgn = sign)
        sign_vector.append((-1)**exponent)
    return sign_vector

# Legacy functions
def sign_probs(gates, input_pauli, noId = False, paulis = ['I', 'X', 'Y', 'Z'], sign = False):
    return protocol_matrix_row(gates, input_pauli, noId = False, paulis = ['I', 'X', 'Y', 'Z'], sign = False)

#_____________________________________________________________________________________________________________
# 2. Devise ACES Protocol

def protocol(dims, rank, p_errs):
    """
    RETURNS:
        Tuple of (protocol, protocol matrix). The protocol is a list of pairs of (circuit, input Pauli),
        with which a line from the protocol matrix is associated. The circuit is a list of Cliffords and
        the input Pauli is a string. E.g.
        protocol((1, 2), 4, ['I', 'Y']) =
            ([(circuit1, 'X'), (circuit2, 'Y'), (circuit3, 'Z'), (circuit4, 'Y')],
                [[1, -1, -1, 1],
                 [1, -1, 1, -1],
                 [1, 1, -1, -1],
                 [1, 1, 1, 1]])
    INPUT:
        | Variable | Type      | Comment
        |__________|___________|___________________________________________________________
        | dims     | list(int) | Dimensions of dims[0]-qubit dims[1]-gate circuit.
        | rank     | int       | Rank we demand the protocol matrix to have
        | p_errs   | str,      | Possible Pauli-errors as str or list of str.
        |          |  list(str)|
    """
    assert rank < ((4**dims[0]-1)**dims[1]+1), 'ERROR: Rank exceeds limit of (4^n-1)^m+1 for circuit dimensions (n, m).'

    pair_1 = (randclifflist(dims), randpauliword(dims[0], noId = True))              # First pair (P, C)
    protocol = [pair_1]                                                              # Begin protocol with 1st pair
    protocol_matrix = [protocol_matrix_row(pair_1[0], pair_1[1], paulis = p_errs)]   # First protocol mat. row

    rk = linalg.matrix_rank(protocol_matrix)                                         # Rank counter
    while rk < rank:
        circ = randclifflist(dims); pauli = randpauliword(dims[0], noId = True)      # New pair
        M = protocol_matrix.copy()                                                   # Test copy of prot. mat
        M.append(protocol_matrix_row(circ, pauli, paulis = p_errs))                  # Add new row from new pair
        new_rank = linalg.matrix_rank(M)                                             # Calculate new rank
        pauli2 = qinf.Pauli(pauli).evolve(circ[0]).to_label()                        # Exclude Paulis evolving to neg Paulis
        if ((circ, pauli) not in protocol) and (new_rank>rk) and (pauli2[0] != '-'): # If rank increases...
            protocol_matrix = M                                                      # Add row to protocol matrix
            rk = new_rank
            protocol.append((circ, pauli))                                           # Add pair to protocol

    return protocol, protocol_matrix

#_____________________________________________________________________________________________________________
# 3. Run Experimental Protocol

def run_protocol(exp_protocol, runs = 512, simulate_noise = True, weights = None, backend = Aer.get_backend('aer_simulator')):
    """
    RETURNS:
        List of circuit eigenvalues from an experimental protocol.
        Prints relative frequencies of experimental outcomes. E.g.
        ptcl = protocol(dims, 4, ['I', 'Y'])
        run_protocol(ptcl[0], weights = ws) =
            [-0.2655527580682079, -1.0, -0.23591457923096884, -0.4157695130902863]
        console:
            {'0': {'1': 0.642, '0': 0.358}, '1': {'0': 0.624, '1': 0.376}}
            {'u': {'-': 1.0}, 'd': {'+': 1.0}}
            {'0': {'0': 0.421, '1': 0.579}, '1': {'1': 0.343, '0': 0.657}}
            {'u': {'-': 0.731, '+': 0.269}, 'd': {'+': 0.684, '-': 0.316}}
            4

    INPUT:
        | Variable       | Type         | Comment
        |________________|______________|______________________________________________________________
        | exp_protocol   | list(tuples) | list of pairs of (circuit, input Pauli).
        | runs           | int          | Rank we demand the protocol matrix to have
        |                |              |   Default: runs = 512
        | simulate_noise | bool         | If = True, noise simulated with weights
        |                |              |   Default: simulate_noise = True
        | weights        | W class      | Noise distribution defined through W class
        |                |              |   Default: weights = None, uniform error distribution
        | backend        | qiskit       | Run experiments on a simulator or on qiskit quantum hardware
        |                |              |   Default: backend = Aer.get_backend('aer_simulator')
    """
    gates = exp_protocol[0][0]
    dims = (len(gates[0].output_dims()), len(gates))

    Lambdas = []; j = 0
    for i in exp_protocol:
        ledg = Ledger(dims[0])
        ledg.perform_Pauli(i[0], i[1], n_runs = runs, noisy = simulate_noise, wghts = weights, simulator = backend)
        print(ledg.rel_freq()[i[1]]); j += 1
        Lambdas.append(ledg.circuit_eig(vars(ledg.RelRes)[i[1]]))
        print(j, end='\r')

    return Lambdas

#_____________________________________________________________________________________________________________
# 4. Calculate Probabilities

def probs_from_protocol(Lambdas, protocol_matrix, dims, p_errs):
    """
    RETURNS:
        A dict with errors and associated probabilities, e.g.
        probs_from_protocol(Lambdas, ptcl[1], ['I', 'Y']) =
            {'I, I': 0.395, 'I, Y': 0.212, 'Y, I': 0.266, 'Y, Y': 0.128}
    INPUT:
        | Variable        | Type       | Comment
        |_________________|____________|__________________________________________________
        | Lambdas         | list(float)| List of circuit eigenvalues from run_protocol
        | protocol_matrix | list(list) | Rank we demand the protocol matrix to have
        | dims            | list(int)  | Dimensions of dims[0]-qubit dims[1]-gate circuit.
        | p_errs          | str,       | Possible Pauli-errors as str or list of str.
        |                 |  list(str) |
    """
    M = protocol_matrix.copy()
    for i in eachindex(Lambdas):
        if Lambdas[i]<0:
            M[i] = [-1*x for x in M[i]]

    errs = allpauliconfigs(dims, paulis = p_errs)
    probs = matmul(linalg.inv(M),Lambdas)
    return {errs[i]: probs[i] for i in eachindex(errs)}
