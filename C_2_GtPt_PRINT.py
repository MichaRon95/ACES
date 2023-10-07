#_____________________________________________________________________________________________________________
# C_2 G-TWISTED PAULI-TWIRLING
#_____________________________________________________________________________________________________________

from C_ACES.C_1_Noise_Simulation_PRINT import *

#_____________________________________________________________________________________________________________
# 1. GtPt

def gates_evolve(gates, pauliword, no_minus = False, qiskit_order = False):
    """
    RETURNS:
        Tuple ofPauli evolved through circuit of Clifford gates, once with option and once without, E.g.
        gates_evolve(randclifflist((4,4)), 'ZYXX', no_minus = True, qiskit_ord = True) = ('YIXY', '-YXIY')
    INPUT:
        |Variable   | Type           | Comment
        |___________|________________|_____________________________________________________
        | gates     | list(Clifford) | List of Clifford gates defining a cicuit.
        | pauliword | str            | String of Paulis, e. g. 'XIY'
    """
    evolved_Pauli = qinf.Pauli(pauliword).evolve(gates[0], frame='s')#; evolved_Pauli
    if len(gates) > 1:
        for i in gates[1::]:
            evolved_Pauli = qinf.Pauli(evolved_Pauli).evolve(i, frame='s')

    after_gates = evolved_Pauli
    if (evolved_Pauli[0] == '-') and (no_minus == True):
        after_gates = evolved_Pauli[1::]
    if qiskit_order == True:
        after_gates = after_gates[::-1]
    return after_gates.to_label(), evolved_Pauli.to_label()

# ___1.0 GtPt single gate

def GtPt_single_gate(G):
    """
    RETURNS:
        Returns G-twisted Pauli-twirled quantum circuit from a single CLifford.
        GtPt_single_gate(randclifflist((2,1))).draw() =
             ┌────────────┐┌───────────┐┌────────────┐
        q_0: ┤0           ├┤0          ├┤0           ├
             │  Pauli(ZZ) ││  Clifford ││  Pauli(IX) │
        q_1: ┤1           ├┤1          ├┤1           ├
             └────────────┘└───────────┘└────────────┘
    INPUT:
        |Variable | Type                | Comment
        |_________|_____________________|__________________________________________________
        | G       |  CircuitInstruction,| Some Clifford gate inside a Circuit
        |         |   QuantumCircuit,   |
        |         |   list(Clifford),   |
        |         |   Clifford          |
    """
    # Possible input types: QuantumCircuit, CircuitInstruction, list(Clifford), Clifford
    if type(G) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        dim = (G.num_qubits, 1)
        cliff = qinf.Clifford(G)
    elif type(G) == qiskit.circuit.quantumcircuitdata.CircuitInstruction:
        cliff = qinf.Clifford(G.operation)
        dim = (len(cliff.output_dims()), 1)
    elif (type(G) == list) and (type(G[0]) == qiskit.quantum_info.operators.symplectic.clifford.Clifford):
        if len(G) > 1:
            print('ERROR: function only takes single gate input')
        elif len(G) == 0:
            print('ERROR: no gate specified')
        else:
            dim = (len(G[0].output_dims()), 1)
            cliff = G[0]
    elif type(G) == qiskit.quantum_info.operators.symplectic.clifford.Clifford:
        dim = (len(G.output_dims()), 1)
        cliff = G
    else:
        print('ERROR: type not matching')

    P_a1 = qinf.Pauli(randpauliconfig(dim))                             # Generate random Pauli word
    GoP_a1 = P_a1.evolve(cliff, frame='s')                              # Evolve Pauli by Gate
    output_circ = QuantumCircuit(dim[0])
    output_circ.append(P_a1, range(0, dim[0]))
    output_circ.append(cliff, range(0, dim[0]))
    output_circ.append(GoP_a1, range(0, dim[0]))

    return output_circ

#_____________________________________________________________________________________________________________
# 2. GtPt Multi Gate Circuit

def GtPt(input_circ):
    """
    RETURNS:
        Returns G-twisted Pauli-twirled quantum circuit of Cliffords.
        GtPt(gates_circ(randclifflist((1,2)))).draw() =
           ┌───┐┌──────────┐┌────┐┌───┐┌──────────┐┌───┐
        q: ┤ X ├┤ Clifford ├┤ -X ├┤ Z ├┤ Clifford ├┤ Y ├
           └───┘└──────────┘└────┘└───┘└──────────┘└───┘

    INPUT:
        |Variable    | Type           | Comment
        |____________|________________|________________________________________________________
        | input_circ | QuantumCircuit | Circuit of Clifford gates.
    """
    if type(input_circ) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        dim = (input_circ[::][0].operation.num_qubits, len(input_circ[::]))
    elif type(input_circ[0]) == qiskit.quantum_info.operators.symplectic.clifford.Clifford:
        dim =  (len(input_circ[0].output_dims()), len(input_circ))
    twirled_circuit = QuantumCircuit(dim[0])
    for i in eachindex(input_circ[::]):
        twirled_circuit.append(GtPt_single_gate(input_circ[::][i]), range(0, dim[0]))

    return twirled_circuit.decompose()
