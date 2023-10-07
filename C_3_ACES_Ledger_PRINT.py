#_____________________________________________________________________________________________________________
# C_3 ACES CIRCUITS AND LEDGER CLASS
#_____________________________________________________________________________________________________________

from C_ACES.C_2_GtPt_PRINT import *

#_____________________________________________________________________________________________________________
# 1. Run State Through Cirucit

def ACES_circ(gates, pauliword, noisy = False, weights = None, qiskit_ord = True):
    """
    RETURNS:
        Tuple (QuantumCircuit, list, str) G-twisted Pauli twirled circuit, the initial
        Pauli eigenstate state and the Pauli measurement basis.
        E.g. ACES_circ(randclifflist((1, 1)), 'Z') =
        (<qiskit.circuit...>,  [[0, -1]], 'Z')
             ┌──────────────────┐┌───┐┌──────────┐┌───┐ ░ ┌─┐
          q: ┤ Initialize(0,-1) ├┤ Y ├┤ Clifford ├┤ X ├─░─┤M├
             └──────────────────┘└───┘└──────────┘└───┘ ░ └╥┘
        c: 1/══════════════════════════════════════════════╩═
                                                           0

    INPUT:
        |Variable    | Type           | Comment
        |____________|________________|___________________________________________________________
        | gates      | list(Clifford) | List of Clifford gates defining a cicuit.
        | pauliword  | str            | String of Paulis, e. g. 'XIY'. RECALL: flipped qubit order
        | noisy      | bool           | Introduces noise to input gates
        |            |                |     default: noisy = False
        | weights    | list(float)    | Prob. distre. for Pauli errors occuring in circuit
        |            |                |     default: weights = None -> uniform prob. dist.
        | qiskit_ord | bool           | If True Pauli(pauliword) in reverse order on qubits.
        |            |                |     default: qiskit_ord = True
        |            |                |                                             ┌───┐
        |            |                |                                        q_0: ┤ X ├
        |            |                |              therefore: Pauli('YX') ->      ├───┤
        |            |                |                                        q_1: ┤ Y ├
        |            |                |                                             └───┘

    """
    dims = (len(gates[0].output_dims()), len(gates))
    pauli_after_gates = gates_evolve(gates, pauliword, no_minus=True, qiskit_order=qiskit_ord)

    qc = QuantumCircuit(dims[0], dims[0])
    ran_eigstat = random_eigstate(pauliword, qiskit_order = qiskit_ord)#;ran_eigstat[0].draw()
    qc.append(ran_eigstat[0], list(range(0, dims[0])))                          # Prepare random eigenstate of Pauli word

    if noisy == False:
        qc.append(GtPt(gates_circ(gates)), range(0, dims[0]))                   # Perform GtPt of circuit
    else:
        w = weights
        qc.append(GtPt(noisy_real(gates, weights = w, asCliffs = True)), range(0, dims[0]))
    #qc.draw()
    qc = qc.decompose()
    qc.barrier()

    for i in range(0, dims[0]):                                                 # Perform measurement in output Pauli base
        if (pauli_after_gates[0][i] == 'I') ^ (pauli_after_gates[0][i] == 'Z'):
            qc.measure(i, i)
        elif pauli_after_gates[0][i] == 'X':
            qc.h(i)
            qc.measure(i, i)
        elif pauli_after_gates[0][i] == 'Y':
            qc.sdg(i); qc.h(i)
            qc.measure(i, i)

    return qc, ran_eigstat[1], pauli_after_gates[1]

#_____________________________________________________________________________________________________________
# 2. Construct Ledger for Experiment

from math import sqrt
from qiskit import IBMQ, Aer

class Ledger:
    """
    ATTRIBUTES:
        | Name           |
        |________________|__________________________________________________________________________________
        | .Arch          | Ledger keepoing record of the outcome of each experiment
        |    .X, .Y, ... | Dicts with the '+' or '-' experiments as lists of '0's and '1's
        | .Res           | Ledger with absolute frequencies of outcomes for each initially chosen Pauli word
        |    .X, .Y, ... | Dicts with the '+' or '-' experiments, each with {'0': ...,'1': ... }
        | .RelRes        | Ledger with relative frequencies of outcomes for each initially chosen Pauli word
        |    .X, .Y, ... |
        | .CircEigVals   | Ledger containing circuit eigenvalues
        |    .X, .Y, ... |

    """
    def __init__(self, l, noIdentity = False):
        """
        INIT:
            Constructs attributes self.Arch, .Res, .RelRes, .CircEigVals containing dict {'+':[], '-':[]} for each Pauli config.
        """
        self.Res = type("Results", (),{})()
        self.Arch = type("Archive", (),{})()
        self.CircEigVals = type("CircuitEigenvalues", (),{})()
        self.RelRes = type("Results", (),{})()
        self.sgn = type("Signs", (),{})()
        configs = pauliconfigs(l, noId=noIdentity)
        self.len = l
        for i in configs:
            setattr(self.Arch, i, {v: [] for v in all_eigenstates(i)})
            setattr(self.sgn, i, 0)

    def new_p_or_m(self, inistates, pauliword, qiskit_ord = True):
        """
        RETURNS:
            String, cathegorizing single qubit Pauli eigenstate as
            '<' or '>'  for 'I'     'u' or 'd'  for 'Y'
            '+' or '-'  for 'X'     '0' or '1'  for 'Z'


        INPUT:
            |Variable    | Type           | Comment
            |____________|________________|____________________________________________________
            | inistates  | list(list)     | Single qubit pauli eigenstate
            | pauliword  | str            | String of Paulis, e. g. 'XIY'
            | qiskit_ord | bool           | If True Pauli(pauliword) in reverse order on qubits.
            |            |                |     default: qiskit_ord = True
        """
        if qiskit_ord == True:
            pauliword = pauliword[::-1]
        res = ''
        for i in eachindex(inistates):
            if pauliword[i] == 'I':
                if inistates[i] == [1,0]:
                    res += '<'
                    #print(key[i])
                elif inistates[i] == [0,1]:
                    res += '>'
            elif pauliword[i] == 'X':
                if inistates[i] == [0.7071067811865475, 0.7071067811865475]:
                    res += '+'
                    #print(key[i])
                elif inistates[i] == [0.7071067811865475, -0.7071067811865475]:
                    res += '-'
            elif pauliword[i] == 'Y':
                if inistates[i] == [0.7071067811865475, 0.7071067811865475j]:
                    res += 'u'
                    #print(key[i])
                elif inistates[i] == [0.7071067811865475, -0.7071067811865475j]:
                    res += 'd'
            if pauliword[i] == 'Z':
                if inistates[i] == [1,0]:
                    res += '0'
                    #print(key[i])
                elif inistates[i] == [0,-1]:
                    res += '1'
        return res

    def outputkey(self, key, pauli_after_gates, qiskit_ord = True):
        """
        RETURNS:
            String, labeling mesurement outcome state accoirding to output Pauli, E.g.
            self.outputkey('01', '-IY', qiskit_ord = True) = '1<'

        INPUT:
            |Variable    | Type           | Comment
            |____________|________________|____________________________________________________
            | inistates  | list(list)     | Single qubit pauli eigenstate
            | pauliword  | str            | String of Paulis, e. g. 'XIY'
            | qiskit_ord | bool           | If True Pauli(pauliword) in reverse order on qubits.
            |            |                |     default: qiskit_ord = True
        """
        new_key = ''
        after_gates = pauli_after_gates
        if pauli_after_gates[0] == '-':                                             # Remove sign from output Pauli
            after_gates = pauli_after_gates[1::]
        for i in eachindex(after_gates):
            if after_gates[i] == 'I':
                if key[i] == '0':
                    new_key += '<'
                    #print(key[i])
                elif key[i] == '1':
                    new_key += '>'
            elif after_gates[i] == 'X':
                if key[i] == '0':
                    new_key += '+'
                    #print(key[i])
                elif key[i] == '1':
                    new_key += '-'
            elif after_gates[i] == 'Y':
                if key[i] == '0':
                    new_key += 'u'
                    #print(key[i])
                elif key[i] == '1':
                    new_key += 'd'
            else:
                new_key += key[i]
        if qiskit_ord == True:
            new_key = new_key[::-1]
        return new_key

    def abs_freq(self):
        """
            Counts absolut number of outputs for each output state label.

        RETURNS:
            Dict with absolute numbers of shots from self.Arch in self.Res
            E.g. self.abs_freq() = {'X': {'+': {'0': 0, '1': 93}, '-': {'0': 91, '1': 0}}, ...}
        """

        configs = list(vars(self.Arch).keys())
        for i in configs:
            setattr(self.Res, i, {k: occurance_counter(v) for k, v in getattr(self.Arch, i).items()})
        return vars(self.Res)
        del configs

    def rel_freq(self, rounded = False):
        """
            Calculates percentage/relative frequencies of output eigenstates for each output Pauli.

        RETURNS:
            Dict with relative frequencies of shots from self.Arch in self.RelRes
            E.g. self.abs_freq() = {'X': {'+': {'0': 0.0, '1': 1.0}, '-': {'0': 1.0, '1': 0.0}}, ...}
        """
        if len(vars(self.Res)) == 0:
            self.abs_freq()
        configs = list(vars(self.Arch).keys())
        for i in eachindex(configs):
            AbsFreq = {v: len(getattr(self.Arch, configs[i])[v]) for v in getattr(self.Arch, configs[i]).keys()}
            if rounded == True:
                out = {v: {m: round((vars(self.Res)[configs[i]][v][m]/AbsFreq[v] if AbsFreq[v]!=0 else 0), ndigits = 2) for m in vars(self.Res)[configs[i]][v].keys()} for v in AbsFreq.keys()}
                setattr(self.RelRes, configs[i], out)
            else:
                out = {v: {m: (vars(self.Res)[configs[i]][v][m]/AbsFreq[v] if AbsFreq[v]!=0 else 0) for m in vars(self.Res)[configs[i]][v].keys()} for v in AbsFreq.keys()}
                setattr(self.RelRes, configs[i], out)
        return vars(self.RelRes)
        del configs, AbsFreq, out

    def perform_Pauli(self, gates, pauliword, n_runs=512, simulator = Aer.get_backend('aer_simulator'), noisy = False, wghts = None):
        """
            Performs ACES experiment with specific input Pauli word and circuit

        INPUT:
            |Variable    | Type           | Comment
            |____________|________________|____________________________________________________
            | gates      | list(Clifford) | List of Clifford gates defining a cicuit.
            | pauliword  | str            | String of Paulis, e. g. 'XIY'
            | n_runs     | int            | Number of runs in experiment
            | simulator  | Qiskit         | Qiskit quantum simulator
            | noisy      | bool           | If = True noise simulated with wghts
            | wghts      | W class        | Noise distribution defined through W class
        """
        dims = (len(gates[0].output_dims()), len(gates))
        for i in range(0, n_runs):
            # Pick random Pauli word
            if noisy == True:
                GtPt_out = ACES_circ(gates, pauliword, noisy=True, weights=wghts)
            else:
                GtPt_out = ACES_circ(gates, pauliword)
            # Run GtPt on Circuit and measure in output Pauli base
            circ = GtPt_out[0]#; circ.draw()
            circ = transpile(circ, simulator)

            # Run and get counts
            result = simulator.run(circ, shots=1).result()
            counts = result.get_counts(circ)
            in_key = self.new_p_or_m(GtPt_out[1], pauliword)
            out_key = self.outputkey(list(counts.keys())[0], GtPt_out[2])
            getattr(self.Arch, pauliword)[in_key].append(out_key)


    def perform_experiment(self, gates, n_runs = 512, simulator = Aer.get_backend('aer_simulator'), noisy = False, wghts = None):
        """
            Performs ACES experiment over whole Pauli group for specific cirucuit

        INPUT:
            |Variable    | Type           | Comment
            |____________|________________|____________________________________________________
            | gates      | list(Clifford) | List of Clifford gates defining a cicuit.
            | n_runs     | int            | Number of runs in experiment
            | simulator  | Qiskit         | Qiskit quantum simulator
            | noisy      | bool           | If = True noise simulated with wghts
            | wghts      | W class        | Noise distribution defined through W class
        """
        dims = (len(gates[0].output_dims()), len(gates))
        for i in vars(self.sgn).keys():
            evopauli = gates_evolve(gates, i)[1]
            if evopauli[0]=='-':
                setattr(self.sgn, i, -1)
            else:
                setattr(self.sgn, i, +1)
        for i in range(0, n_runs):
            # Pick random Pauli word
            pauliword = randpauliconfig((dims[0], 1))
            if noisy == True:
                GtPt_out = ACES_circ(gates, pauliword, noisy=True, weights = wghts)
            else:
                GtPt_out = ACES_circ(gates, pauliword)
            # Run GtPt on Circuit and measure in output Pauli base
            circ = GtPt_out[0]#; circ.draw()
            circ = transpile(circ, simulator)

            # Run and get counts
            result = simulator.run(circ, shots=1).result()
            counts = result.get_counts(circ)
            in_key = self.new_p_or_m(GtPt_out[1], pauliword)
            out_key = self.outputkey(list(counts.keys())[0], GtPt_out[2])
            getattr(self.Arch, pauliword)[in_key].append(out_key)#
        print('Performed ', n_runs, ' runs')

    def sign(self, key):
        """
        RETURNS:
            Eigenvalues of measurement outcome states, e.g. for state labels like '+','-','0' or '1' it returns either +1 or -1.
        """
        if key in ['+', 'u', '0', '<', '>']:
            return +1
        elif key in ['-', 'd', '1']:
            return -1
        else:
            print("""Only Inputs '+','-','0' or '1' allowed!""")

    def keyssign(self, key):
        """
            Maps sign on to the Pauli eigenstates
        """
        return prod(list(map(self.sign, list(key))))

    def circuit_eig(self, sample):
        """
            Calculates circuit eigenvalues for a specific sample associated with an input Pauli word.

        INPUT:
            |Variable    | Type | Comment
            |____________|______|______________________________________________________________
            | sample     | dict | ACES results from a single input Pauli word
            |            |      |   'X': {'+': {'1': 0.693..., '0': 0.306...},
            |            |      |         '-': {'0': 0.705..., '1': 0.294...}}
        """
        res = 0
        for v in sample.keys():
            for w in sample[v].keys():
                res += self.keyssign(v)*self.keyssign(w)*sample[v][w]/(len(sample.keys()))
        return res

    def circuit_eigenvalues(self):
        """
            Calculates circuit eigenvalues for the whole Ledger.
        """
        configs = list(vars(self.Arch).keys())
        global results
        if len(vars(self.RelRes)) == 0:
            results = self.rel_freq()
        else:
            results = vars(self.RelRes)
        for i in results.keys():
            setattr(self.CircEigVals, i, self.circuit_eig(results[i])*getattr(self.sgn, i))
        return vars(self.CircEigVals)
