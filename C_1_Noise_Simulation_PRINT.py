# ___________________________________________________________________________________________
# C_1 PAULI ERROR MODELS AND NOISY REALISATIONS
# ___________________________________________________________________________________________

# ___0.0 Import packages

from C_ACES.C_0_2_Better_Iterators_PRINT import *

# ___________________________________________________________________________________________
# 1. Random Pauli configurations in space or time

# ___1.0 Output all possible Pauli words


def pauliconfigs(l, paulis="IXYZ", noId=False):
    """
    RETURNS:
        All possible Pauli-error configs of length l,
        e.g. pauliconfigs(2, paulis='IZ') = ['II', 'IZ', 'ZI', 'ZZ']

    INPUT:
        |Variable | Type        | Comment
        |_________|_____________|___________________________________________________________
        | l       | int         | Length of the output Pauli-error configurations.
        | paulis  | str,        | Possible choice of Pauli-errors as str or list of str.
        |         |  list(str)  |   Default: paulis = 'IXYZ'
        | noId    | bool        | Decides if 'III...' will be disregarded.
        |         |             |   Default: noId = False

    """
    if type(paulis) == str:
        configs = tuptostr(list(itertools.product(paulis, repeat=l)))
    elif type(paulis) == list:
        configs = tuptostr(list(itertools.product(paulis, repeat=l)))

    if (noId == True) and ("I" in paulis):
        configs.remove("I" * l)
        return configs
    else:
        return configs
    del configs


# ___1.1 Output all length = dim[1] configs of dim[0]-qubit Puali words


def allpauliconfigs(dims, paulis="IXYZ"):
    """
    RETURNS:
        All possible Pauli-error configs for dims (space, time).

    INPUT:
        |Variable | Type        | Comment
        |_________|_____________|___________________________________________________________
        | dims    | list(int)   | Dimensions of auli config in space (qubits) and time.
        | paulis  | str,        | Possible choice of Pauli-errors as str or list of str.
        |         |  list(str)  |   Default: paulis= 'IXYZ'
    """
    n_confs = dims[0] * dims[1]

    t_configs = list(itertools.product(paulis, repeat=dims[1]))
    s_configs = list(map("".join, list(itertools.product(paulis, repeat=dims[0]))))

    if dims[0] == dims[1] == 1:
        s_t_configs = paulis
    elif dims[1] == 1:
        s_t_configs = s_configs
    elif dims[0] == 1:
        s_t_configs = t_configs
    else:
        s_t_configs = list(
            itertools.product(s_configs, repeat=dims[1])
        )  # all possible pauli configurations
        #  for s qubits and t gates
    s_t_configs = list(map(str, s_t_configs))
    resconfigs = []
    for i in s_t_configs:
        pauliconf = []
        for j in i:
            if j.isalpha() or (j == ",") or (j == " "):
                pauliconf.append(j)
        resconfigs.append("".join(pauliconf))
    return resconfigs


# ___1.2 Output random Pauli word according to prob. dist.


def randpauliconfig(dims, paulis="IXYZ", weights=None):
    """
    RETURNS:
        Pauli-error words of dims (space, time) picked out of all possible configs according to
        weights. E.g. randpauliconfig((2,2), paulis=['I', 'Z']) = ('ZZ', 'ZI')

    INPUT:
        |Variable | Type        | Comment
        |_________|_____________|___________________________________________________________
        | dims    | list(int)   | Dimensions of auli config in space (qubits) and time.
        | paulis  | str,        | Possible choice of Pauli-errors as str or list of str.
        |         |  list(str)  |   Default: paulis= 'IXYZ'
        | weights | list(float) | Probability distribution over all possible output paulis.
        |         |             |   Default: weights = None (-> gives uniform prob. dist.)
        |         | __.Weights  | Class with pauli configs for dims as attributes
        |         |             |   E.g. getattr(W, 'Y, Y') = 0.7
        |         | dict        | Dict with pauli config for dims as keyssign
        |         |             |   E.g. W['Y, Y'] = 0.7
        | noId    | bool        | Decides if 'III...' will be in disregarded.
        |         |             |   Default: noId = False

    """
    n_confs = dims[0] * dims[1]
    if weights == None:  # Default: pick uniformly
        res = ""
        for s in range(0, dims[0]):
            res += random.choice(paulis)
        for t in range(1, dims[1]):
            res += ", "
            for s in range(0, dims[0]):
                res += random.choice(paulis)
    elif type(weights) == Weights:
        res = random.choices(
            list(vars(weights.Values).keys()), list(vars(weights.Values).values())
        )[0]
    elif type(weights) == dict:
        res = random.choices(list(weights.keys()), list(weights.values()))[0]
    return res if len(tuple(res.split(", "))) == 1 else res.split(", ")


def randpauliword(l, paulis="IXYZ", noId=False):
    """
    RETURNS:
        Uniformly selected Pauli word of length l picked out of all possible configs.
        E.g. randerrorconfig(3) = 'IXZ'

    INPUT:
        |Variable | Type        | Comment
        |_________|_____________|___________________________________________________________
        | l       | int         | Length of the output Pauli-error configuration.
        | paulis  | str,        | Possible choice of Pauli-errors as str or list of str.
        |         |  list(str)  |   Default: paulis = 'IXYZ'
        | weights | list(float) | Probability distribution over all possible output paulis.
        |         |             |   Default: weights = None (-> gives uniform prob. dist.)
        | noId    | bool        | Decides if configuration 'III...' will be disregarded.
        |         |             |   Default: noId = False

    """
    weights = [1 / (len(paulis) ** l)] * (len(paulis) ** l)  # Default: pick uniformly

    if "I" not in paulis:  # Identity must be in paulis
        if type(paulis) == str:
            paulis = "I" + paulis
        else:
            paulis.insert(0, "I")
    elif paulis[0] != "I":
        if type(paulis) == str:
            paulis = "I" + paulis.replace("I", "")
        else:
            paulis = paulis.remove(paulis.index("I"), paulis)
            paulis.insert(0, "I")

    if type(paulis) == str:
        configs = tuptostr(list(itertools.product(paulis, repeat=l)))
    elif type(paulis) == list:
        configs = list(itertools.product(paulis, repeat=l))

    if noId == True:
        res = random.choices(configs, weights=weights, k=1)[0]
        while ("X" not in res) and ("Y" not in res) and ("Z" not in res):
            res = random.choices(configs, weights=weights, k=1)[0]
        return res
    else:
        return random.choices(configs, weights=weights, k=1)[0]


# ___1.3 Weights class allows simulation of the noise model


class Weights:
    """
        Define error model through setting up a probability distribution for Pauli errors.

    ATTRIBUTES:
        | Name           |
        |________________|__________________________________________________________________________________
        | .Values        | dict with probabilities for errors
        |                |      vars(W.Values) = {'I': 0.3, 'Y': 0.5, 'X': 0.2}
        | .dims          | dimensions of the error model (number of qubits, number of gates)
        |                |      W.dims = (1, 1)
        | .CompleteDist  | Bool: True if probabilities sum up to 1, False if not
    """

    def __init__(self, dims, paulis="IXYZ"):
        self.Values = type("Values", (), {})()
        self.dims = dims
        self.CompleteDist = False

    def setweight(self, config, weight):
        """
            Defines error model by associating errors with probabilities

        RETURNS:
            Print out of list with errors with warning if probs no sum up to 1
            I :  0.3
            Y :  0.5
            X :  0.2

        INPUT:
            |Variable    | Type           | Comment
            |____________|________________|____________________________________________________
            | config     | list(str)      | List of error configs ['XX, YZ', 'XI, XZ', ...]
            | weight     | list(float)    | List of error probabilities [0.2, 0.1, ...]
        """
        # Allow different finput formats
        if (
            type(config) == list and type(weight) == list
        ):  # List of errors and probabilities
            assert len(config) == len(
                weight
            ), "ERROR: number of weights and configs not equal."
            for i in eachindex(weight):
                setattr(self.Values, config[i], weight[i])
        elif (
            type(config) == str and type(weight) == float
        ):  # Single error with probability
            assert (
                len(config.split(", ")) == self.dims[1]
            ), "ERROR: config does not fit dimensions."
            assert (
                len(config.split(", ")[0]) == self.dims[0]
            ), "ERROR: config does not fit dimensions."
            setattr(self.Values, config, weight)
        # Check and warn, if the distribution is complete, i.e. sums up to 1
        if sum(list(vars(self.Values).values())) < 1:
            self.CompleteDist = False
            print("Sum of weights < 1")
        elif sum(list(vars(self.Values).values())) > 1:
            self.CompleteDist = False
            print("Sum of weights > 1")
        else:  # _Prob dist completion marker
            self.CompleteDist = True  #     self.CompleteDist = True
        for i in vars(self.Values).keys():
            print(i, ": ", vars(self.Values)[i])

    def fill_weights(self, paulis="IXYZ", uniform=True):
        """
        Completes prob dist to 1 if only specific error probs are defined, either uniformly or
        by giving identity the remaining prob.

        INPUT:
            |Variable    | Type           | Comment
            |____________|________________|________________________________________________________________
            | config     | list(str)      | List of error configs ['XX, YZ', 'XI, XZ', ...]
            | weight     | list(float)    | List of error probabilities [0.2, 0.1, ...]
            | uniform    | Bool           | Completes prob dist by either
            |            |                |     uniform = True spreading remainder over rest of error space
            |            |                |             = False setting Id to remainder
            | paulis     | str, list(str) | Possible choice of Pauli-errors as str or list of str.
            |            |                |     Default: paulis = 'IXYZ'
        """
        ps = paulis
        confs = allpauliconfigs(
            self.dims, paulis=ps
        )  # Returns set of all Pauli configs for dims.

        if (len(vars(self.Values)) == 0) and (
            uniform == True
        ):  # If we didn't use self.setweight, fill_weights
            self.setweight(
                confs, [1 / len(confs)] * len(confs)
            )  #   creates a uniform Pauli config distribution.
        elif uniform == True:
            defined_errs = list(vars(self.Values).keys())
            remainingconfs = [e for e in confs if e not in defined_errs]
            counterprobability = 1 - sum(list(vars(self.Values).values()))
            self.setweight(
                remainingconfs,
                [counterprobability / len(remainingconfs)] * len(remainingconfs),
            )
        elif uniform == False:
            counterprobability = 1 - sum(list(vars(self.Values).values()))
            Idim = "I" * self.dims[0] + (", " + "I" * self.dims[0]) * (self.dims[1] - 1)
            self.setweight([Idim], [counterprobability])


# Code examples Weghts class
# W1 = Weights((2, 1))
# W1.setweight(['IX', 'YY', 'XZ', 'II'],
#                [0.3, 0.1, 0.2, 0.4])

# W2 = Weights((1, 2))
# W2.setweight(['I, X', 'Y, Y', 'X, Z', 'I, I'],
#                [0.3, 0.1, 0.2, 0.4])

# W3 = Weights((1,2)); W3.setweight('X, Y', 0.3)
# W3.fill_weights(uniform = False)
# W3 = Weights((1,1)); W3.setweight('X', 0.3)
# W3.fill_weights(uniform = True)

# ___1.3 Pick eigenstates of Pauliwords


def random_eigstate(pauliword, qiskit_order=False):
    """
    RETURNS:
        Tuple (QuantumCircuit, List) containing a randomly chosen eigenstate for each element of the input Pauli word.
        random_eigstate('IZ') = (<qiskit.circuit.quantumcircuit.QuantumCircuit at 0x7fcd10ae9280>, [[0, -1], [0, 1]])
              ┌──────────────────┐
         q_0: ┤ Initialize(0,-1) ├
              ├─────────────────┬┘
         q_1: ┤ Initialize(0,1) ├─
              └─────────────────┘

    INPUT:
        |Variable   | Type  | Comment
        |___________|_______|_______________________________________________________________
        | pauliword | str   | String of Paulis, e. g. 'XIY'
    """
    pauliword = pauliword[1::] if (pauliword[0] in ["+", "-"]) else pauliword
    qc = QuantumCircuit(len(pauliword))  # Create a quantum circuit
    reg = []
    if qiskit_order == False:
        indices = list(range(0, len(pauliword)))
    else:
        indices = list(range(0, len(pauliword)))[::-1]
    for i in indices:
        if pauliword[i] == "I":
            x = random.choices([[1, 0], [0, 1]])[0]  # Chose at random from |0>, |1>
            reg.append(x)
        elif pauliword[i] == "Z":
            x = random.choices([[1, 0], [0, -1]])[0]  # Chose at random from |0>, |1>
            reg.append(x)
        elif pauliword[i] == "X":
            x = random.choices(
                [[1 / sqrt(2), 1 / sqrt(2)], [1 / sqrt(2), -1 / sqrt(2)]]
            )[
                0
            ]  # Chose at random from |0>, |1>
            reg.append(x)
        elif pauliword[i] == "Y":
            x = random.choices(
                [[1 / sqrt(2), 1j / sqrt(2)], [1 / sqrt(2), -1j / sqrt(2)]]
            )[
                0
            ]  # Chose at random from |0>, |1>
            reg.append(x)
    for i in range(0, len(reg)):
        qc.initialize(reg[i], i)
    return qc, reg
    del qc, reg, x


def all_eigenstates(pauliword, qiskit_ord=True):
    """
    RETURNS:
        All eigenstates of a given Pauli word as symbols:
        I: '<', '>'     X: '+', '-'     Y: 'u', 'd'     Z: '0', '1'
        all_eigenstates('XY') = ['u+', 'u-', 'd+', 'd-']

    INPUT:
        |Variable   | Type  | Comment
        |___________|_______|_______________________________________________________________
        | pauliword | str   | String of Paulis, e. g. 'XIY'
    """
    res = []
    if qiskit_ord == True:
        pauliword = pauliword[::-1]
    for i in pauliword:
        if i == "I":
            res.append("<>")
        elif i == "X":
            res.append("+-")
        elif i == "Y":
            res.append("ud")
        elif i == "Z":
            res.append("01")
    if len(pauliword) == 1:
        return [i for i in res[0]]
    else:
        confs_res = tuptostr(list(itertools.product(res[0], res[1])))
        for i in res[2::]:
            confs_res = tuptostr(list(itertools.product(confs_res, i)))
        return confs_res


# ___________________________________________________________________________________________
# 2. Random Clifford list


def randclifflist(dims):
    """
    RETURNS:
        List of t = dims[1] uniformly picked Clifford gates on s = dims[0] qubits,
        e.g. randclifflist((2,1)) = [Clifford(array([...]))]
    """
    gates = []
    for i in range(0, dims[1]):
        gates.append(qinf.random_clifford(dims[0]))
    return gates


# ___________________________________________________________________________________________
# 3. Pauli noisy realisations of circuits from Clifford gates

# ___3.0 Turn List of CLiffors into QuantumCircuit


def gates_circ(gates):
    """
    RETURNS:
        QuantumCircuit of Clifford gates.
        gates_circ(randclifflist((1,2))).draw() =
           ┌──────────┐┌──────────┐
        q: ┤ Clifford ├┤ Clifford ├
           └──────────┘└──────────┘
    INPUT:
        |Variable | Type          | Comment
        |_________|_______________|________________________________________________________
        | gates   | list(Clifford)| List of Clifford gates defining a cicuit.
    """
    dims = (len(gates[0].output_dims()), len(gates))
    circ = QuantumCircuit(dims[0])
    for i in eachindex(gates):
        gtilde = gates[i]
        circ.append(gtilde, range(0, dims[0]))
    return circ


# ___3.1 Noisy realization of circuit as list of Cliffords


def noisy_real(gates, paulis="IXYZ", weights=None, asCliffs=False):
    """
    RETURNS:
        Returns quantum circuit with acc. to weights randomly picked Pauli errors before
        each gate in 'gates':
        noisy_real(randclifflist((1,2))).draw() =
           ┌───┐┌──────────┐┌───┐┌──────────┐
        q: ┤ Y ├┤ Clifford ├┤ X ├┤ Clifford ├
           └───┘└──────────┘└───┘└──────────┘

    INPUT:
        |Variable | Type          | Comment
        |_________|_______________|________________________________________________________
        | gates   | list(Clifford)| List of Clifford gates defining a cicuit.
        | paulis  | str,          | Possible choice of Pauli-errors as str or list of str.
        |         |  list(str)    |   Default: paulis= 'IXYZ'
        | weights | list(float)   | Probability distribution over all possible output paulis.
        |         |               |   Default: weights = None (-> gives uniform prob. dist.)

    """
    dims = (len(gates[0].output_dims()), len(gates))
    paulierrors = randpauliconfig(dims, paulis=paulis, weights=weights)

    qc = QuantumCircuit(dims[0])
    noisy_gates = []

    if type(paulierrors) == str:
        qc.append(qinf.Pauli(paulierrors), range(0, dims[0]))
        qc.append(gates[0], range(0, dims[0]))
        if asCliffs == True:
            noisy_gates.append(qinf.Clifford(qc))
    else:
        for i in eachindex(gates):
            qc.append(qinf.Pauli(paulierrors[i]), range(0, dims[0]))
            qc.append(gates[i], range(0, dims[0]))
            if asCliffs == True:
                noisy_gates.append(qinf.Clifford(qc))
                qc = QuantumCircuit(dims[0])
    if asCliffs == True:
        return noisy_gates
    else:
        return qc
