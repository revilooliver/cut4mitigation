#!/usr/bin/env python

import numpy as np
import itertools
from qiskit.algorithms.eigensolvers import NumPyEigensolver
from qiskit import IBMQ, Aer
from qiskit import QuantumRegister, ClassicalRegister, execute
from qiskit.circuit import QuantumCircuit, ParameterVector
from qiskit.utils import QuantumInstance
from qiskit.opflow.primitive_ops import Z2Symmetries, PauliOp


# This file includes most of the functions needed in constructing the Hamiltonian and Overlap Matrix 'H' and 'S' in Quantum Subspace Expansion.

# Here we define the functions to operate Paulis:
def SingleMultip(a, b):
    # a, b here should be tuples with coefficient and a Pauli string and 'a' being in the front.

    if a[1] == b[1]:
        ab = (a[0] * b[0], 'I')

    elif a[1] == 'I':
        ab = (a[0] * b[0], b[1])
    elif b[1] == 'I':
        ab = (a[0] * b[0], a[1])

    elif a[1] == 'X':
        if b[1] == 'Y':
            ab = (a[0] * b[0] * (0. + 1.j), 'Z')
        elif b[1] == 'Z':
            ab = (a[0] * b[0] * (0. - 1.j), 'Y')
        else:
            ab = 'Input does not belong to Pauli space'

    elif a[1] == 'Y':
        if b[1] == 'X':
            ab = (a[0] * b[0] * (0. - 1.j), 'Z')
        elif b[1] == 'Z':
            ab = (a[0] * b[0] * (0. + 1.j), 'X')
        else:
            ab = 'Input does not belong to Pauli space'

    elif a[1] == 'Z':
        if b[1] == 'Y':
            ab = (a[0] * b[0] * (0. - 1.j), 'X')
        elif b[1] == 'X':
            ab = (a[0] * b[0] * (0. + 1.j), 'Y')
        else:
            ab = 'Input does not belong to Pauli space'

    else:
        ab = 'Input does not belong to Pauli space'

    return ab


def MultiMultip(a, b):
    # This is an example of multiplication of Pauli string of n qubits.

    a_string = list(a[1])
    b_string = list(b[1])

    num = 1
    string = ''
    for i in range(len(a_string)):
        intermed = SingleMultip((1. + 0.j, a_string[i]), (1. + 0.j, b_string[i]))
        num = num * intermed[0]
        string += intermed[1]

    num = num * a[0] * b[0]
    ab = (num, string)

    return ab


def MeasureCircuit(Sum_Op):
    # Determine how many commute groups are in the SummedOp
    num_terms = len(Sum_Op)

    # Find the Paulis with least number of I in it(shoud be 0).
    # The problem is here. In this case, not not all the subgroups have a term with no I in it. So we have to loop over all
    # of the terms to construct a Pauli string.
    Pauli = ''

    for i in range(Sum_Op.num_qubits):
        intermed = []
        for j in range(num_terms):
            intermed.append(Sum_Op[j].to_label()[i])
        print(intermed)
        if 'X' in intermed:
            Pauli += 'X'
        elif 'Y' in intermed:
            Pauli += 'Y'
        else:
            Pauli += 'Z'

    if len(Pauli) != Sum_Op.num_qubits:
        raise Exception('The length does not equal, traverse has problem.')

    Pauli_string = Pauli[::-1]  # This has reversed the order.
    # Now Pauli_string is the target that we should use to construct the measurement circuit.

    qc = QuantumCircuit(Sum_Op.num_qubits, Sum_Op.num_qubits)
    qc.barrier()

    for i in range(Sum_Op.num_qubits):
        if Pauli_string[i] == 'X':
            qc.u(np.pi / 2, 0, np.pi, i)
        if Pauli_string[i] == 'Y':
            qc.u(np.pi / 2, 0, np.pi / 2, i)
        else:
            None

    qc.measure(range(Sum_Op.num_qubits), range(Sum_Op.num_qubits));

    return qc


def SvMeasureCircuit(Sum_Op):
    # Determine how many terms are in the SummedOp
    num_terms = len(Sum_Op)

    # Find the Paulis with least number of I in it(shoud be 0).
    # The problem is here. In this case, not not all the subgroups have a term with no I in it. So we have to loop over all
    # of the terms to construct a Pauli string.
    Pauli = ''

    for i in range(Sum_Op.num_qubits):
        intermed = []
        for j in range(num_terms):
            intermed.append(Sum_Op[j].primitive.to_label()[i])
        # print(intermed)
        if 'X' in intermed:
            Pauli += 'X'
        elif 'Y' in intermed:
            Pauli += 'Y'
        else:
            Pauli += 'Z'

    if len(Pauli) != Sum_Op.num_qubits:
        raise Exception('The length does not equal, traverse has problem.')

    Pauli_string = Pauli[::-1]  # This has reversed the order.
    # Now Pauli_string is the target that we should use to construct the measurement circuit.

    qc = QuantumCircuit(Sum_Op.num_qubits, Sum_Op.num_qubits)

    for i in range(Sum_Op.num_qubits):
        if Pauli_string[i] == 'X':
            qc.u(np.pi / 2, 0, np.pi, i)
        if Pauli_string[i] == 'Y':
            qc.u(np.pi / 2, 0, np.pi / 2, i)
        else:
            None

    return qc


def evaluation(d: dict, shots: int, Pauli: str):
    # This Pauli_string is in arbitrary form of I, X, Y and Z.
    # Determine the number of qubits, which is also related to the number of measurement outcomes.
    num_qubits = len(Pauli)
    Pauli_string = ''
    for i in Pauli:
        if i == 'I':
            Pauli_string += i
        else:
            Pauli_string += 'Z'

    def kbits(n):
        result = []
        for k in range(0, n + 1):
            for bits in itertools.combinations(range(n), k):
                s = ['0'] * n
                for bit in bits:
                    s[bit] = '1'
                result.append(''.join(s))
        return result

    # Generate all binary strings of N bits.
    outcomes = kbits(num_qubits)

    def get_from(d: dict, key: str):
        value = 0
        if key in d:
            value = d[key]
        return value

    # Here we compute the expectation value.
    expectation_value = 0
    for i in outcomes:
        intermediate = 0
        for j in range(num_qubits):
            if (Pauli_string[j] == 'Z') and (i[j] == '1'):
                intermediate += 1
            else:
                None

        if (intermediate % 2) == 0:
            expectation_value += get_from(d, i)
        else:
            expectation_value -= get_from(d, i)

    expectation_value = expectation_value / shots

    return expectation_value


def evaluation_IZ(d: dict, Pauli: str):
    # This Pauli_string is in arbitrary form of I and Z.
    if 'X' in Pauli or 'Y' in Pauli:
        raise Exception('This evaluation function is for measurements of terms containing only I and Z.')

    # Determine the number of qubits, which is also related to the number of measurement outcomes.
    num_qubits = len(Pauli)
    Pauli_string = ''
    for i in Pauli:
        if i == 'I':
            Pauli_string += i
        else:
            Pauli_string += 'Z'

    # Generate all binary strings of N bits.
    outcomes = ['0101', '1101', '0111', '1111', '1001', '1011', '0110', '1110', '1010']

    def get_from(d: dict, key: str):
        value = 0
        if key in d:
            value = d[key]
        return value

    # Here we compute the expectation value.
    expectation_value = 0
    total_effective_shots = 0
    for i in outcomes:
        intermediate = 0
        for j in range(num_qubits):
            if (Pauli_string[j] == 'Z') and (i[j] == '1'):
                intermediate += 1
            else:
                None

        if (intermediate % 2) == 0:
            expectation_value += get_from(d, i)
        else:
            expectation_value -= get_from(d, i)

        total_effective_shots += get_from(d, i)

    expectation_value = expectation_value / total_effective_shots
    return expectation_value, total_effective_shots


def Svevaluation(sv, Pauli: str):
    # This Pauli_string is in arbitrary form of I, X, Y and Z.
    # Determine the number of qubits, which is also related to the number of measurement outcomes.
    num_qubits = len(Pauli)
    Pauli_string = ''
    for i in Pauli:
        if i == 'I':
            Pauli_string += i
        else:
            Pauli_string += 'Z'

    def kbits(n):
        result = []
        for k in range(0, n + 1):
            for bits in itertools.combinations(range(n), k):
                s = ['0'] * n
                for bit in bits:
                    s[bit] = '1'
                result.append(''.join(s))
        return result

    # Generate all binary strings of N bits.
    outcomes = kbits(num_qubits)

    # Here we compute the expectation value.
    expectation_value = 0
    for i in outcomes:
        intermediate = 0
        for j in range(num_qubits):
            if (Pauli_string[j] == 'Z') and (i[j] == '1'):
                intermediate += 1
            else:
                None

        j = int(i, 2)
        if (intermediate % 2) == 0:
            expectation_value += np.real(sv[j] * np.conjugate(sv[j]))
        else:
            expectation_value -= np.real(sv[j] * np.conjugate(sv[j]))

    return expectation_value

