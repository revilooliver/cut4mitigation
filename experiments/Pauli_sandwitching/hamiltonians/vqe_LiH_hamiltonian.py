import numpy as np
#standard Qiskit libraries
from qiskit import QuantumCircuit, transpile, Aer, IBMQ
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from qiskit.providers.aer import QasmSimulator
from qiskit_nature.drivers.second_quantization import PySCFDriver
from qiskit_nature.transformers.second_quantization.electronic import FreezeCoreTransformer
from qiskit_nature.mappers.second_quantization import ParityMapper
from qiskit_nature.converters.second_quantization.qubit_converter import QubitConverter
from qiskit_nature.circuit.library import HartreeFock
from qiskit_nature.problems.second_quantization.electronic import ElectronicStructureProblem
from qiskit.circuit import Parameter, QuantumCircuit, QuantumRegister
from qiskit_nature.algorithms.ground_state_solvers.minimum_eigensolver_factories import NumPyMinimumEigensolverFactory
from qiskit_nature.algorithms.ground_state_solvers import GroundStateEigensolver

def exact_diagonalizer(problem, converter):
    solver = NumPyMinimumEigensolverFactory()
    calc = GroundStateEigensolver(converter, solver)
    result = calc.solve(problem)
    return result

def main():
    molecule = 'Li 0.0 0.0 0.0; H 0.0 0.0 0.0'
    driver = PySCFDriver(atom=molecule)
    print("before run")
    qmolecule = driver.run()
    print("finished run")
    fct = [FreezeCoreTransformer(freeze_core=True, remove_orbitals=[3,4])]

    problem = ElectronicStructureProblem(driver, fct)
    second_q_ops = problem.second_q_ops()

    # Hamiltonian
    main_op = second_q_ops[0]

    converter = QubitConverter(mapper=ParityMapper(), two_qubit_reduction=True)

    # Mapping Fermions to Qubits
    #num_particles = (problem.molecule_data_transformed.num_alpha,
                 #problem.molecule_data_transformed.num_beta)
    num_particles = problem.num_particles
    qubit_op = converter.convert(main_op, num_particles = problem.num_particles)


    num_spin_orbitals = problem.num_spin_orbitals
    init_state = HartreeFock(num_spin_orbitals, num_particles, converter)

    # Parameters for q-UCC antatze
    num_particles = problem.num_particles
    num_spin_orbitals = problem.num_spin_orbitals
    n = qubit_op.num_qubits
    qc = QuantumCircuit(qubit_op.num_qubits)

    #the variational parameter
    p=1
    for i in range(n):
        theta = Parameter(f"ry_theta{p}" )
        qc.ry(theta, i)
        p += 1
    qubit_label = 0
    # qc.ry(theta, range(n))
    #qc.rz(theta, range(n))
    for i in range(n-1):
        qc.cz(i, i+1)
    for i in range(n):
        theta = Parameter(f"ry_theta{p}" )
        qc.ry(theta, i)
        p += 1
    #qc.rz(theta, range(n))

    # Add the initial state
    ansatz = qc
    ansatz.compose(init_state, front=True, inplace=True)

    result_exact = exact_diagonalizer(problem, converter)
    exact_energy = np.real(result_exact.eigenenergies[0])
    print("Exact electronic energy", exact_energy)
    print(result_exact)
    return None


if __name__ == "__main__":
    main()