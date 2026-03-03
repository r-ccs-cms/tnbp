from qiskit import QuantumCircuit
from qiskit.qasm2 import dumps
# from qiskit.qasm3 import dumps

qc = QuantumCircuit(3,3)
qc.h(0)
qc.cx(0,1)
qc.rz(3.14159265389793/2, 2)
qc.measure(1,1)

# OpenQASM 2.0
qasm_str = dumps(qc)

# dump file
with open("circuit.qasm", "w") as f:
    f.write(qasm_str)
print("wrote circuit.qasm")
