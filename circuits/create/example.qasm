OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[3];
cz qubits[0],qubits[1];
h qubits[2];
h qubits[0];
cz qubits[1],qubits[2];