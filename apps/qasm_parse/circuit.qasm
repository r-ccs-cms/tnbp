OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
h q[0];
cx q[0],q[1];
rz(1.570796326948965) q[2];
measure q[1] -> c[1];