OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
cz q[0],q[1];
cz q[1],q[2];
cz q[2],q[3];