OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
cz q[5],q[3];
rzz(-0.76459065) q[0],q[3];
rzz(-0.76459065) q[0],q[1];
rx(7.1054101) q[0];
rzz(-0.76459065) q[1],q[2];
rx(7.1054101) q[1];
rzz(-0.76459065) q[2],q[3];
rx(7.1054101) q[2];
rzz(-0.76459065) q[3],q[4];
rx(7.1054101) q[4];
cz q[5],q[3];
rx(7.1054101) q[3];
h q[5];
