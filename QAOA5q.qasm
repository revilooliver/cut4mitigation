OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg meas[5];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
rzz(-0.76459065) q[0],q[3];
rzz(-0.76459065) q[0],q[1];
rzz(-0.76459065) q[1],q[2];
rzz(-0.76459065) q[2],q[3];
rzz(-0.76459065) q[3],q[4];
rx(7.1054101) q[0];
rx(7.1054101) q[1];
rx(7.1054101) q[2];
rx(7.1054101) q[3];
rx(7.1054101) q[4];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
measure q[4] -> meas[4];
