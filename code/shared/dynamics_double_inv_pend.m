function [x_dot, y] = dynamics_double_inv_pend(t, x, u, params)
% --Double Inverted Pendulum Dynamics--

% unpack parameters
mc = params.mc;
m1 = params.m1;
m2 = params.m2;
J1 = params.J1;
J2 = params.J2;
l1 = params.l1;
l2 = params.l2;
L1 = params.L1;
L2 = params.L2;
cc = params.cc;
c1 = params.c1;
c2 = params.c2;
g = params.g;

% unpack states
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);
x6 = x(6);

x1dot = x2;
x3dot = x4;
x5dot = x6;

qddot = compute_q_ddot(x1,x2,x3,x4,x5,x6,u,mc,m1,m2,J1,J2,l1,l2,L1,L2,cc,c1,c2,g);
x2dot = qddot(1);
x4dot = qddot(2);
x6dot = qddot(3);

x_dot = [x1dot; x2dot; x3dot; x4dot; x5dot; x6dot];
y = x_dot;

end