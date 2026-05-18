%% Double Inverted Pendulum Derivations
%% EOMs
clc;clear;close('all');

% cart pos
syms xc theta1 theta2 l1 l2 L1 L2 real
x1 = xc - l1*sin(theta1);
y1 = l1*cos(theta1);

x2 = xc - L1*sin(theta1) - l2*sin(theta2);
y2 = L1*cos(theta1) + l2*cos(theta2);

% cart vel
syms xc_dot theta1_dot theta2_dot real
dx1 = gradient(x1, [xc; theta1; theta2])'*[xc_dot; theta1_dot; theta2_dot];
dy1 = gradient(y1, [xc; theta1; theta2])'*[xc_dot; theta1_dot; theta2_dot];

dx2 = gradient(x2, [xc; theta1; theta2])'*[xc_dot; theta1_dot; theta2_dot];
dy2 = gradient(y2, [xc; theta1; theta2])'*[xc_dot; theta1_dot; theta2_dot];

v1sq = simplify(dx1^2 + dy1^2);
v2sq = simplify(expand(dx2^2 + dy2^2));

% v2sq_hand = simplify(xc_dot^2 -2*L1*xc_dot*theta1_dot*cos(theta1) - 2*l2*xc_dot*theta2_dot*cos(theta2) + L1^2*theta1_dot^2 + 2*L1*l2*theta1_dot*theta2_dot*cos(theta1)*cos(theta2) + 2*L1*l2*theta1_dot*theta2_dot*sin(theta1)*sin(theta2) + l2^2*theta2_dot^2);

% --Lagrangian Mechanics--
syms g mc m1 m2 J1 J2 real

% potential
U = g*(m1*y1 + m2*y2);

% kinetic
Tc = .5*mc*xc_dot^2;
T1 = .5*m1*v1sq + .5*J1*theta1_dot^2;
T2 = .5*m2*v2sq + .5*J2*theta2_dot^2;
T = Tc + T1 + T2;

% lagrangian
L = simplify(expand(T - U));

% -Compute Euler Lagrange Equations-
syms xc_ddot theta1_ddot theta2_ddot real

% EL 1
dL_dxc = diff(L, xc);
dL_dxcdot = diff(L, xc_dot);
dt_dL_dxcdot = gradient(dL_dxcdot, [xc; theta1; theta2; xc_dot; theta1_dot; theta2_dot])'*[xc_dot; theta1_dot; theta2_dot; xc_ddot; theta1_ddot; theta2_ddot];

% EL 2
dL_dtheta1 = diff(L, theta1);
dL_dtheta1dot = diff(L, theta1_dot);
dt_dL_dtheta1dot = gradient(dL_dtheta1dot, [xc; theta1; theta2; xc_dot; theta1_dot; theta2_dot])'*[xc_dot; theta1_dot; theta2_dot; xc_ddot; theta1_ddot; theta2_ddot];

% EL 3
dL_dtheta2 = diff(L, theta2);
dL_dtheta2dot = diff(L, theta2_dot);
dt_dL_dtheta2dot = gradient(dL_dtheta2dot, [xc; theta1; theta2; xc_dot; theta1_dot; theta2_dot])'*[xc_dot; theta1_dot; theta2_dot; xc_ddot; theta1_ddot; theta2_ddot];

% -Compute ODEs w/ Generalized Forces-
syms F cc c1 c2 real
Q_xc = F - cc*xc_dot;
Q_theta1 = -c1*theta1_dot + c2*(theta2_dot - theta1_dot);
Q_theta2 = -c2*(theta2_dot - theta1_dot);

eom1 = dt_dL_dxcdot - dL_dxc == Q_xc;
eom2 = dt_dL_dtheta1dot - dL_dtheta1 == Q_theta1;
eom3 = dt_dL_dtheta2dot - dL_dtheta2 == Q_theta2;

% M*q_ddot = B
q_ddot_sym = [xc_ddot; theta1_ddot; theta2_ddot];
[M, B] = equationsToMatrix([eom1; eom2; eom3], q_ddot_sym);
M = simplify(M);
B = simplify(B);
q_ddot_eom = M\B;

% define states
vars_old = [xc; xc_dot; theta1; theta1_dot; theta2; theta2_dot; F];
syms x1 x2 x3 x4 x5 x6 u real
vars_new = [x1; x2; x3; x4; x5; x6; u];
q_ddot = subs(q_ddot_eom, vars_old, vars_new);


% matlabFunction(q_ddot, 'File', 'compute_q_ddot', 'Vars', [vars_new; mc; m1; m2; J1; J2; l1; l2; L1; L2; cc; c1; c2; g],'Optimize', false)

% linearized state space symbolic
x1dot = x2;
x2dot = q_ddot(1);
x3dot = x4;
x4dot = q_ddot(2);
x5dot = x6;
x6dot = q_ddot(3);
x_dot = [x1dot; x2dot; x3dot; x4dot; x5dot; x6dot];

A = jacobian(x_dot, [x1; x2; x3; x4; x5; x6]);
B = jacobian(x_dot, u);
C = eye(6);

% matlabFunction(A, B, C, 'File', 'compute_linearized_dynamics', 'Vars', [vars_new; mc; m1; m2; J1; J2; l1; l2; L1; L2; cc; c1; c2; g],'Optimize', false)

% ---Passivity Control/ Feedback Linearization---
x_vars = [x1; x2; x3; x4; x5; x6];
g_x = jacobian(x_dot, u);
h = x1;

ydot = gradient(h, x_vars)'*x_dot;
yddot = gradient(ydot,x_vars)'*x_dot;

beta = 1/(gradient(ydot,x_vars)'*g_x);
alpha = -subs(yddot,u,0)*beta;

beta = simplify(beta);
alpha = simplify(alpha);

% matlabFunction(alpha, beta, 'File', 'compute_alpha_beta', 'Vars', [vars_new; mc; m1; m2; J1; J2; l1; l2; L1; L2; cc; c1; c2; g],'Optimize', false)

E_tot = T + U;
E_current = subs(E_tot, vars_old, vars_new);
Edot_current = simplify(gradient(E_current, x_vars)'*x_dot);

% matlabFunction(E_current, 'File', 'compute_energy', 'Vars', [vars_new; mc; m1; m2; J1; J2; l1; l2; L1; L2; cc; c1; c2; g],'Optimize', false)

syms v
E_pend = subs(E_current, x2, 0);
Edot_pend = gradient(E_pend, x_vars)'*x_dot;
Edot_pend = simplify(subs(Edot_pend, u, alpha + beta*v));
