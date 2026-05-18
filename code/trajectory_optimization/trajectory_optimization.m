%% Trajectory Optimization
clc;clear;close('all');

% define constants
mc = 400e-3; params.mc = mc;
m1 = 0.15; params.m1 = m1;
m2 = 0.15; params.m2 = m2;
L1 = 0.5; params.L1 = L1;
L2 = 0.5; params.L2 = L2;
l1 = L1/2; params.l1 = l1;
l2 = L2/2; params.l2 = l2;
J1 = m1*L1^2/12; params.J1 = J1;
J2 = m2*L2^2/12; params.J2 = J2;
cc = 0.0; params.cc = cc;
c1 = 0.0; params.c1 = c1;
c2 = 0.0; params.c2 = c2;
g = 9.81; params.g = g;

% eq pts
x0 = [0.0; 0.0; pi; 0.0; pi; 0]; % down-down
xstar = zeros(6,1); % up-up

n = length(x0);

N = 151;
dt = 0.04;
t = 0:dt:(N-1)*dt;

% Initial Condition Trajectory
X0_ic = zeros(n, N);
U0_ic = zeros(1, N);
% Parameters for the "pump" guess
cart_amplitude = 2.0; % Guess the cart swings 2 meters left/right
num_swings = 1.5;     % 1.5 full back-and-forth cycles

for k = 1:N
    % Fraction goes from 0 to 1 over the trajectory
    fraction = (k-1) / (N-1);
    
    % 1. Cart Guess: A sine wave that naturally dampens out to 0 at the end
    xc_guess = cart_amplitude * sin(2 * pi * num_swings * fraction) * (1 - fraction);
    
    % 2. Pendulum Guess: Sweep linearly from down (pi) to up (0)
    % (Assuming down-down is [pi; pi] and up-up is [0; 0])
    theta1_guess = pi * (1 - fraction);
    theta2_guess = pi * (1 - fraction);
    
    % We leave velocities at 0 for the guess to keep things simple
    X0_ic(:, k) = [xc_guess; 0; theta1_guess; 0; theta2_guess; 0];
end

z0 = [X0_ic(:); U0_ic(:)];


% Inequality Constraints
lb = -inf(size(z0));
ub = inf(size(z0));
lb_x = reshape(lb(1:n*N), n, N); % n x N
ub_x = reshape(ub(1:n*N), n, N);
lb_u = reshape(lb(n*N+1:end), 1, N); % m x N
ub_u = reshape(ub(n*N+1:end), 1, N);

for k = 1:N
    lb_x(1,k) = -8;
    ub_x(1,k) = 8;

    lb_u(1,k) = -50;
    ub_u(1,k) = 50;
end
lb = [lb_x(:); lb_u(:)];
ub = [ub_x(:); ub_u(:)];

options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter', 'MaxFunctionEvaluations', 100000, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
Z_opt = fmincon(@(z)cost_traj(z, n, N), z0, [], [], [], [], lb, ub, @(z)nonlcon(z,n,N,dt,params,x0,xstar), options);

X_star_coarse = reshape(Z_opt(1:n*N), n, N);
U_star_coarse = Z_opt(n*N+1:end)';
save('reference_traj_coarse.mat', 'X_star_coarse', 'U_star_coarse')

%% Interpolate
load('reference_traj_coarse.mat')
t_int = 0:0.02:t(end);
X_star = interp1(t, X_star_coarse', t_int, 'spline')';
U_star = interp1(t, U_star_coarse, t_int, 'spline');
save('reference_traj.mat', 'X_star', 'U_star')

%% Plots
close('all');
load('reference_traj.mat')
dip_plots(t_int, X_star, params, length(t_int), 1, 0.01);

%% Mesh refinement
N = 301;
dt = 0.02;
t = 0:dt:(N-1)*dt;

X0_ic = X_star;
U0_ic = U_star;
z0 = [X0_ic(:); U0_ic(:)];

% Inequality Constraints
lb = -inf(size(z0));
ub = inf(size(z0));
lb_x = reshape(lb(1:n*N), n, N); % n x N
ub_x = reshape(ub(1:n*N), n, N);
lb_u = reshape(lb(n*N+1:end), 1, N); % m x N
ub_u = reshape(ub(n*N+1:end), 1, N);

for k = 1:N
    lb_x(1,k) = -8;
    ub_x(1,k) = 8;

    lb_u(1,k) = -50;
    ub_u(1,k) = 50;
end
lb = [lb_x(:); lb_u(:)];
ub = [ub_x(:); ub_u(:)];

Z_opt = fmincon(@(z)cost_traj(z, n, N), z0, [], [], [], [], lb, ub, @(z)nonlcon(z,n,N,dt,params,x0,xstar), options);

X_star_coarse = reshape(Z_opt(1:n*N), n, N);
U_star_coarse = Z_opt(n*N+1:end)';
save('reference_traj_coarse.mat', 'X_star_coarse', 'U_star_coarse')

%% Interpolate 2
load('reference_traj_coarse.mat')
t_int = 0:0.01:t(end);
X_star = interp1(t, X_star_coarse', t_int, 'spline')';
U_star = interp1(t, U_star_coarse, t_int, 'spline');
save('reference_traj.mat', 'X_star', 'U_star')

%% Plots 2
close('all');
load('reference_traj.mat')
dip_plots(t_int, X_star, params, length(t_int), 1, 0.01);