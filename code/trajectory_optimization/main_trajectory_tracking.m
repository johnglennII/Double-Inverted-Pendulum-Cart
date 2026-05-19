%% Robotics Final Project
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

x0 = [0.0; 0.0; pi; 0.0; pi; 0]; % down-down
% x0 = zeros(6,1); % up-up

x_star = zeros(6,1); % up-up
% xstar = [0; 0; 0; 0; pi; 0]; % up-down
% xstar = [0; 0; pi; 0; 0; 0]; % down-up
u_star = 0;
n = length(x0);

% integrator settings
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
tf = 10;
dt = 0.01;
t = 0:dt:tf;
num_samples = length(t);

load('reference_traj.mat')
N_su = size(X_star,2);
X_star = [X_star, zeros(n, length(t)-N_su)];
U_star = [U_star, zeros(1, length(t)-N_su)];

% === Control Design ===
R = .001;
Q = diag([.4, .01, 15, 5, 15, .5]);

N_mpc = 30;
dt_mpc = 0.01;
[A,B,C] = compute_linearized_dynamics(x_star, u_star, params);
sys_c = ss(A,B,C,0);
sys_d = c2d(sys_c, dt_mpc);
Af = sys_d.A; Bf = sys_d.B;
[~, Qf, ~] = dlqr(Af, Bf, Q, R);

% pre-allocate
x_history = nan(n, num_samples);
u_history = nan(1, num_samples);
x_current = x0;
x_history(:,1) = x_current;
ctrl_mode_history = ones(1, num_samples);

tracking = true;
tic
for k = 1:num_samples-1
    t_current = t(k);
    t_next = t(k+1);
    
    % Tracking Controller
    if tracking
        x_star_k = X_star(:,k);
        u_star_k = U_star(k);
        [A_k,B_k,C_k] = compute_linearized_dynamics(x_star_k, u_star_k, params);
        sys_c_k = ss(A_k,B_k,C_k,0);
        sys_dk = c2d(sys_c_k, dt_mpc);
        Ad_k = sys_dk.A; Bd_k = sys_dk.B;
        [H, F] = setup_LMPC(Ad_k,Bd_k,Q,R,Q,N_mpc);
    elseif ~tracking
        % [A_k,B_k,C_k] = compute_linearized_dynamics(x_star, u_star, params);
        % sys_c_k = ss(A_k,B_k,C_k,0);
        % sys_dk = c2d(sys_c_k, dt_mpc);
        % Ad_k = sys_dk.A; Bd_k = sys_dk.B;
        [H, F] = setup_LMPC(Af,Bf,Q,R,Qf,N_mpc);
        ctrl_mode_history(k+1) = 2;
    end

    du_current = solve_LMPC(x_current-x_star_k, H, F, -50*ones(N_mpc,1), 50*ones(N_mpc,1));
    u_current = u_star_k + du_current;

    [~, x_out] = ode45(@(t,x) dynamics_double_inv_pend(t,x,u_current,params), [t_current, t_next], x_current, options);
    x_current = x_out(end,:)';

    % store data
    x_history(:,k+1) = x_current;
    u_history(:,k+1) = u_current;

    if abs(x_current(1)) > 20
        warning('Cart pos exceeded limits at t=%.2fs', t(k))
        num_samples = k;
        break
    end

    x_wrapped = x_current; x_wrapped(3) = wrapToPi(x_current(3)); x_wrapped(5) = wrapToPi(x_current(5));
    if abs(x_wrapped(3)) < 1*pi/180 && abs(x_wrapped(5)) < 1*pi/180
        tracking = false;
    end
end
t_sim = toc;
fprintf('True time: %.2fs\n', tf);
fprintf('Simulation time: %.2fs\n', t_sim);

%% ---Plots---
close('all');
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 'remove')

export_plt = input('Export Plots? y/n: ', 's');

refresh_rate = 1;

[f1, f2, f3] = dip_plots(t, x_history, params, num_samples, refresh_rate, u_history, ctrl_mode_history);
t_switch = t(find(ctrl_mode_history==2, 1, 'first'));

% Reference vs Actual: States
f4 = figure('WindowState', 'maximized');
tl = tiledlayout(6,1, 'Padding', 'compact');
title(tl, sprintf('State Trajectory: Actual vs Reference'), 'interpreter', 'latex');
for l = 1:n
    ax(l) = nexttile;
    hold on; grid on;
    h = plot(t, [x_history(l,:); X_star(l,:)]);
    ylabel(sprintf('$x_{%i}$', l))
end
xlabel('time (s)')
lgd = legend(ax(1), h(1:2), {'$x_{actual}$', '$X^*$'});
lgd.Location = 'northeast';

% Tracking Error: States
tracking_err_x = x_history - X_star;
f5 = figure('WindowState', 'maximized');
tl = tiledlayout(6,1, 'Padding', 'compact');
title(tl, sprintf('Tracking Error: States'), 'interpreter', 'latex');
for l = 1:n
    ax(l) = nexttile;
    hold on; grid on;
    h = plot(t, tracking_err_x(l,:));
    yline(0, 'k--')
    xline(t_switch, 'g--')
    ylabel(sprintf('$x_{%i}$ error', l))
end
lgd = legend(ax(1), {'', '','switch time'});
lgd.Location = 'northeast';
xlabel('time (s)')

% Reference vs Actual: Inputs
tracking_err_u = u_history - U_star;
f6 = figure('WindowState', 'maximized');
% both
subplot(2,1,1)
hold on; grid on;
plot(t, [u_history; U_star]);
title('Input Trajectory: Actual vs Reference');
legend('$u_{actual}$', '$U^*$')

% u error
subplot(2,1,2)
hold on; grid on;
plot(t, tracking_err_u);
yline(0, 'k--')
xline(t_switch, 'g--')
title('Tracking Error: Input');
legend('', '','switch time');
xlabel('time (s)')

if export_plt == 'y'
    currentPath = fileparts(mfilename('fullpath'));
    mediaPath = fullfile(currentPath, '..', '..', 'media', 'trajectory_optimization');
    exportgraphics(f2, fullfile(mediaPath, 'pos_u.png'), 'Resolution', 300);
    exportgraphics(f3, fullfile(mediaPath, 'states.png'), 'Resolution', 300);
    exportgraphics(f4, fullfile(mediaPath, 'state_traj.png'), 'Resolution', 300);
    exportgraphics(f5, fullfile(mediaPath, 'state_error.png'), 'Resolution', 300);
    exportgraphics(f6, fullfile(mediaPath, 'input_traj.png'), 'Resolution', 300);
end