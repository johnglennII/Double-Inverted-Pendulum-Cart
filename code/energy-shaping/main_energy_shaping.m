%% Robotics Final Project
clc;clear;close('all');

% Define constants
mc = 400e-3;
m1 = 0.15;
m2 = 0.15;
L1 = 0.5;
L2 = 0.5;
l1 = L1/2;
l2 = L2/2;
J1 = m1*L1^2/12;
J2 = m2*L2^2/12;
cc = 0.0;
c1 = 0.0;
c2 = 0.0;
g = 9.81;

params.mc = mc;
params.m1 = m1;
params.m2 = m2;
params.L1 = L1;
params.L2 = L2;
params.l1 = l1;
params.l2 = l2;
params.J1 = J1;
params.J2 = J2;
params.cc = cc;
params.c1 = c1;
params.c2 = c2;
params.g = g;

% x0 = [0; 0; 10*pi/180; 0; -10*pi/180; 0];
x0 = [0.0; 0.01; pi; 0.0; pi; 0]; % down-down
% x0 = zeros(6,1); % up-up

xstar = zeros(6,1); % up-up
% xstar = [0; 0; 0; 0; pi; 0]; % up-down
% xstar = [0; 0; pi; 0; 0; 0]; % down-up
ustar = 0;
n = length(x0);

% control design
[A,B,C] = compute_linearized_dynamics(xstar, ustar, params);
Q = diag([.1, .01, 15, 5, 15, .5]); % WORKING
R = 1;

% WORKING GAINS:
K_gain = 5;
k_p = 5;
k_d = 3;

E_desired = compute_energy(xstar, params);

[K_lqr,~,poles] = lqr(A,B,Q,R);

dt_mpc = 0.01;
sys_c = ss(A,B,C,0);
sys_d = c2d(sys_c, dt_mpc);
Ad = sys_d.A; Bd = sys_d.B;
[~, Qf, ~] = dlqr(Ad, Bd, Q, R);
N = 30;
[H, F] = setup_LMPC(Ad,Bd,Q,R,Qf,N);


% integrator settings
options = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
tf = 30;
dt = 0.001;
t = 0:dt:tf;
num_samples = length(t);

% pre-allocate
x_history = nan(n, num_samples);
u_history = nan(1, num_samples);
v_history = nan(1, num_samples);
E_history = nan(1, num_samples);
ctrl_mode_history = ones(1, num_samples);

x_current = x0;
x_history(:,1) = x_current;

switched = false;
tic
for k = 1:num_samples-1
    t_current = t(k);
    t_next = t(k+1);

    x_wrapped = x_current; x_wrapped(3) = wrapToPi(x_current(3)); x_wrapped(5) = wrapToPi(x_current(5));
    
    if abs(x_wrapped(3)) < 15*pi/180 && abs(x_wrapped(5)) < 15*pi/180 && abs(x_current(6)) < 2 || switched
        % LMPC Controller
        u_current = solve_LMPC(x_wrapped, H, F, -75*ones(N,1), 75*ones(N,1));
        ctrl_mode_history(k+1) = 2;
        switched = true;
        
        % u_current = -K_lqr*(x_wrapped-xstar);

    elseif ~switched
        % Swing Up
        x_fake = x_current;
        x_fake(2) = 0;
        E_pend = compute_energy(x_fake, params);
        W = (m1*l1 + m2*L1)*x_current(4)*cos(x_current(3)) + m2*l2*x_current(6)*cos(x_current(5));

        if abs(W) < 0.05 || true
            v = -K_gain*(E_pend - E_desired)*W - k_p*x_current(1) - k_d*x_current(2);
        elseif abs(W) > 0.05 % damping compensation, iffy numerically
            v = -K_gain*(E_pend - E_desired)*W - k_p*x_current(1) - k_d*x_current(2) + ((c1+c2)*x_current(4)^2 + c2*x_current(6)^2 - 2*c2*x_current(4)*x_current(6))/W;
        end
        v = max(min(v, 75), -75);

        % v = sin(5*t_current);
        [alpha, beta] = compute_alpha_beta(x_current, params);
        u_current = alpha + beta*v;
    end

    [~, x_out] = ode45(@(t,x) dynamics_double_inv_pend(t,x,u_current,params), [t_current, t_next], x_current, options);
    x_current = x_out(end,:)';

    % store data
    x_history(:,k+1) = x_wrapped;
    u_history(:,k+1) = u_current;
    E_history(:,k+1) = E_pend;
    v_history(:,k+1) = v;

    if abs(x_current(1)) > 20
        warning('Cart pos exceeded limits at t=%.2fs', t(k))
        num_samples = k;
        break
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

refresh_rate = 4;

x1_history = x_history(1,:);
x2_history = x_history(2,:);
x3_history = x_history(3,:);
x4_history = x_history(4,:);
x5_history = x_history(5,:);
x6_history = x_history(6,:);

% joint coordinates (x,y) using derived eqns
j1_x = x1_history;
j1_y = zeros(1, num_samples);

j2_x = j1_x - L1.*sin(x3_history);
j2_y = L1.*cos(x3_history);

ee_x = j2_x - L2.*sin(x5_history);
ee_y = j2_y + L2.*cos(x5_history);

% --Animation--
f1 = figure('Name', 'Double Inverted Pendulum', 'Color', 'w');
hold on; grid on; axis equal;

% axis limits
xlim([-2.5, 2.5])
ylim([-(L1+L2+.25), (L1+L2+.25)])
xlabel('Position (m)'); ylabel('Height (m)')
title('Double Inverted Pendulum Cart')

% plot first frame
cart_width = 0.2;
cart_height = 0.1;

h_cart = rectangle('Position', [j1_x(1)-.5*cart_width, -.5*cart_height, cart_width, cart_height]);
h_link1 = plot([j1_x(1), j2_x(1)], [j1_y(1), j2_y(1)], 'color', [0, 0.4470, 0.7410]);
h_link2 = plot([j2_x(1), ee_x(1)], [j2_y(1), ee_y(1)], 'color', [0, 0.4470, 0.7410]);
h_joints = plot([j1_x(1), j2_x(1)], [j1_y(1), j2_y(1)], 'o', 'MarkerSize', 3, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'Color', [0, 0.4470, 0.7410]);
h_ctrl_mode = text(2, (L1+L2+.45), 'Swing-up', 'Color', '#D95319');
time = text(2, 2, 't=0.0 s');
% time = annotation('textbox', [0.75 0.85 0.2 0.1], 'String', 't=0.0 s', 'EdgeColor', 'none');

pause(.75)
idx = 1;
startTime = tic;
% animation loop
while idx < num_samples
    elapsed = toc(startTime);

    % Find the index in the simulation time vector 't' that matches real-world time
    % We use a simple search here; for very high-frequency data, use binary search
    while idx < num_samples && t(idx) < elapsed
        idx = idx + 1;
    end

    if ctrl_mode_history(idx) == 1
        set(h_ctrl_mode, 'String', 'Swing-up', 'Color', '#D95319');
    elseif ctrl_mode_history(idx) == 2
        set(h_ctrl_mode, 'String', 'MPC', 'Color', '#77AC30');
    end
    % update positions
    set(h_cart, 'Position', [j1_x(idx)-.5*cart_width, -.5*cart_height, cart_width, cart_height])
    set(h_link1, 'XData', [j1_x(idx), j2_x(idx)], 'YData', [j1_y(idx), j2_y(idx)])
    set(h_link2, 'XData', [j2_x(idx), ee_x(idx)], 'YData', [j2_y(idx), ee_y(idx)])
    set(h_joints, 'XData', [j1_x(idx), j2_x(idx)], 'Ydata', [j1_y(idx), j2_y(idx)])
    set(time, 'string', sprintf('t=%.1f s', t(idx)));
    drawnow
end
t_anim = toc;
fprintf('Animation time: %.2fs\n', t_anim);


% --theta, u, E vs time--
f2 = figure;
subplot(3,1,1)
hold on; grid on;
title('$\theta_1$ and $\theta_2$')
plot(t, [x3_history; x5_history])
yline(0, 'k--')
legend('$\theta_1$', '$\theta_2$')

subplot(3,1,2)
hold on; grid on;
title('Control Input: u')
plot(t, u_history);

subplot(3,1,3)
hold on; grid on;
title('Pendulum Energy: $E_{pend}$')
plot(t, E_history);
yline(E_desired, 'k--')

% --cart, v, u vs time--
f3 = figure;
subplot(3,1,1)
hold on; grid on;
title('Cart States: $x_1, x_2$')
h = plot(t, [x1_history; x2_history]);
legend('$x_c$', '$\dot{x}_c$')

subplot(3,1,2)
hold on; grid on;
title('Virtual Input: v')
plot(t, v_history)
legend('v')

subplot(3,1,3)
hold on; grid on;
title('Control Input: u')
plot(t, u_history)
legend('u')

% states vs time
f4 = figure;
subplot(1,3,1)
hold on; grid on;
title('Cart')
plot(t, [x1_history; x2_history])
legend('$x_c$', '$\dot{x}_c$')

subplot(1,3,2)
hold on; grid on;
title('Pendulum 1')
plot(t, [x3_history; x4_history])
legend('$\theta_1$', '$\dot{\theta}_1$')

subplot(1,3,3)
hold on; grid on;
title('Pendulum 2')
plot(t, [x5_history; x6_history])
legend('$\theta_2$', '$\dot{\theta}_2$')

if export_plt == 'y'
    exportgraphics(f2, fullfile('Report', 'E_shaping', 'theta_u_E.png'), 'Resolution', 300);
    exportgraphics(f3, fullfile('Report', 'E_shaping', 'cart_v_u.png'), 'Resolution', 300);
    exportgraphics(f4, fullfile('Report', 'E_shaping', 'states.png'), 'Resolution', 300);
end