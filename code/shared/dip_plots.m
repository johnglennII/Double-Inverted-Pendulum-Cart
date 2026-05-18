function [f1, varargout] = dip_plots(t, x_history, params, num_samples, refresh_rate, u_history, ctrl_mode_history)
%Double Inverted Pendulum Animation
% INPUTS:
%   x_history: n x N

arguments
    t
    x_history
    params
    num_samples
    refresh_rate
    u_history = nan(1,length(t))
    ctrl_mode_history = ones(1, num_samples)
end

% num_samples = length(t);

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
f1 = figure('Name', 'Double Inverted Pendulum', 'Color', 'w', 'WindowState','maximized');
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
h_ctrl_mode = text(2, (L1+L2+.35), 'Tracking', 'Color', '#D95319');
time = text(2, (L1+L2+.45), 't=0.0 s');

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
        set(h_ctrl_mode, 'String', 'Tracking', 'Color', '#D95319');
    elseif ctrl_mode_history(idx) == 2
        set(h_ctrl_mode, 'String', 'Regulation', 'Color', '#77AC30');
    end

    % update positions
    set(h_cart, 'Position', [j1_x(idx)-.5*cart_width, -.5*cart_height, cart_width, cart_height])
    set(h_link1, 'XData', [j1_x(idx), j2_x(idx)], 'YData', [j1_y(idx), j2_y(idx)])
    set(h_link2, 'XData', [j2_x(idx), ee_x(idx)], 'YData', [j2_y(idx), ee_y(idx)])
    set(h_joints, 'XData', [j1_x(idx), j2_x(idx)], 'Ydata', [j1_y(idx), j2_y(idx)])
    set(time, 'string', sprintf('t=%.1f s', t(idx)));
    drawnow
end
t_anim = toc(startTime);
fprintf('Animation time: %.2fs\n', t_anim);


if nargout > 1
    % --xc, theta, u vs time--
    varargout{1} = figure;
    subplot(3,1,1)
    hold on; grid on;
    title('$x_c$')
    plot(t, x1_history)

    subplot(3,1,2)
    hold on; grid on;
    title('$\theta_1$ and $\theta_2$')
    plot(t, [x3_history; x5_history])
    legend('$\theta_1$', '$\theta_2$')
    
    subplot(3,1,3)
    hold on; grid on;
    title('Control Input: u')
    plot(t, u_history);
    xlabel('time (s)')
    
    % states vs time
    varargout{2} = figure;
    subplot(3,1,1)
    hold on; grid on;
    title('Cart')
    plot(t, [x1_history; x2_history])
    legend('$x_c$', '$\dot{x}_c$')
    
    subplot(3,1,2)
    hold on; grid on;
    title('Pendulum 1')
    plot(t, [x3_history; x4_history])
    legend('$\theta_1$', '$\dot{\theta}_1$')
    
    subplot(3,1,3)
    hold on; grid on;
    title('Pendulum 2')
    plot(t, [x5_history; x6_history])
    legend('$\theta_2$', '$\dot{\theta}_2$')
    xlabel('time (s)')
end


end