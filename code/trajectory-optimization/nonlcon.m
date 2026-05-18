function [g, h, dg, dh] = nonlcon(z, n, N, dt, params, x0, xf)
%Nonlinear Constraints
% INPUTS:
%   z: [x1, x2, ... xN, u1, u2, ..., uN]' (n+1)*N x 1
%   n: num states
%   N: num time steps

X = reshape(z(1:n*N), n, N); % n x N states
U = reshape(z(n*N+1:end), 1, N); % 1 x N inputs

h = zeros(n*(N-1), 1);
g = [];

num_con = n*(N-1) + 2*n;
num_vars = N*(n+1);

% Pre-allocate the sparse analytical Jacobian
    if nargout > 2
        dg = [];
        dh = spalloc(num_vars, num_con, num_con*14);
    end

idx = 1;
for k = 1:N-1
    x_k = X(:,k);
    x_k1 = X(:,k+1);
    u_k = U(k);
    u_k1 = U(k+1);

    f_k = dynamics_double_inv_pend(0, x_k, u_k, params);
    f_k1 = dynamics_double_inv_pend(0, x_k1, u_k1, params);

    % h_k = x_k1 - x_k - dt*f_k1; % backwards euler int
    h_k = x_k1 - x_k - .5*dt*(f_k + f_k1); % trapezoidal int
    h(idx:idx+n-1) = h_k;

    % Construct Analytical Jacobian if requested by fmincon
    if nargout > 2
        % Get continuous Jacobians from your existing MPC function!
        [A_k, B_k, ~] = compute_linearized_dynamics(x_k, u_k, params);
        [A_k1, B_k1, ~] = compute_linearized_dynamics(x_k1, u_k1, params);
        
        % Derivatives of h_k with respect to x_k, u_k, x_k1, u_k1
        dx_k  = -eye(n) - 0.5*dt*A_k;
        du_k  = -0.5*dt*B_k;
        dx_k1 = eye(n) - 0.5*dt*A_k1;
        du_k1 = -0.5*dt*B_k1;
        
        % Indices for variables in the z vector
        idx_x_k  = (k-1)*n + 1 : k*n;
        idx_x_k1 = k*n + 1 : (k+1)*n;
        idx_u_k  = n*N + k;
        idx_u_k1 = n*N + k + 1;
        
        % Indices for constraints in the h vector
        idx_h = idx : idx+n-1;
        
        % Assign to transposed Jacobian (rows = vars, cols = constraints)
        dh(idx_x_k, idx_h)  = dx_k';
        dh(idx_x_k1, idx_h) = dx_k1';
        dh(idx_u_k, idx_h)  = du_k';
        dh(idx_u_k1, idx_h) = du_k1';
    end

    idx = idx + n;
end

% IC & BC constraint
h = [h; X(:,1)-x0; X(:,end)-xf];

% Jacobians for IC & BC constraint
if nargout > 2
    % Initial state constraint Jacobians (I)
    idx_h = idx : idx+n-1;
    dh(1:n, idx_h) = eye(n);
    idx = idx + n;
    
    % Final state constraint Jacobians (I)
    idx_h = idx : idx+n-1;
    dh(n*(N-1)+1 : n*N, idx_h) = eye(n);
end

end