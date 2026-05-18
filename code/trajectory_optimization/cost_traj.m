function [J_z, grad_J] = cost_traj(z, n, N)
%NLP Cost Function
% INPUTS:
%   n: num states
%   N: num time steps
%   z: [x1, x2, ... xN, u1, u2, ..., uN]' (n+1)*N x 1

u = z(n*N+1:end);

J_z = norm(u, 2)^2;

if nargout > 1
    grad_J = zeros(size(z));
    grad_J(n*N+1:end) = 2*u;

end