function [H, F] = setup_LMPC(A, B, Q, R, Qf, N)
%

n = size(A,1); % num states
m = size(B,2); % num inputs

% lifted dynamics
A_lifted = zeros(n*(N+1), n);
B_lifted = zeros(n*(N+1), m*N);
A_lifted(1:n,:) = eye(n);
for i = 1:N
    row_idx = i*n+1:i*n+n;
    A_lifted(row_idx, :) = A^i;

    for j = 1:i
        col_idx = (j-1)*m+1:j*m;
        B_lifted(row_idx, col_idx) = (A^(i-j))*B;
    end
end

% Q and R
Q_cells = repmat({Q}, 1, N); % 1xN cells of nxn Q
R_cells = repmat({R}, 1, N); % 1xN cells of mxm R

Q_tilde = blkdiag(Q_cells{:}, Qf); % n*(N+1) square
R_tilde = blkdiag(R_cells{:}); % m*N square

% objective function matrices
H = B_lifted'*Q_tilde*B_lifted + R_tilde;
H = (H+H')/2;
F = B_lifted'*Q_tilde*A_lifted;


end