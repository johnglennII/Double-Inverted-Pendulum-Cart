function u_current = solve_LMPC(x_k, H, F, lb, ub)
%Solve LMPC

f = F*x_k;

options = optimoptions('quadprog', 'Display', 'off');

U_opt = quadprog(2*H, 2*f, [], [], [], [], lb, ub, [], options);
u_current = U_opt(1);

end