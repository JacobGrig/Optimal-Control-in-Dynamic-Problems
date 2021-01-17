function f_dot = odefun_first(t, x, k_1, k_2, alpha, psi_1)
x_1 = x(1);
x_2 = x(2);
psi_2 = x(3);
u_1 = (psi_2 - alpha > 0) .* (psi_2 - alpha) ./ 2;
u_2 = (psi_1 + psi_2 .* x_2 > 0) .* k_2 + (psi_1 + psi_2 .* x_2 < 0) .* k_1;
f_dot = ...
    [...
    x_2 + u_2; ...
    u_1 + u_2 .* x_2;
    -psi_1 - psi_2 .* u_2 + 0 .* (t + x_1)
    ];
end