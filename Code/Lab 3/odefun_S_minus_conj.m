function dfunc = odefun_S_minus_conj(t, x, alpha)
dfunc = [x(2) + 0 .* t; x(2) - 2 * x(1) - 3 * x(1) * sin(x(1) ^ 2) + 2 * x(1) ^ 2 * cos(x(1)) - alpha; ...
    x(4) * (3 * sin(x(1) ^ 2) + 2 * x(1) ^ 2 * sin(x(1)) + 6 * x(1) ^ 2 * cos(x(1) ^ 2) - 4 * x(1) * cos(x(1)) + 2); -x(3) - x(4)];
end