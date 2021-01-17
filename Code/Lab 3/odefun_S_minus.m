function dfunc = odefun_S_minus(t, x, alpha)
dfunc = [x(2) + 0 .* t; x(2) - 2 * x(1) - 3 * x(1) * sin(x(1) ^ 2) + 2 * x(1) ^ 2 * cos(x(1)) - alpha];
end