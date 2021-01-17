function func = fsolve_func_1_case_1(x_1, x_2, T, psi, S, L, epsilon)
func(1) = L + epsilon - x_1(T, psi(1), psi(2));
func(2) = S - epsilon - x_2(T, psi(1), psi(2));
end