function func = fsolve_func_1_case_2(x_1, x_2, T, tau1, psi_1, psi_2_0, S, L, epsilon)
func(1) = L + epsilon - x_1(T, psi_1, psi_2_0, tau1);
func(2) = S - epsilon - x_2(T, psi_1, psi_2_0, tau1);
end