function res = ode_fund(t, fund, A_matr)
fund_matr = reshape(fund, 2, 2);
res_matr = A_matr(t) * fund_matr;
res = res_matr(:);
end