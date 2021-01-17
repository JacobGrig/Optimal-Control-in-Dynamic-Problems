function draw_psi1_t(psi_opt, t_psi_opt, psi1_curs)

psi1_t_plot = figure('name', 'Time Dependence of psi_1(t)');
title('Time Dependence of \psi_1(t)');
hold on;
xlabel('t');
ylabel('\psi_1');

grid;

if (~isempty(psi_opt))
    [plot_psi1_t, plot_psi1_t_opt] = draw_time(psi1_t_plot, psi_opt(:, 1), t_psi_opt, psi1_curs); 
    legend([plot_psi1_t(1), plot_psi1_t_opt(1)], ...
        '\psi_1(t)', 'optimal \psi_1(t)');
end
end