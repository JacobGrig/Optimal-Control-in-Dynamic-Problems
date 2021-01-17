function draw_psi2_t(psi_opt, t_psi_opt, psi2_curs)

psi2_t_plot = figure('name', 'Time Dependence of psi_2(t)');
title('Time Dependence of \psi_2(t)');
hold on;
xlabel('t');
ylabel('\psi_2');

grid;

if (~isempty(psi_opt))
    [plot_psi2_t, plot_psi2_t_opt] = draw_time(psi2_t_plot, psi_opt(:, 2), t_psi_opt, psi2_curs); 
    legend([plot_psi2_t(1), plot_psi2_t_opt(1)], ...
        '\psi_2(t)', 'optimal \psi_2(t)');
end
end