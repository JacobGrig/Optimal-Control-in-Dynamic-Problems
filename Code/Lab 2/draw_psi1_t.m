function draw_psi1_t(psi_opt, switches)

psi1_t_plot = figure('name', 'Time Dependence of psi_1(t)');
title('Time Dependence of \psi_1(t)');
hold on;
xlabel('t');
ylabel('\psi_1');

grid;

[plot_psi1_t, plot_psi1_t_switches] = draw_time(psi1_t_plot, psi_opt, 1, switches); 

if (~isempty(switches))
    
    legend([plot_psi1_t(1), plot_psi1_t_switches(1)], ...
        '\psi_1(t)', 'switches');
end
end