function draw_psi2_t(psi_opt, switches)

psi2_t_plot = figure('name', 'Time Dependence of psi_2(t)');
title('Time Dependence of \psi_2(t)');
hold on;
xlabel('t');
ylabel('\psi_2');

grid;

[plot_psi2_t, plot_psi2_t_switches] = draw_time(psi2_t_plot, psi_opt, 2, switches); 

if (~isempty(switches))
    
    legend([plot_psi2_t(1), plot_psi2_t_switches(1)], ...
        '\psi_2(t)', 'switches');
end
end