function draw_conjugate(psi_opt, switches)

psi1_psi2_plot = figure('name', 'Conjugate Function');
title('Conjugate Function');
hold on;
xlabel('\psi_1');
ylabel('\psi_2');

[plot_psi, plot_psi_switches] = draw_phase(psi1_psi2_plot, psi_opt, switches);

if (~isempty(switches))
    legend([plot_psi(1), plot_psi_switches(1)], ...
            'conjugate functions', 'switches');
end
end