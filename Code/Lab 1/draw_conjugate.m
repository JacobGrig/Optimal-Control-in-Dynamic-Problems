function draw_conjugate(psi_opt, psi1_curs, psi2_curs)

psi1_psi2_plot = figure('name', 'Conjugate Function');
title('Conjugate Function');
hold on;
xlabel('\psi_1');
ylabel('\psi_2');

grid;

[plot_psi, plot_psi_opt] = draw_phase(psi1_psi2_plot, psi_opt, psi1_curs, psi2_curs, false);

if (~isempty(psi_opt))
    legend([plot_psi(1), plot_psi_opt(1)], ...
            'conjugate functions', 'optimal conjugate function');
end
end