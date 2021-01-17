function draw_u2_t(u_opt, t_psi_opt, u2_curs)

u2_t_plot = figure('name', 'Time Dependence of u_2(t)');
title('Time Dependence of u_2(t)');
hold on;
xlabel('t');
ylabel('u_1');

grid;

if (~isempty(u_opt))
    [plot_u2_t, plot_u2_t_opt] = draw_time(u2_t_plot, u_opt(2, :).', t_psi_opt, u2_curs);

    legend([plot_u2_t(1), plot_u2_t_opt(1)], ...
        'u_2(t)', 'optimal u_2(t)');
end
end