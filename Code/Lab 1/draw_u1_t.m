function draw_u1_t(u_opt, t_psi_opt, u1_curs)

u1_t_plot = figure('name', 'Time Dependence of u_1(t)');
title('Time Dependence of u_1(t)');
hold on;
xlabel('t');
ylabel('u_1');

grid;

if (~isempty(u_opt))
    [plot_u1_t, plot_u1_t_opt] = draw_time(u1_t_plot, u_opt(1, :).', t_psi_opt, u1_curs);

    legend([plot_u1_t(1), plot_u1_t_opt(1)], ...
        'u_1(t)', 'optimal u_1(t)');
end
end