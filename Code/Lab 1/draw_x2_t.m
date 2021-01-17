function draw_x2_t(x_opt, t_traj_opt, x2_curs)
x2_t_plot = figure('name', 'Time Depencence of x_2(t)');
title('Time Depencence of x_2(t)');
hold on;
xlabel('t');
ylabel('x_2');

grid;

if (~isempty(x_opt))
    [plot_x2_t, plot_x2_t_opt] = draw_time(x2_t_plot, x_opt(:, 2), t_traj_opt, x2_curs);
    legend([plot_x2_t(1), plot_x2_t_opt(1)], ...
        'x_2(t)', 'optimal x_2(t)');
end
end