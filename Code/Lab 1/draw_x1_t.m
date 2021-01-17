function draw_x1_t(x_opt, t_traj_opt, x1_curs)
x1_t_plot = figure('name', 'Time Depencence of x_1(t)');
title('Time Depencence of x_1(t)');
hold on;
xlabel('t');
ylabel('x_1');

grid;

if (~isempty(x_opt))
    [plot_x1_t, plot_x1_t_opt] = draw_time(x1_t_plot, x_opt(:, 1), t_traj_opt, x1_curs);
    legend([plot_x1_t(1), plot_x1_t_opt(1)], ...
        'x_1(t)', 'optimal x_1(t)');
end
end