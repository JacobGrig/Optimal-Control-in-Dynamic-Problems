function draw_x1_t(x_opt, switches)

x1_t_plot = figure('name', 'Time Dependence of x_1(t)');
title('Time Dependence of x_1(t)');
hold on;
xlabel('t');
ylabel('x_1');

grid;

[plot_x1_t, plot_x1_t_switches] = draw_time(x1_t_plot, x_opt, 1, switches);

if (~isempty(switches))

    legend([plot_x1_t(1), plot_x1_t_switches(1)], ...
        'x_1(t)', 'switches');
end
end