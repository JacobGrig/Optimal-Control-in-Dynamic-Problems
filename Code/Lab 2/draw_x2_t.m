function draw_x2_t(x_opt, switches)

x2_t_plot = figure('name', 'Time Dependence of x_2(t)');
title('Time Dependence of x_2(t)');
hold on;
xlabel('t');
ylabel('x_2');

grid;

[plot_x2_t, plot_x2_t_switches] = draw_time(x2_t_plot, x_opt, 2, switches);

if (~isempty(switches))

    legend([plot_x2_t(1), plot_x2_t_switches(1)], ...
        'x_2(t)', 'switches');
end
end