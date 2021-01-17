function draw_u2_t(u_opt, switches)

u2_t_plot = figure('name', 'Time Dependence of u_2(t)');
title('Time Dependence of u_2(t)');
hold on;
xlabel('t');
ylabel('u_2');

grid;

[plot_u2_t, plot_u2_t_switches] = draw_time(u2_t_plot, u_opt, 2, switches);

if (~isempty(switches))

    legend([plot_u2_t(1), plot_u2_t_switches(1)], ...
        'u_2(t)', 'switches');
end
end