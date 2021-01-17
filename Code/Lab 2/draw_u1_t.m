function draw_u1_t(u_opt, switches)

u1_t_plot = figure('name', 'Time Dependence of u_1(t)');
title('Time Dependence of u_1(t)');
hold on;
xlabel('t');
ylabel('u_1');

grid;

[plot_u1_t, plot_u1_t_switches] = draw_time(u1_t_plot, u_opt, 1, switches);

if (~isempty(switches))

    legend([plot_u1_t(1), plot_u1_t_switches(1)], ...
        'u_1(t)', 'switches');
end
end