function draw_traj(x_opt, switches, S, L, epsilon)

x1_x2_plot = figure('name', 'Optimal Trajectory');
title('Optimal Trajectory');
hold on;
xlabel('x_1');
ylabel('x_2');

x1_left = L - epsilon;
x1_right = L + epsilon;

x2_lower = S - epsilon;
x2_upper = S + epsilon;

plot([x1_left, x1_right], [x2_lower, x2_lower], 'k');
plot([x1_left, x1_right], [x2_upper, x2_upper], 'k');
plot([x1_left, x1_left], [x2_lower, x2_upper], 'k');
plot([x1_right, x1_right], [x2_lower, x2_upper], 'k');

start_plot = plot(0, 0, 'ok');

[plot_x, plot_x_switches] = draw_phase(x1_x2_plot, x_opt, switches);

if (~isempty(switches))
    legend([plot_x(1), plot_x_switches(1), start_plot], ...
        'optimal trajectory', 'switches', 'start point');
end
end