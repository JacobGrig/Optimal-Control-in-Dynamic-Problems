function draw_control(u_opt, u1_curs, u2_curs, a, b, c, number_of_points_for_splitting_P)

u1_u2_plot = figure('name', 'Optimal Control');
title('Optimal Control');
hold on;
xlabel('u_1');
ylabel('u_2');

grid;

u1_splitting = linspace(-sqrt(b ./ (a + c)), sqrt(b ./ (a + c)), ...
    number_of_points_for_splitting_P);
plot(u1_splitting, a .* (u1_splitting .^ 2), 'k');
plot(u1_splitting, b - c .* (u1_splitting .^ 2), 'k');

[plot_u, plot_u_opt] = draw_phase(u1_u2_plot, u_opt.', u1_curs, u2_curs, true);
if (~isempty(u_opt))
    legend([plot_u(1), plot_u_opt(1)], ...
        'sub-optimal controls', 'optimal control');
end
end