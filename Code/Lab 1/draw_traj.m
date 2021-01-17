function draw_traj(x_opt, x1_curs, x2_curs, r, d, k, ...
    trimmer, psi_cur_opt, psi_trans_opt, ...
    number_of_points_for_splitting_X0, number_of_points_for_splitting_X1)

x1_x2_plot = figure('name', 'Optimal and Trial Trajectories');
title('Optimal and Trial Trajectories');
hold on;
xlabel('x_1');
ylabel('x_2');

x1_max = max(r, sqrt(k)) + trimmer;
x2_max = max(r, abs(d) + k) + trimmer;

x_max = max(x1_max, x2_max);

axis([-x_max, x_max, -x_max, x_max]);

grid;

alpha = linspace(0, 2 * pi, number_of_points_for_splitting_X0);
plot(r .* cos(alpha), r .* sin(alpha), 'k');

x1_splitting = linspace(-sqrt(k), sqrt(k), number_of_points_for_splitting_X1);
plot(x1_splitting, d + k - x1_splitting .^ 2, 'k');
plot(x1_splitting, d - k + x1_splitting .^ 2, 'k');

[plot_x, plot_x_opt] = draw_phase(x1_x2_plot, x_opt, x1_curs, x2_curs, false);

if (~isempty(x_opt))
    final_point = x_opt(end, :);
    psi_cur_opt = psi_cur_opt ./ sqrt(psi_cur_opt(1) .^ 2 + psi_cur_opt(2) .^ 2);
    psi_final = final_point + psi_cur_opt.';
    psi_trans_final = final_point + psi_trans_opt.';

    plot_psi_cur = plot([final_point(1), psi_final(1)], ...
                        [final_point(2), psi_final(2)], 'g');
    plot_psi_trans = plot([final_point(1), psi_trans_final(1)], ...
                          [final_point(2), psi_trans_final(2)], 'm');
    legend([plot_x(1), plot_x_opt(1), plot_psi_cur(1), plot_psi_trans(1)], ...
        'trial trajectories', 'optimal trajectory', '-\psi(t_1)', '\psi^1');
end
end