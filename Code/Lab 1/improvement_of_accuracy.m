function [t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt, is_solvable] ...
 = improvement_of_accuracy(number_of_improvements, ... % number of improvements
    A_matr, B_matr, ...                                % parameters of initial task
    t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...       % user's parameters
    a, b, c, r, d, k, ...                              % parameters of the sets P, X0, X1
    number_of_points_for_splitting, ...                % auxiliary parameters
    number_X0, number_X1, trimmer)

% initial interval
interval = [0, 2 * pi];

% the flag, if the task is solvable
is_solvable = false;

% repeat this number of times
for i = 1 : number_of_improvements
    [t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt] = ...
        main_algorithm(A_matr, B_matr, ...
        t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
        a, b, c, r, d, k, ...
        interval, number_of_points_for_splitting);
    draw_traj(x_opt, x1_curs, x2_curs, r, d, k, trimmer, ...
        psi_cur_opt, psi_trans_opt, number_X0, number_X1);
    if (~isempty(u_opt))
        % local improvement
        t1_max = t1_min;
        len = interval(2) - interval(1);
        interval = [phi_opt - len / 6, phi_opt + len / 6];
        is_solvable = true;
    else
        % global improvement
        number_of_points_for_splitting = number_of_points_for_splitting * 3;
    end
end
end