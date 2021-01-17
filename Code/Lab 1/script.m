% Laboratory Work #1

%% System #1

clear;
clc;

A = @(t) [0 .* t, -3 + 0 .* t; 2 + 0 .* t, 2 + 0 .* t];
B = @(t) [1 + 0 .* t, 2 + 0 .* t; 0 .* t, 4 + 0 .* t];
a = 1;
b = 1;
c = 2;
k = 0.5;
r = 0.5;
d = 3;
number = 61;
number_X0 = 60;
number_X1 = 60;
number_P = 60;
t1_max = 1;
eps_X1 = 1e-3;
eps_B = 1e-3;
eps_rel = 1e-8;
eps_abs = 1e-10;

trimmer = 1;

[t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt, is_solvable] = ...
    improvement_of_accuracy(3, A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
    a, b, c, r, d, k, number, number_X0, number_X1, trimmer);

%%

[t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt] = ...
    main_algorithm(A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
    a, b, c, r, d, k, [0, 2*pi], number);

trimmer = 1;

draw_traj(x_opt, x1_curs, x2_curs, r, d, k, trimmer, ...
    psi_cur_opt, psi_trans_opt, number_X0, number_X1);
draw_x1_t(x_opt, t_traj_opt, x1_curs);
draw_x2_t(x_opt, t_traj_opt, x2_curs);

draw_control(u_opt, u1_curs, u2_curs, a, b, c, number_P);
draw_u1_t(u_opt, t_psi_opt, u1_curs);
draw_u2_t(u_opt, t_psi_opt, u2_curs);

draw_conjugate(psi_opt, psi1_curs, psi2_curs);
draw_psi1_t(psi_opt, t_psi_opt, psi1_curs);
draw_psi2_t(psi_opt, t_psi_opt, psi2_curs);

%% System #2

clear;
clc;

A = @(t) [1 + 0 .* t, 4 + 0 .* t; -1 + 0 .* t, -3 + 0 .* t];
B = @(t) [-1 + 0 .* t, 0 .* t; 2 + 0 .* t, -2 + 0 .* t];
a = 0.3;
b = 1;
c = 0.3;
k = 1;
r = 0.5;
d = -2;
number = 61;
number_X0 = 60;
number_X1 = 60;
number_P = 60;
t1_max = 1.2;
eps_X1 = 1e-3;
eps_B = 1e-3;
eps_rel = 1e-8;
eps_abs = 1e-10;

trimmer = 1;

[t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt, is_solvable] = ...
    improvement_of_accuracy(3, A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
    a, b, c, r, d, k, number, number_X0, number_X1, trimmer);
%%
trimmer = 1;

draw_traj(x_opt, x1_curs, x2_curs, r, d, k, trimmer, ...
    psi_cur_opt, psi_trans_opt, number_X0, number_X1);

%%

[t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt] = ...
    main_algorithm(A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
    a, b, c, r, d, k, [0, 2*pi], number);

trimmer = 1;

draw_traj(x_opt, x1_curs, x2_curs, r, d, k, trimmer, ...
    psi_cur_opt, psi_trans_opt, number_X0, number_X1);
draw_x1_t(x_opt, t_traj_opt, x1_curs);
draw_x2_t(x_opt, t_traj_opt, x2_curs);

draw_control(u_opt, u1_curs, u2_curs, a, b, c, number_P);
draw_u1_t(u_opt, t_psi_opt, u1_curs);
draw_u2_t(u_opt, t_psi_opt, u2_curs);

draw_conjugate(psi_opt, psi1_curs, psi2_curs);
draw_psi1_t(psi_opt, t_psi_opt, psi1_curs);
draw_psi2_t(psi_opt, t_psi_opt, psi2_curs);

%% System #3

clear;
clc;

A = @(t) [cos(t .^ 2), 0 .* t; -1 + 0 .* t, -2 .* cos(4 .* t .^ 2)];
B = @(t) [2 .* sin(t), -1 ./ (t + 1); 0 .* t, -sin(3 * t)];
a = 1;
b = 5;
c = 2;
k = 0.5;
r = 0.5;
d = 2;
number = 30;
number_X0 = 60;
number_X1 = 60;
number_P = 60;
t1_max = 5;
eps_X1 = 1e-3;
eps_B = 1e-3;
eps_rel = 1e-8;
eps_abs = 1e-10;
[t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt] = ...
    main_algorithm(A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
    a, b, c, r, d, k, [0, 2*pi], number);

trimmer = 1;

draw_traj(x_opt, x1_curs, x2_curs, r, d, k, trimmer, ...
    psi_cur_opt, psi_trans_opt, number_X0, number_X1);
draw_x1_t(x_opt, t_traj_opt, x1_curs);
draw_x2_t(x_opt, t_traj_opt, x2_curs);

draw_control(u_opt, u1_curs, u2_curs, a, b, c, number_P);
draw_u1_t(u_opt, t_psi_opt, u1_curs);
draw_u2_t(u_opt, t_psi_opt, u2_curs);

draw_conjugate(psi_opt, psi1_curs, psi2_curs);
draw_psi1_t(psi_opt, t_psi_opt, psi1_curs);
draw_psi2_t(psi_opt, t_psi_opt, psi2_curs);

%% System #4

clear;
clc;

A = @(t) [0.2 .* sin(t), 0 .* t; 0 .* t, 0.2 .* cos(t)];
B = @(t) [0 .* t, 0.2 + 0 .* t; -0.3 + 0 .* t, 0 .* t];
a = 1;
b = 1;
c = 2;
k = 0.87;
r = 0.5;
d = -3;
number = 61;
number_X0 = 60;
number_X1 = 60;
number_P = 60;
t1_max = 20;
eps_X1 = 1e-3;
eps_B = 1e-3;
eps_rel = 1e-8;
eps_abs = 1e-10;

trimmer = 1;
% [t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt] = ...
%     main_algorithm(A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
%     a, b, c, r, d, k, [0, 2*pi], number, number_X0, number_X1, number_P);
% [t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
%     x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
%     t_psi_opt, psi_cur_opt, psi_trans_opt] = ...
%     main_algorithm(A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
%     a, b, c, r, d, k, [0, 2*pi], number);
[t1_min, u_opt, x_opt, phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, psi1_curs, psi2_curs, ...
    t_psi_opt, psi_cur_opt, psi_trans_opt, is_solvable] = ...
    improvement_of_accuracy(3, A, B, t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...
    a, b, c, r, d, k, number, number_X0, number_X1, trimmer);

draw_traj(x_opt, x1_curs, x2_curs, r, d, k, trimmer, ...
    psi_cur_opt, psi_trans_opt, number_X0, number_X1);
