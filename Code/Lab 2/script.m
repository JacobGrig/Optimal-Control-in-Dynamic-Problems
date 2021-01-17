% Laboratory Work #2

%% P = P_1

%% System #1

% optimal solution (J = 0)

k_1 = -15;
k_2 = 10;

alpha = 8;

T = 1;

S = 2;
L = -10;
epsilon = 4;

delta = 1e-1;

number_of_points_for_splitting = 10;
number_of_points_for_enum = 10;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = first_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #2

% without switches

k_1 = -2;
k_2 = 0.1;

alpha = 8;

T = 1;

S = 5;
L = 0;
epsilon = 4;

delta = 1e-1;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 100;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = first_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #3

% one switch

k_1 = -15;
k_2 = 1;

alpha = 8;

T = 1;

S = 5;
L = 1.58;
epsilon = 4;

delta = 1e-1;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 100;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = first_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

disp(['Optimal value of functional is ', num2str(functional_min), '!']);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #4

% two switches

k_1 = -1;
k_2 = 1;

alpha = 10;

T = 1;

S = 0.1849;
L = 6.898;
epsilon = 3;

delta = 1e-2;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 100;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = first_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

disp(['Optimal value of functional is ', num2str(functional_min), '!']);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% P = P_2

%% System #5

% optimal solution (J = -alpha ^ 2 / 4 * T)

k_1 = -15;
k_2 = 10;

alpha = 8;

T = 1;

S = 2;
L = -10;
epsilon = 4;

delta = 1e-1;

number_of_points_for_splitting = 10;
number_of_points_for_enum = 10;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #6

% without switches

k_1 = -2;
k_2 = 0.1;

alpha = 8;

T = 1;

S = 5;
L = 0;
epsilon = 4;

delta = 1e-1;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 20;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #7

% without switches (in difference with P_1, where we had one switch)

k_1 = -15;
k_2 = 1;

alpha = 8;

T = 1;

S = 5;
L = 1.58;
epsilon = 4;

delta = 1e-1;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 20;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #8

% without switches (in P_1 we had two switches!)

k_1 = -1;
k_2 = 1;

alpha = 10;

T = 1;

S = 0.1849;
L = 6.898;
epsilon = 3;

delta = 1e-2;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 20;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #9

% without switches (in P_1 no solution!)

k_1 = -1;
k_2 = 1;

alpha = 10;

T = 1;

S = -0.1572;
L = 1.8788;
epsilon = 0.5;

delta = 1e-2;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 20;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #10

% one switch (in P_1 we had two switches!)

k_1 = -1;
k_2 = 1;

alpha = 10;

T = 1;

S = 3.4865;
L = 0.3001;
epsilon = 0.1;

delta = 1e-1;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 20;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #11

% without switches (in P_1 no solution!)

k_1 = -1;
k_2 = 1;

alpha = 10;

T = 1;

S = 0.4505;
L = 0.8178;
epsilon = 0.1;

delta = 1e-2;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 100;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = first_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end

%% System #12

% one switch (in P_1 we had two switches!)

k_1 = -0.5;
k_2 = 0.05;

alpha = 1;

T = 1;

S = 0.2276;
L = -0.3578;
epsilon = 0.1;

delta = 1e-3;

number_of_points_for_splitting = 100;
number_of_points_for_enum = 20;

trimmer = 0.1;

%%

[u_opt, x_opt, psi_opt, functional_min, switches] = second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, delta, ...
    number_of_points_for_splitting, number_of_points_for_enum, trimmer);

if (~isempty(u_opt))
    
    disp(['Optimal value of functional is ', num2str(functional_min), '!']);
    
    draw_conjugate(psi_opt, switches);
    draw_traj(x_opt, switches, S, L, epsilon);

    draw_psi1_t(psi_opt, switches);
    draw_psi2_t(psi_opt, switches);

    draw_x1_t(x_opt, switches);
    draw_x2_t(x_opt, switches);

    draw_u1_t(u_opt, switches);
    draw_u2_t(u_opt, switches);

end