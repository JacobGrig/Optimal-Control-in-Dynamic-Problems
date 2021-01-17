function [t1_min, u_opt, x_opt, ...
    phi_opt, t_traj_opt, trans, psi_opt, ...
    x1_curs, x2_curs, u1_curs, u2_curs, ...
    psi1_curs, psi2_curs, t_psi_opt, ...
    psi_cur_opt, psi_trans_opt] = main_algorithm( ...
    A_matr, B_matr, ...                                   % parameters of initial task
    t1_max, eps_X1, eps_B, eps_rel, eps_abs, ...          % user's parameters
    a, b, c, r, d, k, ...                                 % parameters of the sets P, X0, X1
    phi_interval, number_of_points_for_splitting_phi) ... % auxiliary parameters

u_opt = [];
x_opt = [];
phi_opt = Inf;
t_traj_opt = [];
psi_opt = [];
t_psi_opt = [];
psi_cur_opt = [];
psi_trans_opt = [];

if ((abs(d - k) <= r) || ...
    (abs(d + k) <= r) || ...
    ((d + k > r) && (d - k < -r)))
    disp('The sets have non-empty intersection!');
    x1_curs = {};
    x2_curs = {};
    u1_curs = {};
    u2_curs = {};
    psi1_curs = {};
    psi2_curs = {};
    t1_min = 0;
    trans = 1;
    return;
end

% initialization
t0 = 0;
trans = -1;
x1_curs = cell(1, number_of_points_for_splitting_phi - 1);
x2_curs = cell(1, number_of_points_for_splitting_phi - 1);
u1_curs = cell(1, number_of_points_for_splitting_phi - 1);
u2_curs = cell(1, number_of_points_for_splitting_phi - 1);
psi1_curs = cell(1, number_of_points_for_splitting_phi - 1);
psi2_curs = cell(1, number_of_points_for_splitting_phi - 1);
phi_splitting = linspace(phi_interval(1), phi_interval(2), number_of_points_for_splitting_phi);

% calculation of the shift of the B_matr for a given eps_B (<=> calculating delta)
[~, fund] = ode45(@(t, fund) ode_fund(t, fund, A_matr), [t0, t1_max], [1; 0; 0; 1]);
fund_max = max(abs(fund), [], 1);
fund_max_sum = sqrt(2) .* max(fund_max(1) + fund_max(2), fund_max(3) + fund_max(4));
u_max = max(b, sqrt(b ./ (a + c)));
delta = eps_B ./ (u_max .* fund_max_sum .* t1_max);

% setting the event handler and tolerances
event_hndl = @(t, x) event_X1(t, x, d, k, eps_X1);
opts_traj = odeset('Events', event_hndl, 'RelTol', eps_rel, 'AbsTol', eps_abs);
opts_conj = odeset('RelTol', eps_rel, 'AbsTol', eps_abs);

% searching in all directions of psi_0
for i = 1 : (number_of_points_for_splitting_phi - 1)
    psi_0 = [cos(phi_splitting(i)); sin(phi_splitting(i))];
    
    % using the first transversality condition
    [~, x_0_cur] = support_X0(psi_0, r);
    
    % solving the conjugate system and interpolating
    [t_psi, psi] = ode45(@(t, psi) -A_matr(t).' * psi, [t0, t1_max], psi_0, opts_conj);
    psi1_hndl = @(t) interp1(t_psi, psi(:, 1), t, 'spline');
    psi2_hndl = @(t) interp1(t_psi, psi(:, 2), t, 'spline');
    psi_hndl = @(t) [psi1_hndl(t); psi2_hndl(t)];
    
    psi1_curs{i} = [t_psi.'; psi(:, 1).'];
    psi2_curs{i} = [t_psi.'; psi(:, 2).'];
    
%     psi1_curs{i} = psi1_hndl;
%     psi2_curs{i} = psi2_hndl;
    
    % finding the sub-optimal control, using the maximum condition
    len_psi = length(t_psi);
    u_cur = zeros(2, len_psi);
    for j = 1 : len_psi
        B_matr_cur = B_matr(t_psi(j));
        psi_val = psi(j, :);
        indicator = B_matr_cur.' * psi_val.';
        
        % defining, if we can unquely find the argmax of support_P
        if ((indicator(1) ~= 0) || (indicator(2) ~= 0))
            [~, u_cur(:, j)] = support_P(indicator, a, b, c);
        else
            [~, u_cur(:, j)] = support_P(delta .* psi_val.', a, b, c);
        end
    end

    % interpolating founded sub-optimal control to any point
    u1_cur_hndl = @(t) interp1(t_psi, u_cur(1, :), t, 'spline');
    u2_cur_hndl = @(t) interp1(t_psi, u_cur(2, :), t, 'spline');
    u_cur_hndl = @(t) [u1_cur_hndl(t); u2_cur_hndl(t)];

    u1_curs{i} = [t_psi.'; u_cur(1, :)];
    u2_curs{i} = [t_psi.'; u_cur(2, :)];
    
%     u1_curs{i} = u1_cur_hndl;
%     u2_curs{i} = u2_cur_hndl;
    
    % solving the ODE of initial task
    [t_traj, x_cur, t_event, x_cur_event, ie] = ode45( ...
        @(t, x) A_matr(t) * x + B_matr(t) * u_cur_hndl(t), [0, t1_max], x_0_cur, opts_traj);
    
    x1_curs{i} = [t_traj.'; x_cur(:, 1).'];
    x2_curs{i} = [t_traj.'; x_cur(:, 2).'];
    
%     x1_curs{i} = x1_cur_hndl;
%     x2_curs{i} = x2_cur_hndl;
    
    % defining, if we arrived to X1 faster, than earlier
    if (ie)
        if (t_event < t1_max)
            t1_max = t_event;
            u_opt = u_cur;
            x_opt = x_cur;
            t_traj_opt = t_traj;
            t_psi_opt = t_psi;
            phi_opt = phi_splitting(i);
            psi_opt = psi;

            % checking the second transversality condition
            psi_cur = -psi_hndl(t_event);
            if (x_cur_event(2) < d)
                psi_trans = sqrt(1 ./ (1 + 4 .* (x_cur_event(1) .^ 2))) .* ...
                    [2 .* x_cur_event(1); -1];
            else
            if (x_cur_event(2) > d)
                psi_trans = sqrt(1 ./ (1 + 4 .* (x_cur_event(1) .^ 2))) .* ...
                    [2 .* x_cur_event(1); 1];
            else
            if ((abs(psi_cur(2)) < abs(psi_cur(1)) / (2 .* sqrt(k))) && ...
                    sign(psi_cur(1)) == sign(x_cur_event(1)))
                psi_trans = psi_cur;
            else
                psi_trans_tmp = [sign(x_cur_event(1)) .* 2 .* sqrt(k), sign(psi_cur(2))];
                psi_trans = psi_trans_tmp ./ ...
                    sqrt(psi_trans_tmp(1) .^ 2 + psi_trans_tmp(2) .^ 2);
            end
            end
            end
            trans = dot(psi_cur, psi_trans) ./ sqrt(psi_cur(1) .^ 2 + psi_cur(2) .^ 2);
            psi_trans_opt = psi_trans;
            psi_cur_opt = psi_cur;
        end
    end
end
t1_min = t1_max;
if (~isempty(u_opt))
    
    disp(['cos(psi(t1), psi^1) = ', num2str(trans)]);
    disp(['t1_min = ', num2str(t1_min)]);
    disp(' ');
end
end