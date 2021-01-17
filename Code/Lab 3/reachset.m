function [X, Y, switches_line] = reachset(alpha, t, N, plot_traj, dist, K)

X = [];
Y = [];

fzero_opt = optimset('Display', 'none');

switches_line = [];

x_2_zero_opts = odeset('Events', @event_func_x_2, 'RelTol', 2.22045e-14);
psi_2_zero_opts = odeset('Events', @event_func_psi_2, 'RelTol', 2.22045e-14);

% S+

[times, x_start_plus] = ode45(@(t, x) odefun_S_plus(t, x, alpha), [0, t], ...
    [0, 0], x_2_zero_opts);

switches_line = [switches_line; x_start_plus];

if (N > 0)

    x_start_plus_fun = @(t) interp1(times, x_start_plus, t, 'spline');

    time_switch = times(end);

    time_splitting = linspace(0, time_switch, N);

    for index = 1 : N
        time = time_splitting(index);
        cur_trajectory = x_start_plus_fun(time_splitting(1 : index));
        while time < t
            [cur_times, cur_x] = ode45(@(t, x) odefun_S_minus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), 1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            if (time >= t)
                break;
            end

            [cur_times, cur_x] = ode45(@(t, x) odefun_S_plus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), -1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            if (time < t)
                switches_line = [switches_line; cur_x(end, 1 : 2)];
            end
        end

%         if index == N
%             switches_line = cur_trajectory(end : -1 : 1, 1 : 2);
%         end
        
        if (plot_traj)
            plot(cur_trajectory(:, 1), cur_trajectory(:, 2), 'm');
        end
        
        X = [X, cur_trajectory(end, 1)];
        Y = [Y, cur_trajectory(end, 2)];

    end
    
    while true
        index = 1;
        while true

            index = index + 1;

            cur_dist = norm([X(index) - X(index - 1), Y(index) - Y(index - 1)]);

            if (cur_dist > dist)
                
                X_cur = [];
                Y_cur = [];

                t1_cur = time_splitting(index - 1);
                t2_cur = time_splitting(index);

                time_splitting_cur = linspace(t1_cur, t2_cur, K);

                for index_cur = 1 : K
                    time = time_splitting_cur(index_cur);
                    cur_trajectory = x_start_plus_fun(time_splitting_cur(1 : index_cur));
                    while time < t
                        [cur_times, cur_x] = ode45(@(t, x) odefun_S_minus_conj(t, x, alpha), ...
                            [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), 1, 0], ...
                            psi_2_zero_opts);

                        cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

                        time = cur_times(end);
                        if (time >= t)
                            break;
                        end

                        [cur_times, cur_x] = ode45(@(t, x) odefun_S_plus_conj(t, x, alpha), ...
                            [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), -1, 0], ...
                            psi_2_zero_opts);

                        cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

                        time = cur_times(end);
                        
                        if (time < t)
                            switches_line = [switches_line; cur_x(end, 1 : 2)];
                        end
                    end

%                     if (index == N)
%                         switches_line = [switches_line; cur_trajectory];
%                     end

                    if (plot_traj)
                        plot(cur_trajectory(:, 1), cur_trajectory(:, 2), 'm');
                    end
                    
                    X_cur = [X_cur; cur_trajectory(end, 1)];
                    Y_cur = [Y_cur; cur_trajectory(end, 2)];
                end
                
                X = [X(1 : (index - 1)), X_cur.', X(index : end)];
                Y = [Y(1 : (index - 1)), Y_cur.', Y(index : end)];
                time_splitting = [time_splitting(1 : (index - 1)), ...
                    time_splitting_cur, time_splitting(index : end)];
                break;
            end
            
            if (index == length(X))
                break;
            end

        end
        if index == length(X) 
            break;
        end
    end
    
else
    for index = 1 : length(times)
        time = times(index);
        cur_trajectory = x_start_plus(1 : index, :);
        while time < t
            [cur_times, cur_x] = ode45(@(t, x) odefun_S_minus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), 1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            if (time >= t)
                break;
            end

            [cur_times, cur_x] = ode45(@(t, x) odefun_S_plus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), -1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            
            if (time < t)
                switches_line = [switches_line; cur_x(end, 1 : 2)];
            end
        end

%         if index == length(times)
%             switches_line = cur_trajectory(end : -1 : 1, 1 : 2);
%         end
        
        if (plot_traj)
            plot(cur_trajectory(:, 1), cur_trajectory(:, 2), 'm');
        end
        
        X = [X, cur_trajectory(end, 1)];
        Y = [Y, cur_trajectory(end, 2)];

    end
    
end

[times, x_start_minus] = ode45(@(t, x) odefun_S_minus(t, x, alpha), [0, t], ...
    [0, 0], x_2_zero_opts);

switches_line = [switches_line; x_start_minus];

if (N > 0)
    
    X_minus = [];
    Y_minus = [];

    x_start_minus_fun = @(t) interp1(times, x_start_minus, t, 'spline');

    time_switch = times(end);
    time_splitting = linspace(0, time_switch, N);

    for index = 1 : N
        time = time_splitting(index);
        cur_trajectory = x_start_minus_fun(time_splitting(1 : index));
        while time < t
            [cur_times, cur_x] = ode45(@(t, x) odefun_S_plus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), -1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            if (time >= t)
                break;
            end

            [cur_times, cur_x] = ode45(@(t, x) odefun_S_minus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), 1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            
            if (time < t)
                switches_line = [switches_line; cur_x(end, 1 : 2)];
            end
        end

%         if (index == N)
%             switches_line = [switches_line; cur_trajectory];
%         end

        if (plot_traj)
            plot(cur_trajectory(:, 1), cur_trajectory(:, 2), 'm');
        end
        
        X_minus = [X_minus, cur_trajectory(end, 1)];
        Y_minus = [Y_minus, cur_trajectory(end, 2)];

    end
    
    while true
        index = 1;
        while true

            index = index + 1;

            cur_dist = norm([X_minus(index) - X_minus(index - 1), Y_minus(index) - Y_minus(index - 1)]);

            if (cur_dist > dist)
                
                X_cur = [];
                Y_cur = [];

                t1_cur = time_splitting(index - 1);
                t2_cur = time_splitting(index);

                time_splitting_cur = linspace(t1_cur, t2_cur, K);

                for index_cur = 1 : K
                    time = time_splitting_cur(index_cur);
                    cur_trajectory = x_start_minus_fun(time_splitting_cur(1 : index_cur));
                    while time < t
                        [cur_times, cur_x] = ode45(@(t, x) odefun_S_plus_conj(t, x, alpha), ...
                            [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), -1, 0], ...
                            psi_2_zero_opts);

                        cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

                        time = cur_times(end);
                        if (time >= t)
                            break;
                        end

                        [cur_times, cur_x] = ode45(@(t, x) odefun_S_minus_conj(t, x, alpha), ...
                            [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), 1, 0], ...
                            psi_2_zero_opts);

                        cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

                        time = cur_times(end);
                        
                        if (time < t)
                            switches_line = [switches_line; cur_x(end, 1 : 2)];
                        end
                    end

%                     if (index == N)
%                         switches_line = [switches_line; cur_trajectory];
%                     end

                    if (plot_traj)
                        plot(cur_trajectory(:, 1), cur_trajectory(:, 2), 'm');
                    end
                    
                    X_cur = [X_cur; cur_trajectory(end, 1)];
                    Y_cur = [Y_cur; cur_trajectory(end, 2)];
                end
                
                X_minus = [X_minus(1 : (index - 1)), X_cur.', X_minus(index : end)];
                Y_minus = [Y_minus(1 : (index - 1)), Y_cur.', Y_minus(index : end)];
                time_splitting = [time_splitting(1 : (index - 1)), ...
                    time_splitting_cur, time_splitting(index : end)];
                
                break;
            end
            
            if (index == length(X_minus))
                break;
            end

        end
        if index == length(X_minus) 
            break;
        end
    end
    
    X = [X, X_minus];
    Y = [Y, Y_minus];
    
else
    for index = 1 : length(times)
        time = times(index);
        cur_trajectory = x_start_minus(1 : index, :);
        while time < t
            [cur_times, cur_x] = ode45(@(t, x) odefun_S_plus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), -1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            if (time >= t)
                break;
            end

            [cur_times, cur_x] = ode45(@(t, x) odefun_S_minus_conj(t, x, alpha), ...
                [time, t], [cur_trajectory(end, 1), cur_trajectory(end, 2), 1, 0], ...
                psi_2_zero_opts);

            cur_trajectory = [cur_trajectory; cur_x(2 : end, 1 : 2)];

            time = cur_times(end);
            
            if (time < t)
                switches_line = [switches_line; cur_x(end, 1 : 2)];
            end
        end

%         if (index == length(times))
%             switches_line = [switches_line; cur_trajectory];
%         end

        if (plot_traj)
            plot(cur_trajectory(:, 1), cur_trajectory(:, 2), 'm');
        end
        
        X = [X, cur_trajectory(end, 1)];
        Y = [Y, cur_trajectory(end, 2)];

    end
    
end

x2_spec = 0;
func_spec_1 = @(x1) - 2 * x1 - 3 * x1 * sin(x1 ^ 2) + 2*x1^2*cos(x1) + alpha;
func_spec_2 = @(x1) - 2 * x1 - 3 * x1 * sin(x1 ^ 2) + 2*x1^2*cos(x1) - alpha;

x1_min = min(X);
x1_max = max(X);

x1_split = linspace(x1_min, x1_max, N);

for x1_cur_0 = x1_split
    x1_zero = fzero(func_spec_1, x1_cur_0, fzero_opt);
    
    plot(x1_zero, x2_spec, 'og', 'LineWidth', 2);
end

for x1_cur_0 = x1_split
    x1_zero = fzero(func_spec_2, x1_cur_0, fzero_opt);
    
    plot(x1_zero, x2_spec, 'om', 'LineWidth', 2);
end

%[X, Y] = delete_cross(X, Y);
end