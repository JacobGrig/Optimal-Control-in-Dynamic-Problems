function [u_opt, x_opt, ...
    psi_opt, functional_min, ...
    switches] ...
= second_solution ( ...
    T, k_1, k_2, L, S, epsilon, alpha, ...                            % parameters of initial task
    delta, ...                                                        % user's parameters
    number_of_points_for_splitting, number_of_points_for_enum, ...
    trimmer) ...                                                      % auxiliary parameters
    
functional_min = Inf;

switches = [];

fsolve_opts = optimoptions('fsolve', 'Display','none');

t_opt = linspace(0, T, number_of_points_for_splitting);

% trying to get an optimal solution

f = @(k) k .* T + alpha ./ 2 .* (k .* T - exp(k .* T) + 1) ./ (k .^ 2);
g = @(k) -alpha ./ 2 .* (exp(k .* T) - 1) ./ k;

k_splitting = linspace(k_1, k_2, number_of_points_for_enum);

found = false;

for k_cur = k_splitting
    if (S - epsilon <= g(k_cur) && ...
        g(k_cur) <= S + epsilon && ...
        L - epsilon <= f(k_cur) && ...
        f(k_cur) <= L + epsilon)
        
        found = true;
        
        u1_opt = -alpha ./ 2 .* ones(1, number_of_points_for_splitting);
        x_2_t = @(t) -alpha ./ 2 .* (exp(k_cur .* t) - 1) ./ k_cur;
        x_1_t = @(t) k_cur .* t + alpha ./ 2 .* (k_cur .* t - exp(k_cur .* T) + 1) ./ (k_cur .^ 2);
        
        x2_opt = x_2_t(t_opt);
        x1_opt = x_1_t(t_opt);
        
        u2_opt = k_cur .* ones(1, number_of_points_for_splitting);

        psi1_opt = zeros(1, number_of_points_for_splitting);
        psi2_opt = zeros(1, number_of_points_for_splitting);
        
        functional_min = -alpha .^ 2 ./ 4 .* T;
        
    end
end

if (~found)
    
    % first case: psi_1 = 0, psi_2_0 (only without switches)
    
    u2_cur = k_2;
    
    x_1_t_0 = @(t, psi_2_0) k_2.*t + (alpha./2 - psi_2_0./2 - (alpha.*exp(k_2.*t))./2 + (psi_2_0.*exp(-k_2.*t))./2 + k_2.*((alpha.*t)./2 + (psi_2_0.*t)./2))./k_2.^2;
    x_2_t_0 = @(t, psi_2_0) (exp(-k_2.*t).*(exp(k_2.*t) - 1).*(psi_2_0 - alpha.*exp(k_2.*t)))./(2.*k_2);
    
    psi_2_t_0 = @(t, psi_2_0) psi_2_0.*exp(-k_2.*t);
    
    u_1_t_0 = @(t, psi_2_0) (psi_2_t_0(t, psi_2_0) - alpha) ./ 2;
    
    psi_2_0_cur = (2.*k_2.*exp(T.*k_2).*(S + epsilon + (alpha.*(exp(T.*k_2) - 1))./(2.*k_2)))./(exp(T.*k_2) - 1);
    
    if (psi_2_0_cur < 0)
    
        psi_2_t = @(t) psi_2_t_0(t, psi_2_0_cur);

        x_1_t = @(t) x_1_t_0(t, psi_2_0_cur);
        x_2_t = @(t) x_2_t_0(t, psi_2_0_cur);
        
        if (abs(L - x_1_t(T)) < epsilon + delta)
            u_1_t = @(t) u_1_t_0(t, psi_2_0_cur);
            int_func = @(t) u_1_t(t) .^ 2 + alpha .* u_1_t(t);
            
            functional_cur = integral(int_func, 0, T);
            
            if (functional_cur < functional_min)
                psi1_opt = zeros(1, number_of_points_for_splitting);
                psi2_opt = psi_2_t(t_opt);
                
                x1_opt = x_1_t(t_opt);
                x2_opt = x_2_t(t_opt);
                
                u1_opt = u_1_t(t_opt);
                u2_opt = u2_cur .* ones(1, number_of_points_for_splitting);
                
                functional_min = functional_cur;
            end
        end
    end
    
    % second case: psi_1_0 = 0, psi_2_0 > alpha
    
    psi1_cur = 0;
    u2_first = k_2;

    psi_2_t_0 = @(t, psi_2_0) exp(-k_2 .* t) .* psi_2_0;

    u_1_t_0 = @(t, psi_2_0) (psi_2_t_0(t, psi_2_0) - alpha) ./ 2;

    x_2_t_0_first = @(t, psi_2_0) (psi_2_0.*(exp(k_2.*t)./2 - exp(-k_2.*t)./2))./(2.*k_2) - (alpha.*(exp(k_2.*t) - 1))./(2.*k_2);

    % without switches

    psi_2_0_cur = fsolve(@(psi_2_0) x_2_t_0_first(T, psi_2_0) - (S - epsilon), alpha, fsolve_opts);

    if (psi_2_0_cur > alpha)
        u_1_t = @(t) u_1_t_0(t, psi_2_0_cur);

        x_2_t = @(t) x_2_t_0_first(t, psi_2_0_cur);
        
        if (prod(x_2_t(t_opt) >= 0))

            int_func = @(t) u_1_t(t) .^ 2 + alpha .* u_1_t(t);

            functional_cur = integral(int_func, 0, T);

            x_1_t = @(t) - ((exp(-k_2.*t).*(alpha.*exp(T.*k_2) - alpha.*exp(2.*T.*k_2) + alpha.*exp(k_2.*t) - alpha.*exp(2.*k_2.*t) - 2.*alpha.*exp(T.*k_2).*exp(k_2.*t) + alpha.*exp(T.*k_2).*exp(2.*k_2.*t) + alpha.*exp(2.*T.*k_2).*exp(k_2.*t)))./2 - (k_2.*exp(-k_2.*t).*(2.*S.*exp(T.*k_2) - 2.*epsilon.*exp(T.*k_2) - alpha.*t.*exp(k_2.*t) - 4.*S.*exp(T.*k_2).*exp(k_2.*t) + 2.*S.*exp(T.*k_2).*exp(2.*k_2.*t) + 4.*epsilon.*exp(T.*k_2).*exp(k_2.*t) - 2.*epsilon.*exp(T.*k_2).*exp(2.*k_2.*t) + alpha.*t.*exp(2.*T.*k_2).*exp(k_2.*t)))./2)./(k_2.^2.*(exp(2.*T.*k_2) - 1)) - (k_2.*exp(-k_2.*t).*(2.*t.*exp(k_2.*t) - 2.*t.*exp(2.*T.*k_2).*exp(k_2.*t)))./(2.*(exp(2.*T.*k_2) - 1));

            if (functional_cur < functional_min && ... 
                abs(L - x_1_t(T)) <= epsilon + delta)

                psi_2_t = @(t) psi_2_t_0(t, psi_2_0_cur);

                psi1_opt = psi1_cur .* ones(1, number_of_points_for_splitting);
                psi2_opt = psi_2_t(t_opt);

                x1_opt = x_1_t(t_opt);
                x2_opt = x_2_t(t_opt);

                u1_opt = u_1_t(t_opt);
                u2_opt = u2_first .* ones(1, number_of_points_for_splitting);

                functional_min = functional_cur;
            end
        end
    end
    
    % one switch
    
    tau_cur = log((exp(-T.*k_2).*(2.*S.*k_2 - 2.*epsilon.*k_2 + alpha.*exp(T.*k_2) - 2.*(k_2.*(S - epsilon).*(S.*k_2 - epsilon.*k_2 + alpha.*exp(T.*k_2))).^(1./2)))./alpha)./k_2;

    psi_2_0_cur_tau1 = @(tau1) (2.*alpha.*exp(k_2.*tau1))./(exp(k_2.*tau1) + 1);
    
    psi_2_0_cur = psi_2_0_cur_tau1(tau_cur);
    
    if (psi_2_0_cur > alpha && ...
        0 < tau_cur && ...
        tau_cur < T)
        
        [t_cur, x_cur] = ode45(@(t, x) odefun_second(t, x, k_1, k_2, alpha, 0), [0, T], [0, 0, psi_2_0_cur]);
        x_1_vect = x_cur(:, 1);
        x_2_vect = x_cur(:, 2);
        psi_2_vect = x_cur(:, 3);
        if (abs(L - x_1_vect(end)) <= epsilon + delta)

            u_1_vect = (psi_2_vect - alpha) ./ 2;

            functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

            if (functional_cur < functional_min)

                u_2_vect = (psi_1_cur_1 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                           (psi_1_cur_1 + psi_2_vect .* x_2_vect < 0) .* k_1;

                t_opt = t_cur.';

                psi1_opt = psi_1_cur_1 .* ones(1, length(t_cur));
                psi2_opt = psi_2_vect.';

                x1_opt = x_1_vect.';
                x2_opt = x_2_vect.';

                u1_opt = u_1_vect.';
                u2_opt = u_2_vect.';

                functional_min = functional_cur;

                switches = [tau_cur];
            end
        end
    end
    
    % third case: psi_1 = 0, 0 < psi_2_0 < alpha
    
    psi1_cur = 0;
    
    u2_first = k_1;

    psi_2_t_0 = @(t, psi_2_0) exp(-k_1 .* t) .* psi_2_0;

    u_1_t_0 = @(t, psi_2_0) (psi_2_t_0(t, psi_2_0) - alpha) ./ 2;

    x_2_t_0_first = @(t, psi_2_0) (psi_2_0.*(exp(k_1.*t)./2 - exp(-k_1.*t)./2))./(2.*k_1) - (alpha.*(exp(k_1.*t) - 1))./(2.*k_1);

    % without switches

    psi_2_0_cur = fsolve(@(psi_2_0) x_2_t_0_first(T, psi_2_0) - (S - epsilon), alpha, fsolve_opts);

    if (psi_2_0_cur < alpha)
        u_1_t = @(t) u_1_t_0(t, psi_2_0_cur);

        x_2_t = @(t) x_2_t_0_first(t, psi_2_0_cur);
        
        if (prod(x_2_t(t_opt) <= 0))

            int_func = @(t) u_1_t(t) .^ 2 + alpha .* u_1_t(t);

            functional_cur = integral(int_func, 0, T);

            x_1_t = @(t) - ((exp(-k_1.*t).*(alpha.*exp(T.*k_1) - alpha.*exp(2.*T.*k_1) + alpha.*exp(k_1.*t) - alpha.*exp(2.*k_1.*t) - 2.*alpha.*exp(T.*k_1).*exp(k_1.*t) + alpha.*exp(T.*k_1).*exp(2.*k_1.*t) + alpha.*exp(2.*T.*k_1).*exp(k_1.*t)))./2 - (k_1.*exp(-k_1.*t).*(2.*S.*exp(T.*k_1) - 2.*epsilon.*exp(T.*k_1) - alpha.*t.*exp(k_1.*t) - 4.*S.*exp(T.*k_1).*exp(k_1.*t) + 2.*S.*exp(T.*k_1).*exp(2.*k_1.*t) + 4.*epsilon.*exp(T.*k_1).*exp(k_1.*t) - 2.*epsilon.*exp(T.*k_1).*exp(2.*k_1.*t) + alpha.*t.*exp(2.*T.*k_1).*exp(k_1.*t)))./2)./(k_1.^2.*(exp(2.*T.*k_1) - 1)) - (k_1.*exp(-k_1.*t).*(2.*t.*exp(k_1.*t) - 2.*t.*exp(2.*T.*k_1).*exp(k_1.*t)))./(2.*(exp(2.*T.*k_1) - 1));

            if (functional_cur < functional_min && ... 
                abs(L - x_1_t(T)) <= epsilon + delta)
            
                t_opt = linspace(0, T, number_of_points_for_splitting);

                psi_2_t = @(t) psi_2_t_0(t, psi_2_0_cur);

                psi1_opt = psi1_cur .* ones(1, number_of_points_for_splitting);
                psi2_opt = psi_2_t(t_opt);

                x1_opt = x_1_t(t_opt);
                x2_opt = x_2_t(t_opt);

                u1_opt = u_1_t(t_opt);
                u2_opt = u2_first .* ones(1, number_of_points_for_splitting);

                functional_min = functional_cur;
            end
        end
    end
    
    % one switch
    
    tau_cur = log((exp(-T.*k_1).*(2.*S.*k_1 - 2.*epsilon.*k_1 + alpha.*exp(T.*k_1) - 2.*(k_1.*(S - epsilon).*(S.*k_1 - epsilon.*k_1 + alpha.*exp(T.*k_1))).^(1./2)))./alpha)./k_1;

    psi_2_0_cur_tau1 = @(tau1) (2.*alpha.*exp(k_1.*tau1))./(exp(k_1.*tau1) + 1);
    
    psi_2_0_cur = psi_2_0_cur_tau1(tau_cur);
    
    if (psi_2_0_cur < alpha && ...
        0 < tau_cur && ...
        tau_cur < T)
        
        [t_cur, x_cur] = ode45(@(t, x) odefun_second(t, x, k_2, k_1, alpha, 0), [0, T], [0, 0, psi_2_0_cur]);
        x_1_vect = x_cur(:, 1);
        x_2_vect = x_cur(:, 2);
        psi_2_vect = x_cur(:, 3);
        if (abs(L - x_1_vect(end)) <= epsilon + delta)

            u_1_vect = (psi_2_vect - alpha) ./ 2;

            functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

            if (functional_cur < functional_min)

                u_2_vect = (psi_1_cur_1 + psi_2_vect .* x_2_vect > 0) .* k_1 + ...
                           (psi_1_cur_1 + psi_2_vect .* x_2_vect < 0) .* k_2;

                t_opt = t_cur.';

                psi1_opt = psi_1_cur_1 .* ones(1, length(t_cur));
                psi2_opt = psi_2_vect.';

                x1_opt = x_1_vect.';
                x2_opt = x_2_vect.';

                u1_opt = u_1_vect.';
                u2_opt = u_2_vect.';

                functional_min = functional_cur;

                switches = [tau_cur];
            end
        end
    end
    
    % fourth case: psi_1 = 0, psi_2_0 == alpha

    if (functional_min > 0 && abs(S - epsilon) <= delta && ...
            T >= max((L + epsilon + delta) ./ k_1, (L - epsilon - delta) ./ k_2))
        
        t_opt = linspace(0, T, number_of_points_for_splitting);
    
        u1_opt = zeros(1, number_of_points_for_splitting);
        x2_opt = zeros(1, number_of_points_for_splitting);

        psi1_opt = zeros(1, number_of_points_for_splitting);
        psi2_opt = alpha .* ones(1, number_of_points_for_splitting);

        if (L - epsilon - delta <= k_1 .* T && k_1 .* T <= L + epsilon + delta)
            u2_opt = k_1 .* ones(1, number_of_points_for_splitting);
            x1_opt = k_1 .* t_opt;
        else
            if (L - epsilon - delta <= k_2 .* T && k_2 .* T <= L + epsilon + delta)
                u2_opt = k_2 .* ones(1, number_of_points_for_splitting);
                x1_opt = k_2 .* t_opt;
            else
                u2_opt = L ./ T .* ones(1, number_of_points_for_splitting);
                x1_opt = L ./ T .* t_opt;
            end
        end
        functional_min = 0;
    end
    
    % psi_1 = 0 is done!
    
    % fifth case: psi_1 < 0
    
    u2_first = k_1;
    
    psi_2_t_1_0_first = @(t, psi_1, psi_2_0) psi_2_0.*exp(-k_1.*t) + (psi_1.*(exp(-k_1.*t) - 1))./k_1;
    u_1_t_1_0_first = @(t, psi_1, psi_2_0) (psi_2_t_1_0_first(t, psi_1, psi_2_0) - alpha) ./ 2;
    
    x_1_t_1_0_first = @(t, psi_1, psi_2_0) (2.*alpha.*k_1 - 2.*k_1.*psi_2_0 - psi_1.*exp(k_1.*t) + psi_1.*exp(-k_1.*t) + 4.*k_1.^4.*t + 2.*k_1.*psi_1.*t - 2.*alpha.*k_1.*exp(k_1.*t) + k_1.*psi_2_0.*exp(k_1.*t) + k_1.*psi_2_0.*exp(-k_1.*t) + 2.*alpha.*k_1.^2.*t)./(4.*k_1.^3);
    x_2_t_1_0_first = @(t, psi_1, psi_2_0) (exp(-k_1.*t).*(exp(k_1.*t) - 1).*(psi_1 + k_1.*psi_2_0 - psi_1.*exp(k_1.*t) - 2.*alpha.*k_1.*exp(k_1.*t) + k_1.*psi_2_0.*exp(k_1.*t)))./(4.*k_1.^2);    
    
    % without switches
    
    % psi_2(T) < 0
    
    psi_1_cur = -(k_1.*(2.*alpha - 2.*S.*k_1 - 2.*epsilon.*k_1 - 2.*alpha.*exp(T.*k_1) - 2.*L.*k_1.^2 + 2.*T.*k_1.^3 - 2.*epsilon.*k_1.^2 + 2.*T.*k_1.^3.*exp(T.*k_1) - 2.*epsilon.*k_1.^2.*exp(T.*k_1) + T.*alpha.*k_1 + 2.*S.*k_1.*exp(T.*k_1) + 2.*epsilon.*k_1.*exp(T.*k_1) - 2.*L.*k_1.^2.*exp(T.*k_1) + T.*alpha.*k_1.*exp(T.*k_1)))./(T.*k_1 - 2.*exp(T.*k_1) + T.*k_1.*exp(T.*k_1) + 2);
    psi_2_0_cur = (2.*S.*k_1 - 2.*alpha + 2.*epsilon.*k_1 + 4.*alpha.*exp(T.*k_1) - 2.*alpha.*exp(2.*T.*k_1) + 2.*L.*k_1.^2 - 2.*T.*k_1.^3 + 2.*epsilon.*k_1.^2 + 4.*T.*k_1.^3.*exp(T.*k_1) - 2.*T.*k_1.^3.*exp(2.*T.*k_1) - 4.*epsilon.*k_1.^2.*exp(T.*k_1) + 2.*epsilon.*k_1.^2.*exp(2.*T.*k_1) - T.*alpha.*k_1 - 2.*S.*k_1.*exp(2.*T.*k_1) - 2.*epsilon.*k_1.*exp(2.*T.*k_1) - 4.*L.*k_1.^2.*exp(T.*k_1) + 2.*L.*k_1.^2.*exp(2.*T.*k_1) + T.*alpha.*k_1.*exp(2.*T.*k_1) + 4.*S.*T.*k_1.^2.*exp(T.*k_1) + 4.*T.*epsilon.*k_1.^2.*exp(T.*k_1))./((exp(T.*k_1) - 1).*(T.*k_1 - 2.*exp(T.*k_1) + T.*k_1.*exp(T.*k_1) + 2));
    
    if (psi_1_cur < 0)
        u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        psi_2_t_first = @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        
        if (psi_2_t_first(T) < 0)
        
            x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

            if (prod(psi_1_cur + psi_2_t_first(t_opt) .* x_2_t_first(t_opt) <= 0))
                int_func = @(t) u_1_t_first(t) .^ 2 + alpha .* u_1_t_first(t);

                functional_cur = integral(int_func, 0, T);

                if (functional_cur < functional_min)
                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_cur .* ones(1, number_of_points_for_splitting);
                    psi2_opt = psi_2_t_first(t_opt);

                    x1_opt = x_1_t_first(t_opt);
                    x2_opt = x_2_t_first(t_opt);

                    u1_opt = u_1_t_first(t_opt);
                    u2_opt = u2_first .* ones(1, number_of_points_for_splitting);

                    functional_min = functional_cur;
                end
            end
        end
    end
    
    % psi_2(T) > 0
    
    psi_1_cur = -(k_1.*(2.*alpha - 2.*S.*k_1 + 2.*epsilon.*k_1 - 2.*alpha.*exp(T.*k_1) - 2.*L.*k_1.^2 + 2.*T.*k_1.^3 - 2.*epsilon.*k_1.^2 + 2.*T.*k_1.^3.*exp(T.*k_1) - 2.*epsilon.*k_1.^2.*exp(T.*k_1) + T.*alpha.*k_1 + 2.*S.*k_1.*exp(T.*k_1) - 2.*epsilon.*k_1.*exp(T.*k_1) - 2.*L.*k_1.^2.*exp(T.*k_1) + T.*alpha.*k_1.*exp(T.*k_1)))./(T.*k_1 - 2.*exp(T.*k_1) + T.*k_1.*exp(T.*k_1) + 2);
    psi_2_0_cur = -(2.*alpha - 2.*S.*k_1 + 2.*epsilon.*k_1 - 4.*alpha.*exp(T.*k_1) + 2.*alpha.*exp(2.*T.*k_1) - 2.*L.*k_1.^2 + 2.*T.*k_1.^3 - 2.*epsilon.*k_1.^2 - 4.*T.*k_1.^3.*exp(T.*k_1) + 2.*T.*k_1.^3.*exp(2.*T.*k_1) + 4.*epsilon.*k_1.^2.*exp(T.*k_1) - 2.*epsilon.*k_1.^2.*exp(2.*T.*k_1) + T.*alpha.*k_1 + 2.*S.*k_1.*exp(2.*T.*k_1) - 2.*epsilon.*k_1.*exp(2.*T.*k_1) + 4.*L.*k_1.^2.*exp(T.*k_1) - 2.*L.*k_1.^2.*exp(2.*T.*k_1) - T.*alpha.*k_1.*exp(2.*T.*k_1) - 4.*S.*T.*k_1.^2.*exp(T.*k_1) + 4.*T.*epsilon.*k_1.^2.*exp(T.*k_1))./((exp(T.*k_1) - 1).*(T.*k_1 - 2.*exp(T.*k_1) + T.*k_1.*exp(T.*k_1) + 2));

    if (psi_1_cur < 0)
        u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        psi_2_t_first = @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        
        if (psi_2_t_first(T) > 0)
        
            x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

            if (prod(psi_1_cur + psi_2_t_first(t_opt) .* x_2_t_first(t_opt) <= 0))
                int_func = @(t) u_1_t_first(t) .^ 2 + alpha .* u_1_t_first(t);

                functional_cur = integral(int_func, 0, T);

                if (functional_cur < functional_min)
                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_cur .* ones(1, number_of_points_for_splitting);
                    psi2_opt = psi_2_t_first(t_opt);

                    x1_opt = x_1_t_first(t_opt);
                    x2_opt = x_2_t_first(t_opt);

                    u1_opt = u_1_t_first(t_opt);
                    u2_opt = u2_first .* ones(1, number_of_points_for_splitting);

                    functional_min = functional_cur;
                end
            end
        end
    end
    
    % one switch
    
    u2_second = k_2;
    
    psi_2_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) exp(-k_2.*(t - tau1)).*(psi_2_0.*exp(-k_1.*tau1) + (psi_1.*(exp(-k_1.*tau1) - 1))./k_1) + (psi_1.*(exp(-k_2.*(t - tau1)) - 1))./k_2;
    
    u_1_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) (psi_2_t_1_0_tau1_second(t, psi_1, psi_2_0, tau1) - alpha) ./ 2;
    
    x_2_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) (psi_1.*exp(- k_2.*t - k_2.*tau1).*(exp(2.*k_2.*t) - exp(2.*k_2.*tau1)))./(4.*k_2.^2) - (psi_1.*(exp(k_2.*(t - tau1)) - 1))./(2.*k_2.^2) - (alpha.*(exp(k_2.*(t - tau1)) - 1))./(2.*k_2) + (psi_2_0.*exp(- k_2.*t - k_2.*tau1).*exp(-k_1.*tau1).*(exp(2.*k_2.*t) - exp(2.*k_2.*tau1)))./(4.*k_2) - (psi_1.*exp(- k_2.*t - k_2.*tau1).*(exp(2.*k_2.*t) - exp(2.*k_2.*tau1)))./(4.*k_1.*k_2) + (exp(-k_1.*tau1).*exp(k_2.*(t - tau1)).*(exp(k_1.*tau1) - 1).*(psi_1 + k_1.*psi_2_0 - psi_1.*exp(k_1.*tau1) - 2.*alpha.*k_1.*exp(k_1.*tau1) + k_1.*psi_2_0.*exp(k_1.*tau1)))./(4.*k_1.^2) + (psi_1.*exp(- k_2.*t - k_2.*tau1).*exp(-k_1.*tau1).*(exp(2.*k_2.*t) - exp(2.*k_2.*tau1)))./(4.*k_1.*k_2);

    x_1_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) k_1.*tau1 + k_2.*(t - tau1) + (alpha.*tau1)./(2.*k_1) + (psi_1.*tau1)./(2.*k_1.^2) + (alpha.*(t - tau1))./(2.*k_2) + (psi_1.*(t - tau1))./(2.*k_2.^2) + (psi_1 - psi_1.*exp(-k_2.*(t - tau1)))./(4.*k_1.*k_2.^2) - (alpha.*(exp(k_1.*tau1) - 1))./(2.*k_1.^2) - (psi_1.*(exp(k_1.*tau1) - 1))./(4.*k_1.^3) + (psi_1.*(exp(-k_1.*tau1) - 1))./(4.*k_1.^3) + (psi_2_0.*(exp(k_1.*tau1) - 1))./(4.*k_1.^2) + (psi_2_0.*(exp(-k_1.*tau1) - 1))./(4.*k_1.^2) - (alpha.*(exp(k_2.*(t - tau1)) - 1))./(2.*k_2.^2) - (psi_1.*(exp(k_2.*(t - tau1)) - 1))./(4.*k_2.^3) + (psi_1.*exp(2.*k_2.*tau1).*(exp(-k_2.*(t + tau1)) - exp(-2.*k_2.*tau1)))./(4.*k_2.^3) + (psi_2_0.*exp(-tau1.*(k_1 - 2.*k_2)).*(exp(-k_2.*(t + tau1)) - exp(-2.*k_2.*tau1)))./(4.*k_2.^2) + (psi_2_0.*exp(-k_1.*tau1).*(exp(k_2.*(t - tau1)) - 1))./(4.*k_2.^2) + (alpha.*(exp(k_2.*(t - tau1)) - 1))./(2.*k_1.*k_2) - (psi_1.*(exp(k_2.*(t - tau1)) - 1))./(4.*k_1.*k_2.^2) + (psi_1.*(exp(k_2.*(t - tau1)) - 1))./(2.*k_1.^2.*k_2) + (psi_1.*exp(-tau1.*(k_1 - 2.*k_2)).*(exp(-k_2.*(t + tau1)) - exp(-2.*k_2.*tau1)))./(4.*k_1.*k_2.^2) - (alpha.*exp(k_1.*tau1).*(exp(k_2.*(t - tau1)) - 1))./(2.*k_1.*k_2) - (psi_1.*exp(k_1.*tau1).*(exp(k_2.*(t - tau1)) - 1))./(4.*k_1.^2.*k_2) + (psi_1.*exp(-k_1.*tau1).*(exp(k_2.*(t - tau1)) - 1))./(4.*k_1.*k_2.^2) - (psi_1.*exp(-k_1.*tau1).*(exp(k_2.*(t - tau1)) - 1))./(4.*k_1.^2.*k_2) + (psi_2_0.*exp(k_1.*tau1).*(exp(k_2.*(t - tau1)) - 1))./(4.*k_1.*k_2) - (psi_2_0.*exp(-k_1.*tau1).*(exp(k_2.*(t - tau1)) - 1))./(4.*k_1.*k_2);

    % psi_2(T) < 0
    
    psi_1_cur_tau1 = @(tau1) (2.*k_1.*k_2.*(4.*alpha.*k_1.^3.*exp(-k_1.*tau1) + 4.*alpha.*k_2.^3.*exp(k_2.*(T - tau1)) - 4.*S.*k_1.*k_2.^3 - 2.*alpha.*k_1.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*alpha.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.*k_2.^2 + 2.*alpha.*k_1.^2.*k_2 - 4.*epsilon.*k_1.*k_2.^3 + 2.*S.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + alpha.*k_1.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*alpha.*k_1.^2.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*epsilon.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*S.*k_1.*k_2.^3.*exp(k_1.*tau1) + 2.*S.*k_1.*k_2.^3.*exp(-k_1.*tau1) - 4.*S.*k_1.^3.*k_2.*exp(-k_1.*tau1) - alpha.*k_1.*k_2.^2.*exp(k_1.*tau1) + alpha.*k_1.^2.*k_2.*exp(k_1.*tau1) - alpha.*k_1.*k_2.^2.*exp(-k_1.*tau1) - 3.*alpha.*k_1.^2.*k_2.*exp(-k_1.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(k_1.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(-k_1.*tau1) - 4.*epsilon.*k_1.^3.*k_2.*exp(-k_1.*tau1) - 2.*L.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*L.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^2.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^2.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^3.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 3.*alpha.*k_1.*k_2.^2.*exp(k_2.*(T - tau1)) - alpha.*k_1.^2.*k_2.*exp(k_2.*(T - tau1)) + alpha.*k_1.*k_2.^2.*exp(-k_2.*(T - tau1)) - alpha.*k_1.^2.*k_2.*exp(-k_2.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.^2.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 4.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*k_1.^4.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.^2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*S.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 2.*S.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) + T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + T.*alpha.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - T.*alpha.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + alpha.*k_1.*k_2.^3.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - alpha.*k_1.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^3.*k_2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2)))./(4.*k_1.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 4.*k_1.^4.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 8.*k_1.^4.*exp(-k_1.*tau1) - 8.*k_2.^4.*exp(k_2.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.^3.*k_2.*exp(k_1.*tau1) + 4.*k_1.^3.*k_2.*exp(-k_1.*tau1) + 4.*k_1.*k_2.^3.*exp(k_2.*(T - tau1)) - 4.*k_1.*k_2.^3.*exp(-k_2.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*(T - tau1)) - k_1.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 3.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^4.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1));
    psi_2_0_cur_tau1 = @(tau1) -(2.*(4.*alpha.*k_1.^4 - 4.*L.*k_1.^4.*k_2.^2 - 2.*alpha.*k_1.^4.*exp(k_2.*(T - tau1)) - 2.*alpha.*k_1.^4.*exp(-k_2.*(T - tau1)) + 4.*alpha.*k_2.^4.*exp(k_2.*(T - tau1)) - 4.*S.*k_1.^2.*k_2.^3 + 4.*S.*k_1.^3.*k_2.^2 + 4.*T.*k_1.^4.*k_2.^3 + 8.*alpha.*k_1.^2.*k_2.^2 - 4.*epsilon.*k_1.^2.*k_2.^3 + 4.*epsilon.*k_1.^3.*k_2.^2 - 4.*epsilon.*k_1.^4.*k_2.^2 - 4.*k_1.^4.*k_2.^3.*tau1 + 4.*k_1.^5.*k_2.^2.*tau1 - 2.*alpha.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 8.*alpha.*k_1.^3.*k_2 + 4.*T.*k_1.^2.*k_2.^5.*exp(k_2.*(T - tau1)) - 2.*T.*k_1.^3.*k_2.^4.*exp(k_2.*(T - tau1)) - 2.*T.*k_1.^4.*k_2.^3.*exp(k_2.*(T - tau1)) + 2.*T.*k_1.^3.*k_2.^4.*exp(-k_2.*(T - tau1)) - 2.*T.*k_1.^4.*k_2.^3.*exp(-k_2.*(T - tau1)) - 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_2.*(T - tau1)) - 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_2.*(T - tau1)) + 4.*epsilon.*k_1.^2.*k_2.^3.*exp(k_2.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*(T - tau1)) - 4.*epsilon.*k_1.^2.*k_2.^4.*exp(k_2.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(-k_2.*(T - tau1)) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_2.*(T - tau1)) + 2.*epsilon.*k_1.^4.*k_2.^2.*exp(k_2.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(-k_2.*(T - tau1)) + 2.*epsilon.*k_1.^4.*k_2.^2.*exp(-k_2.*(T - tau1)) - 4.*k_1.^2.*k_2.^5.*tau1.*exp(k_2.*(T - tau1)) + 6.*k_1.^3.*k_2.^4.*tau1.*exp(k_2.*(T - tau1)) - 2.*k_1.^5.*k_2.^2.*tau1.*exp(k_2.*(T - tau1)) - 2.*k_1.^3.*k_2.^4.*tau1.*exp(-k_2.*(T - tau1)) + 4.*k_1.^4.*k_2.^3.*tau1.*exp(-k_2.*(T - tau1)) - 2.*k_1.^5.*k_2.^2.*tau1.*exp(-k_2.*(T - tau1)) - alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 3.*alpha.*k_1.^3.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) + 2.*alpha.*k_1.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*alpha.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*S.*k_1.*k_2.^4.*exp(k_1.*tau1) + 2.*S.*k_1.*k_2.^4.*exp(-k_1.*tau1) + alpha.*k_1.*k_2.^3.*exp(k_1.*tau1) + 4.*alpha.*k_1.^3.*k_2.*exp(k_1.*tau1) - alpha.*k_1.*k_2.^3.*exp(-k_1.*tau1) + 4.*alpha.*k_1.^3.*k_2.*exp(-k_1.*tau1) + 4.*S.*T.*k_1.^4.*k_2.^2 - 2.*epsilon.*k_1.*k_2.^4.*exp(k_1.*tau1) + 2.*epsilon.*k_1.*k_2.^4.*exp(-k_1.*tau1) - 2.*S.*k_1.^4.*k_2.*exp(k_2.*(T - tau1)) + 2.*S.*k_1.^4.*k_2.*exp(-k_2.*(T - tau1)) + 2.*L.*k_1.^2.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*L.*k_1.^2.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*S.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^2.*k_2.^5.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^2.*k_2.^5.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^3.*k_2.^4.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*T.*epsilon.*k_1.^4.*k_2.^2 - 2.*alpha.*k_1.*k_2.^3.*exp(k_2.*(T - tau1)) + 5.*alpha.*k_1.^3.*k_2.*exp(k_2.*(T - tau1)) + 2.*alpha.*k_1.*k_2.^3.*exp(-k_2.*(T - tau1)) + 3.*alpha.*k_1.^3.*k_2.*exp(-k_2.*(T - tau1)) - 2.*epsilon.*k_1.^4.*k_2.*exp(k_2.*(T - tau1)) + 2.*epsilon.*k_1.^4.*k_2.*exp(-k_2.*(T - tau1)) + 3.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) + 2.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*S.*k_1.^2.*k_2.^4.*tau1 - 4.*S.*k_1.^4.*k_2.^2.*tau1 - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*epsilon.*k_1.^2.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*alpha.*k_1.^2.*k_2.^3.*tau1 + 2.*alpha.*k_1.^3.*k_2.^2.*tau1 + 2.*k_1.^2.*k_2.^5.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^5.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 4.*k_1.^3.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^4.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*k_1.^4.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*epsilon.*k_1.^2.*k_2.^4.*tau1 - 4.*epsilon.*k_1.^4.*k_2.^2.*tau1 + 2.*S.*k_1.^2.*k_2.^3.*exp(k_1.*tau1) + 2.*S.*k_1.^2.*k_2.^3.*exp(-k_1.*tau1) - 4.*S.*k_1.^3.*k_2.^2.*exp(-k_1.*tau1) - 4.*L.*k_1.^2.*k_2.^4.*exp(k_2.*(T - tau1)) + 2.*L.*k_1.^3.*k_2.^3.*exp(k_2.*(T - tau1)) + 2.*L.*k_1.^4.*k_2.^2.*exp(k_2.*(T - tau1)) - 2.*L.*k_1.^3.*k_2.^3.*exp(-k_2.*(T - tau1)) + 2.*L.*k_1.^4.*k_2.^2.*exp(-k_2.*(T - tau1)) - 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) - 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(-k_1.*tau1) - 4.*epsilon.*k_1.^3.*k_2.^2.*exp(-k_1.*tau1) + 4.*S.*k_1.^2.*k_2.^3.*exp(k_2.*(T - tau1)) - 2.*S.*k_1.^3.*k_2.^2.*exp(k_2.*(T - tau1)) - 2.*S.*k_1.^3.*k_2.^2.*exp(-k_2.*(T - tau1)) + T.*alpha.*k_1.^4.*k_2.*exp(k_2.*(T - tau1)) - T.*alpha.*k_1.^4.*k_2.*exp(-k_2.*(T - tau1)) - T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - T.*alpha.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - alpha.*k_1.^4.*k_2.*tau1.*exp(k_2.*(T - tau1)) + alpha.*k_1.^4.*k_2.*tau1.*exp(-k_2.*(T - tau1)) + alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*T.*alpha.*k_1.^2.*k_2.^3.*exp(k_2.*(T - tau1)) - 3.*T.*alpha.*k_1.^3.*k_2.^2.*exp(k_2.*(T - tau1)) + T.*alpha.*k_1.^3.*k_2.^2.*exp(-k_2.*(T - tau1)) - alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*(T - tau1)) + 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_2.*(T - tau1)) + alpha.*k_1.^2.*k_2.^3.*tau1.*exp(-k_2.*(T - tau1)) - 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(-k_2.*(T - tau1)) + alpha.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - alpha.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1)))./(4.*k_1.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 4.*k_1.^4.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 8.*k_1.^4.*exp(-k_1.*tau1) - 8.*k_2.^4.*exp(k_2.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.^3.*k_2.*exp(k_1.*tau1) + 4.*k_1.^3.*k_2.*exp(-k_1.*tau1) + 4.*k_1.*k_2.^3.*exp(k_2.*(T - tau1)) - 4.*k_1.*k_2.^3.*exp(-k_2.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*(T - tau1)) - k_1.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 3.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^4.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1));
    
    psi_2_t_tau1_second = @(t, tau1) psi_2_t_1_0_tau1_second(t, psi_1_cur_tau1(tau1), psi_2_0_cur_tau1(tau1), tau1);
    x_2_t_tau1_second = @(t, tau1) x_2_t_1_0_tau1_second(t, psi_1_cur_tau1(tau1), psi_2_0_cur_tau1(tau1), tau1);
    
    tau1_cur = fsolve(@(tau1) psi_1_cur_tau1(tau1) + psi_2_t_tau1_second(tau1, tau1) .* x_2_t_tau1_second(tau1, tau1), T ./ 2, fsolve_opts);

    if (0 < tau1_cur && tau1_cur < T)
        psi_1_cur = psi_1_cur_tau1(tau1_cur);
        psi_2_0_cur = psi_2_0_cur_tau1(tau1_cur);
        
        psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);

        if (psi_1_cur < 0 && ...
            psi_2_t_second(T) < 0)

            u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
            u_1_t_second = @(t) u_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);
            
            psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);
            x_2_t_second = @(t) x_2_t_tau1_second(t, tau1_cur);
            
            index = find(t_opt - tau1_cur >= 0, 1);
            
            if (prod(psi_1_cur + psi_2_t_second(t_opt(index : end)) .* x_2_t_second(t_opt(index : end)) >= 0))

                functional_cur = integral(u_1_t_first, 0, tau1_cur) + integral(u_1_t_second, tau1_cur, T);

                if (functional_cur < functional_min)

                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    index = find(t_opt - tau1_cur >= 0, 1);

                    psi_2_t_first =  @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
                    x_1_t_second = @(t) x_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);

                    x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_final .* ones(1, number_of_points_for_splitting);
                    psi2_opt = [psi_2_t_first(t_opt(1 : (index - 1))), psi_2_t_second(t_opt(index : end))];

                    x1_opt = [x_1_t_first(t_opt(1 : (index - 1))), x_1_t_second(t_opt(index : end))];
                    x2_opt = [x_2_t_first(t_opt(1 : (index - 1))), x_2_t_second(t_opt(index : end))];

                    u1_opt = [u_1_t_first(t_opt(1 : (index - 1))), u_1_t_second(t_opt(index : end))];
                    u2_opt = [u2_first .* ones(1, index - 1), u2_second .* ones(1, number_of_points_for_splitting - index + 1)];

                    functional_min = functional_cur;

                    switches = [tau1_cur];

                end
            end
        end
    end
    
    % psi_2(T) > 0
    
    psi_1_cur_tau1 = @(tau1) -(2.*k_1.*k_2.*(4.*S.*k_1.*k_2.^3 - 4.*alpha.*k_2.^3.*exp(k_2.*(T - tau1)) - 4.*alpha.*k_1.^3.*exp(-k_1.*tau1) + 2.*alpha.*k_1.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*alpha.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.*k_2.^2 - 2.*alpha.*k_1.^2.*k_2 - 4.*epsilon.*k_1.*k_2.^3 - 2.*S.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - alpha.*k_1.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*alpha.*k_1.^2.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*epsilon.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*S.*k_1.*k_2.^3.*exp(k_1.*tau1) - 2.*S.*k_1.*k_2.^3.*exp(-k_1.*tau1) + 4.*S.*k_1.^3.*k_2.*exp(-k_1.*tau1) + alpha.*k_1.*k_2.^2.*exp(k_1.*tau1) - alpha.*k_1.^2.*k_2.*exp(k_1.*tau1) + alpha.*k_1.*k_2.^2.*exp(-k_1.*tau1) + 3.*alpha.*k_1.^2.*k_2.*exp(-k_1.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(k_1.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(-k_1.*tau1) - 4.*epsilon.*k_1.^3.*k_2.*exp(-k_1.*tau1) + 2.*L.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^2.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^2.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^3.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 3.*alpha.*k_1.*k_2.^2.*exp(k_2.*(T - tau1)) + alpha.*k_1.^2.*k_2.*exp(k_2.*(T - tau1)) - alpha.*k_1.*k_2.^2.*exp(-k_2.*(T - tau1)) + alpha.*k_1.^2.*k_2.*exp(-k_2.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*k_1.^2.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 4.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.^4.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^4.*k_2.^2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*S.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) - 2.*S.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) - T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.^2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - T.*alpha.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - alpha.*k_1.*k_2.^3.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^3.*k_2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2)))./(4.*k_1.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 4.*k_1.^4.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 8.*k_1.^4.*exp(-k_1.*tau1) - 8.*k_2.^4.*exp(k_2.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.^3.*k_2.*exp(k_1.*tau1) + 4.*k_1.^3.*k_2.*exp(-k_1.*tau1) + 4.*k_1.*k_2.^3.*exp(k_2.*(T - tau1)) - 4.*k_1.*k_2.^3.*exp(-k_2.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*(T - tau1)) - k_1.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 3.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^4.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1));
    psi_2_0_cur_tau1 = @(tau1) (2.*(4.*L.*k_1.^4.*k_2.^2 - 4.*alpha.*k_1.^4 + 2.*alpha.*k_1.^4.*exp(k_2.*(T - tau1)) + 2.*alpha.*k_1.^4.*exp(-k_2.*(T - tau1)) - 4.*alpha.*k_2.^4.*exp(k_2.*(T - tau1)) + 4.*S.*k_1.^2.*k_2.^3 - 4.*S.*k_1.^3.*k_2.^2 - 4.*T.*k_1.^4.*k_2.^3 - 8.*alpha.*k_1.^2.*k_2.^2 - 4.*epsilon.*k_1.^2.*k_2.^3 + 4.*epsilon.*k_1.^3.*k_2.^2 + 4.*epsilon.*k_1.^4.*k_2.^2 + 4.*k_1.^4.*k_2.^3.*tau1 - 4.*k_1.^5.*k_2.^2.*tau1 + 2.*alpha.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 8.*alpha.*k_1.^3.*k_2 - 4.*T.*k_1.^2.*k_2.^5.*exp(k_2.*(T - tau1)) + 2.*T.*k_1.^3.*k_2.^4.*exp(k_2.*(T - tau1)) + 2.*T.*k_1.^4.*k_2.^3.*exp(k_2.*(T - tau1)) - 2.*T.*k_1.^3.*k_2.^4.*exp(-k_2.*(T - tau1)) + 2.*T.*k_1.^4.*k_2.^3.*exp(-k_2.*(T - tau1)) + 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_2.*(T - tau1)) + 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_2.*(T - tau1)) + 4.*epsilon.*k_1.^2.*k_2.^3.*exp(k_2.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*(T - tau1)) + 4.*epsilon.*k_1.^2.*k_2.^4.*exp(k_2.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(-k_2.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_2.*(T - tau1)) - 2.*epsilon.*k_1.^4.*k_2.^2.*exp(k_2.*(T - tau1)) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(-k_2.*(T - tau1)) - 2.*epsilon.*k_1.^4.*k_2.^2.*exp(-k_2.*(T - tau1)) + 4.*k_1.^2.*k_2.^5.*tau1.*exp(k_2.*(T - tau1)) - 6.*k_1.^3.*k_2.^4.*tau1.*exp(k_2.*(T - tau1)) + 2.*k_1.^5.*k_2.^2.*tau1.*exp(k_2.*(T - tau1)) + 2.*k_1.^3.*k_2.^4.*tau1.*exp(-k_2.*(T - tau1)) - 4.*k_1.^4.*k_2.^3.*tau1.*exp(-k_2.*(T - tau1)) + 2.*k_1.^5.*k_2.^2.*tau1.*exp(-k_2.*(T - tau1)) + alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) + 3.*alpha.*k_1.^3.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 2.*alpha.*k_1.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*alpha.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*S.*k_1.*k_2.^4.*exp(k_1.*tau1) - 2.*S.*k_1.*k_2.^4.*exp(-k_1.*tau1) - alpha.*k_1.*k_2.^3.*exp(k_1.*tau1) - 4.*alpha.*k_1.^3.*k_2.*exp(k_1.*tau1) + alpha.*k_1.*k_2.^3.*exp(-k_1.*tau1) - 4.*alpha.*k_1.^3.*k_2.*exp(-k_1.*tau1) - 4.*S.*T.*k_1.^4.*k_2.^2 - 2.*epsilon.*k_1.*k_2.^4.*exp(k_1.*tau1) + 2.*epsilon.*k_1.*k_2.^4.*exp(-k_1.*tau1) + 2.*S.*k_1.^4.*k_2.*exp(k_2.*(T - tau1)) - 2.*S.*k_1.^4.*k_2.*exp(-k_2.*(T - tau1)) - 2.*L.*k_1.^2.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^2.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*S.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*T.*k_1.^2.*k_2.^5.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^2.*k_2.^5.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^3.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^4.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*T.*epsilon.*k_1.^4.*k_2.^2 + 2.*alpha.*k_1.*k_2.^3.*exp(k_2.*(T - tau1)) - 5.*alpha.*k_1.^3.*k_2.*exp(k_2.*(T - tau1)) - 2.*alpha.*k_1.*k_2.^3.*exp(-k_2.*(T - tau1)) - 3.*alpha.*k_1.^3.*k_2.*exp(-k_2.*(T - tau1)) - 2.*epsilon.*k_1.^4.*k_2.*exp(k_2.*(T - tau1)) + 2.*epsilon.*k_1.^4.*k_2.*exp(-k_2.*(T - tau1)) - 3.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 4.*S.*k_1.^2.*k_2.^4.*tau1 + 4.*S.*k_1.^4.*k_2.^2.*tau1 - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*epsilon.*k_1.^2.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 2.*alpha.*k_1.^2.*k_2.^3.*tau1 - 2.*alpha.*k_1.^3.*k_2.^2.*tau1 - 2.*k_1.^2.*k_2.^5.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^5.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 4.*k_1.^3.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^4.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.^4.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^4.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*epsilon.*k_1.^2.*k_2.^4.*tau1 - 4.*epsilon.*k_1.^4.*k_2.^2.*tau1 - 2.*S.*k_1.^2.*k_2.^3.*exp(k_1.*tau1) - 2.*S.*k_1.^2.*k_2.^3.*exp(-k_1.*tau1) + 4.*S.*k_1.^3.*k_2.^2.*exp(-k_1.*tau1) + 4.*L.*k_1.^2.*k_2.^4.*exp(k_2.*(T - tau1)) - 2.*L.*k_1.^3.*k_2.^3.*exp(k_2.*(T - tau1)) - 2.*L.*k_1.^4.*k_2.^2.*exp(k_2.*(T - tau1)) + 2.*L.*k_1.^3.*k_2.^3.*exp(-k_2.*(T - tau1)) - 2.*L.*k_1.^4.*k_2.^2.*exp(-k_2.*(T - tau1)) + 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(-k_1.*tau1) - 4.*epsilon.*k_1.^3.*k_2.^2.*exp(-k_1.*tau1) - 4.*S.*k_1.^2.*k_2.^3.*exp(k_2.*(T - tau1)) + 2.*S.*k_1.^3.*k_2.^2.*exp(k_2.*(T - tau1)) + 2.*S.*k_1.^3.*k_2.^2.*exp(-k_2.*(T - tau1)) - T.*alpha.*k_1.^4.*k_2.*exp(k_2.*(T - tau1)) + T.*alpha.*k_1.^4.*k_2.*exp(-k_2.*(T - tau1)) + T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.^3.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + alpha.*k_1.^4.*k_2.*tau1.*exp(k_2.*(T - tau1)) - alpha.*k_1.^4.*k_2.*tau1.*exp(-k_2.*(T - tau1)) - alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*alpha.*k_1.^2.*k_2.^3.*exp(k_2.*(T - tau1)) + 3.*T.*alpha.*k_1.^3.*k_2.^2.*exp(k_2.*(T - tau1)) - T.*alpha.*k_1.^3.*k_2.^2.*exp(-k_2.*(T - tau1)) + alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*(T - tau1)) - 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_2.*(T - tau1)) - alpha.*k_1.^2.*k_2.^3.*tau1.*exp(-k_2.*(T - tau1)) + 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(-k_2.*(T - tau1)) - alpha.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1)))./(4.*k_1.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 4.*k_1.^4.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_2.^4.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 8.*k_1.^4.*exp(-k_1.*tau1) - 8.*k_2.^4.*exp(k_2.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.^3.*k_2.*exp(k_1.*tau1) + 4.*k_1.^3.*k_2.*exp(-k_1.*tau1) + 4.*k_1.*k_2.^3.*exp(k_2.*(T - tau1)) - 4.*k_1.*k_2.^3.*exp(-k_2.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*(T - tau1)) - k_1.*k_2.^3.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_2 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.*k_2.^3.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 3.*k_1.^3.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - k_1.^3.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^4.*k_2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^2.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 + k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(k_2.*tau1 - k_1.*tau1 - T.*k_2) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_2 - k_1.*tau1 - k_2.*tau1));
    
    tau1_cur = fsolve(@(tau1) psi_1_cur_tau1(tau1) + psi_2_t_tau1_second(tau1, tau1) .* x_2_t_tau1_second(tau1, tau1), 0, fsolve_opts);

    if (0 < tau1_cur && tau1_cur < T)
        psi_1_cur = psi_1_cur_tau1(tau1_cur);
        psi_2_0_cur = psi_2_0_cur_tau1(tau1_cur);
        
        psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);

        if (psi_1_cur < 0 && ...
            psi_2_t_second(T) > 0)

            u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
            u_1_t_second = @(t) u_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);
            
            psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);
            x_2_t_second = @(t) x_2_t_tau1_second(t, tau1_cur);
            
            index = find(t_opt - tau1_cur >= 0, 1);
            
            if (prod(psi_1_cur + psi_2_t_second(t_opt(index : end)) .* x_2_t_second(t_opt(index : end)) >= 0))

                functional_cur = integral(u_1_t_first, 0, tau1_cur) + integral(u_1_t_second, tau1_cur, T);

                if (functional_cur < functional_min)

                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    index = find(t_opt - tau1_cur >= 0, 1);

                    psi_2_t_first =  @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
                    x_1_t_second = @(t) x_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);

                    x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_cur .* ones(1, number_of_points_for_splitting);
                    psi2_opt = [psi_2_t_first(t_opt(1 : (index - 1))), psi_2_t_second(t_opt(index : end))];

                    x1_opt = [x_1_t_first(t_opt(1 : (index - 1))), x_1_t_second(t_opt(index : end))];
                    x2_opt = [x_2_t_first(t_opt(1 : (index - 1))), x_2_t_second(t_opt(index : end))];

                    u1_opt = [u_1_t_first(t_opt(1 : (index - 1))), u_1_t_second(t_opt(index : end))];
                    u2_opt = [u2_first .* ones(1, index - 1), u2_second .* ones(1, number_of_points_for_splitting - index + 1)];

                    functional_min = functional_cur;

                    switches = [tau1_cur];

                end
            end
        end
    end
    
    % two or more switches
    
    psi_1_cur_1_tau1_tau2_0 = @(tau1, tau2, psi_2_0) (k_1.*(2.*psi_2_0 + exp(k_1.*tau1).*(4.*alpha.^2 - 16.*k_1.^2.*psi_2_0 - 4.*alpha.*psi_2_0 - 16.*alpha.^2.*exp(k_1.*tau1) + 24.*alpha.^2.*exp(2.*k_1.*tau1) - 16.*alpha.^2.*exp(3.*k_1.*tau1) + 4.*alpha.^2.*exp(4.*k_1.*tau1) + 16.*k_1.^4.*exp(2.*k_1.*tau1) + psi_2_0.^2 - 4.*psi_2_0.^2.*exp(k_1.*tau1) + 6.*psi_2_0.^2.*exp(2.*k_1.*tau1) - 4.*psi_2_0.^2.*exp(3.*k_1.*tau1) + psi_2_0.^2.*exp(4.*k_1.*tau1) + 16.*alpha.*k_1.^2.*exp(k_1.*tau1) - 32.*alpha.*k_1.^2.*exp(2.*k_1.*tau1) + 16.*alpha.*k_1.^2.*exp(3.*k_1.*tau1) + 24.*k_1.^2.*psi_2_0.*exp(k_1.*tau1) - 8.*k_1.^2.*psi_2_0.*exp(3.*k_1.*tau1) + 16.*alpha.*psi_2_0.*exp(k_1.*tau1) - 24.*alpha.*psi_2_0.*exp(2.*k_1.*tau1) + 16.*alpha.*psi_2_0.*exp(3.*k_1.*tau1) - 4.*alpha.*psi_2_0.*exp(4.*k_1.*tau1)).^(1./2) - 2.*alpha.*exp(k_1.*tau1) + 4.*alpha.*exp(2.*k_1.*tau1) - 2.*alpha.*exp(3.*k_1.*tau1) - 3.*psi_2_0.*exp(k_1.*tau1) + psi_2_0.*exp(3.*k_1.*tau1) - 4.*k_1.^2.*exp(2.*k_1.*tau1)))./(2.*(exp(k_1.*tau1) - 1).^3);
    psi_1_cur_2_tau1_tau2_0 = @(tau1, tau2, psi_2_0) -(k_1.*(exp(k_1.*tau1).*(4.*alpha.^2 - 16.*k_1.^2.*psi_2_0 - 4.*alpha.*psi_2_0 - 16.*alpha.^2.*exp(k_1.*tau1) + 24.*alpha.^2.*exp(2.*k_1.*tau1) - 16.*alpha.^2.*exp(3.*k_1.*tau1) + 4.*alpha.^2.*exp(4.*k_1.*tau1) + 16.*k_1.^4.*exp(2.*k_1.*tau1) + psi_2_0.^2 - 4.*psi_2_0.^2.*exp(k_1.*tau1) + 6.*psi_2_0.^2.*exp(2.*k_1.*tau1) - 4.*psi_2_0.^2.*exp(3.*k_1.*tau1) + psi_2_0.^2.*exp(4.*k_1.*tau1) + 16.*alpha.*k_1.^2.*exp(k_1.*tau1) - 32.*alpha.*k_1.^2.*exp(2.*k_1.*tau1) + 16.*alpha.*k_1.^2.*exp(3.*k_1.*tau1) + 24.*k_1.^2.*psi_2_0.*exp(k_1.*tau1) - 8.*k_1.^2.*psi_2_0.*exp(3.*k_1.*tau1) + 16.*alpha.*psi_2_0.*exp(k_1.*tau1) - 24.*alpha.*psi_2_0.*exp(2.*k_1.*tau1) + 16.*alpha.*psi_2_0.*exp(3.*k_1.*tau1) - 4.*alpha.*psi_2_0.*exp(4.*k_1.*tau1)).^(1./2) - 2.*psi_2_0 + 2.*alpha.*exp(k_1.*tau1) - 4.*alpha.*exp(2.*k_1.*tau1) + 2.*alpha.*exp(3.*k_1.*tau1) + 3.*psi_2_0.*exp(k_1.*tau1) - psi_2_0.*exp(3.*k_1.*tau1) + 4.*k_1.^2.*exp(2.*k_1.*tau1)))./(2.*(exp(k_1.*tau1) - 1).^3);

    solving_equat = @(tau1, tau2, psi_1, psi_2_0) psi_1 - (exp(k_2.*(tau1 - tau2)).*(psi_2_0.*exp(-k_1.*tau1) + (psi_1.*(exp(-k_1.*tau1) - 1))./k_1) + (psi_1.*(exp(k_2.*(tau1 - tau2)) - 1))./k_2).*((alpha.*(exp(-k_2.*(tau1 - tau2)) - 1))./(2.*k_2) + (psi_1.*(exp(-k_2.*(tau1 - tau2)) - 1))./(2.*k_2.^2) + (psi_1.*exp(- k_2.*tau1 - k_2.*tau2).*(exp(2.*k_2.*tau1) - exp(2.*k_2.*tau2)))./(4.*k_2.^2) + (psi_2_0.*exp(- k_2.*tau1 - k_2.*tau2).*exp(-k_1.*tau1).*(exp(2.*k_2.*tau1) - exp(2.*k_2.*tau2)))./(4.*k_2) - (psi_1.*exp(- k_2.*tau1 - k_2.*tau2).*(exp(2.*k_2.*tau1) - exp(2.*k_2.*tau2)))./(4.*k_1.*k_2) - (exp(-k_1.*tau1).*exp(-k_2.*(tau1 - tau2)).*(exp(k_1.*tau1) - 1).*(psi_1 + k_1.*psi_2_0 - psi_1.*exp(k_1.*tau1) - 2.*alpha.*k_1.*exp(k_1.*tau1) + k_1.*psi_2_0.*exp(k_1.*tau1)))./(4.*k_1.^2) + (psi_1.*exp(- k_2.*tau1 - k_2.*tau2).*exp(-k_1.*tau1).*(exp(2.*k_2.*tau1) - exp(2.*k_2.*tau2)))./(4.*k_1.*k_2));
    
    psi_2_0_cur_11_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_1_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), alpha, fsolve_opts);
    
    psi_2_0_cur_12_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_1_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), -alpha, fsolve_opts);
    
    psi_2_0_cur_21_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_2_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), alpha, fsolve_opts);
    
    psi_2_0_cur_22_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_2_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), -alpha, fsolve_opts);

    tau1_splitting = linspace(0 + trimmer, T - trimmer, number_of_points_for_enum);
        
    counter = 0;
    for tau1_cur = tau1_splitting
        counter = counter + 1;
        tau2_splitting = tau1_splitting((counter + 1) : end);
        for tau2_cur = tau2_splitting

            psi_2_0_cur_11 = psi_2_0_cur_11_tau1_tau2(tau1_cur, tau2_cur);
            psi_2_0_cur_12 = psi_2_0_cur_12_tau1_tau2(tau1_cur, tau2_cur);
            psi_2_0_cur_21 = psi_2_0_cur_21_tau1_tau2(tau1_cur, tau2_cur);
            psi_2_0_cur_22 = psi_2_0_cur_22_tau1_tau2(tau1_cur, tau2_cur);
            
            psi_1_cur_11 = psi_1_cur_1_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_11);
            psi_1_cur_12 = psi_1_cur_1_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_12);
            psi_1_cur_21 = psi_1_cur_2_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_21);
            psi_1_cur_22 = psi_1_cur_2_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_22);

            if (imag(psi_1_cur_11) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_1, k_2, alpha, psi_1_cur_11), [0, T], [0, 0, psi_2_0_cur_11]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_11 < 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_11 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_11 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_11 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;

                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);
                    end
                end
            end    

            if (imag(psi_1_cur_12) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_1, k_2, alpha, psi_1_cur_12), [0, T], [0, 0, psi_2_0_cur_12]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_12 < 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_12 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_12 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_12 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;

                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);
                    end
                end
            end    
            
            if (imag(psi_1_cur_21) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_1, k_2, alpha, psi_1_cur_11), [0, T], [0, 0, psi_2_0_cur_21]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_21 < 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_21 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_21 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_21 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;

                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);
                    end
                end
            end    
            
            if (imag(psi_1_cur_22) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_1, k_2, alpha, psi_1_cur_22), [0, T], [0, 0, psi_2_0_cur_22]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_22 < 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_22 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_22 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_22 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;

                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);
                    end
                end
            end    
        end
    end
    
    % sixth case: psi_1 > 0 (the last!)
    
    u2_first = k_2;
    
    psi_2_t_1_0_first = @(t, psi_1, psi_2_0) psi_2_0.*exp(-k_2.*t) + (psi_1.*(exp(-k_2.*t) - 1))./k_2;
    u_1_t_1_0_first = @(t, psi_1, psi_2_0) (psi_2_t_1_0_first(t, psi_1, psi_2_0) - alpha) ./ 2;
    
    x_1_t_1_0_first = @(t, psi_1, psi_2_0) (2.*alpha.*k_2 - 2.*k_2.*psi_2_0 - psi_1.*exp(k_2.*t) + psi_1.*exp(-k_2.*t) + 4.*k_2.^4.*t + 2.*k_2.*psi_1.*t - 2.*alpha.*k_2.*exp(k_2.*t) + k_2.*psi_2_0.*exp(k_2.*t) + k_2.*psi_2_0.*exp(-k_2.*t) + 2.*alpha.*k_2.^2.*t)./(4.*k_2.^3);
    x_2_t_1_0_first = @(t, psi_1, psi_2_0) (exp(-k_2.*t).*(exp(k_2.*t) - 1).*(psi_1 + k_2.*psi_2_0 - psi_1.*exp(k_2.*t) - 2.*alpha.*k_2.*exp(k_2.*t) + k_2.*psi_2_0.*exp(k_2.*t)))./(4.*k_2.^2);
    
    % without switches
    
    % psi_2(T) < 0
    
    psi_1_cur = -(k_2.*(2.*alpha - 2.*S.*k_2 - 2.*epsilon.*k_2 - 2.*alpha.*exp(T.*k_2) - 2.*L.*k_2.^2 + 2.*T.*k_2.^3 + 2.*epsilon.*k_2.^2 + 2.*T.*k_2.^3.*exp(T.*k_2) + 2.*epsilon.*k_2.^2.*exp(T.*k_2) + T.*alpha.*k_2 + 2.*S.*k_2.*exp(T.*k_2) + 2.*epsilon.*k_2.*exp(T.*k_2) - 2.*L.*k_2.^2.*exp(T.*k_2) + T.*alpha.*k_2.*exp(T.*k_2)))./(T.*k_2 - 2.*exp(T.*k_2) + T.*k_2.*exp(T.*k_2) + 2);
    psi_2_0_cur = -(2.*alpha - 2.*S.*k_2 - 2.*epsilon.*k_2 - 4.*alpha.*exp(T.*k_2) + 2.*alpha.*exp(2.*T.*k_2) - 2.*L.*k_2.^2 + 2.*T.*k_2.^3 + 2.*epsilon.*k_2.^2 - 4.*T.*k_2.^3.*exp(T.*k_2) + 2.*T.*k_2.^3.*exp(2.*T.*k_2) - 4.*epsilon.*k_2.^2.*exp(T.*k_2) + 2.*epsilon.*k_2.^2.*exp(2.*T.*k_2) + T.*alpha.*k_2 + 2.*S.*k_2.*exp(2.*T.*k_2) + 2.*epsilon.*k_2.*exp(2.*T.*k_2) + 4.*L.*k_2.^2.*exp(T.*k_2) - 2.*L.*k_2.^2.*exp(2.*T.*k_2) - T.*alpha.*k_2.*exp(2.*T.*k_2) - 4.*S.*T.*k_2.^2.*exp(T.*k_2) - 4.*T.*epsilon.*k_2.^2.*exp(T.*k_2))./((exp(T.*k_2) - 1).*(T.*k_2 - 2.*exp(T.*k_2) + T.*k_2.*exp(T.*k_2) + 2));
    
    if (psi_1_cur > 0)
        u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        psi_2_t_first = @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        
        if (psi_2_t_first(T) < 0)
        
            x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

            if (prod(psi_1_cur + psi_2_t_first(t_opt) .* x_2_t_first(t_opt) >= 0))
                int_func = @(t) u_1_t_first(t) .^ 2 + alpha .* u_1_t_first(t);

                functional_cur = integral(int_func, 0, T);

                if (functional_cur < functional_min)
                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_cur .* ones(1, number_of_points_for_splitting);
                    psi2_opt = psi_2_t_first(t_opt);

                    x1_opt = x_1_t_first(t_opt);
                    x2_opt = x_2_t_first(t_opt);

                    u1_opt = u_1_t_first(t_opt);
                    u2_opt = u2_first .* ones(1, number_of_points_for_splitting);

                    functional_min = functional_cur;
                end
            end
        end
    end
    
    % psi_2(T) > 0
    
    psi_1_cur = -(k_2.*(2.*alpha - 2.*S.*k_2 + 2.*epsilon.*k_2 - 2.*alpha.*exp(T.*k_2) - 2.*L.*k_2.^2 + 2.*T.*k_2.^3 + 2.*epsilon.*k_2.^2 + 2.*T.*k_2.^3.*exp(T.*k_2) + 2.*epsilon.*k_2.^2.*exp(T.*k_2) + T.*alpha.*k_2 + 2.*S.*k_2.*exp(T.*k_2) - 2.*epsilon.*k_2.*exp(T.*k_2) - 2.*L.*k_2.^2.*exp(T.*k_2) + T.*alpha.*k_2.*exp(T.*k_2)))./(T.*k_2 - 2.*exp(T.*k_2) + T.*k_2.*exp(T.*k_2) + 2);
    psi_2_0_cur = -(2.*alpha - 2.*S.*k_2 + 2.*epsilon.*k_2 - 4.*alpha.*exp(T.*k_2) + 2.*alpha.*exp(2.*T.*k_2) - 2.*L.*k_2.^2 + 2.*T.*k_2.^3 + 2.*epsilon.*k_2.^2 - 4.*T.*k_2.^3.*exp(T.*k_2) + 2.*T.*k_2.^3.*exp(2.*T.*k_2) - 4.*epsilon.*k_2.^2.*exp(T.*k_2) + 2.*epsilon.*k_2.^2.*exp(2.*T.*k_2) + T.*alpha.*k_2 + 2.*S.*k_2.*exp(2.*T.*k_2) - 2.*epsilon.*k_2.*exp(2.*T.*k_2) + 4.*L.*k_2.^2.*exp(T.*k_2) - 2.*L.*k_2.^2.*exp(2.*T.*k_2) - T.*alpha.*k_2.*exp(2.*T.*k_2) - 4.*S.*T.*k_2.^2.*exp(T.*k_2) + 4.*T.*epsilon.*k_2.^2.*exp(T.*k_2))./((exp(T.*k_2) - 1).*(T.*k_2 - 2.*exp(T.*k_2) + T.*k_2.*exp(T.*k_2) + 2));
    
    if (psi_1_cur > 0)
        u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        psi_2_t_first = @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
        
        if (psi_2_t_first(T) > 0)
        
            x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

            if (prod(psi_1_cur + psi_2_t_first(t_opt) .* x_2_t_first(t_opt) >= 0))
                int_func = @(t) u_1_t_first(t) .^ 2 + alpha .* u_1_t_first(t);

                functional_cur = integral(int_func, 0, T);

                if (functional_cur < functional_min)
                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_cur .* ones(1, number_of_points_for_splitting);
                    psi2_opt = psi_2_t_first(t_opt);

                    x1_opt = x_1_t_first(t_opt);
                    x2_opt = x_2_t_first(t_opt);

                    u1_opt = u_1_t_first(t_opt);
                    u2_opt = u2_first .* ones(1, number_of_points_for_splitting);

                    functional_min = functional_cur;
                end
            end
        end
    end
    
    % one switch
    
    u2_second = k_1;
    
    psi_2_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) exp(-k_1.*(t - tau1)).*(psi_2_0.*exp(-k_2.*tau1) + (psi_1.*(exp(-k_2.*tau1) - 1))./k_2) + (psi_1.*(exp(-k_1.*(t - tau1)) - 1))./k_1;
    
    u_1_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) (psi_2_t_1_0_tau1_second(t, psi_1, psi_2_0, tau1) - alpha) ./ 2;
    
    x_2_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) (psi_1.*exp(- k_1.*t - k_1.*tau1).*(exp(2.*k_1.*t) - exp(2.*k_1.*tau1)))./(4.*k_1.^2) - (psi_1.*(exp(k_1.*(t - tau1)) - 1))./(2.*k_1.^2) - (alpha.*(exp(k_1.*(t - tau1)) - 1))./(2.*k_1) + (psi_2_0.*exp(- k_1.*t - k_1.*tau1).*exp(-k_2.*tau1).*(exp(2.*k_1.*t) - exp(2.*k_1.*tau1)))./(4.*k_1) - (psi_1.*exp(- k_1.*t - k_1.*tau1).*(exp(2.*k_1.*t) - exp(2.*k_1.*tau1)))./(4.*k_1.*k_2) + (exp(-k_2.*tau1).*exp(k_1.*(t - tau1)).*(exp(k_2.*tau1) - 1).*(psi_1 + k_2.*psi_2_0 - psi_1.*exp(k_2.*tau1) - 2.*alpha.*k_2.*exp(k_2.*tau1) + k_2.*psi_2_0.*exp(k_2.*tau1)))./(4.*k_2.^2) + (psi_1.*exp(- k_1.*t - k_1.*tau1).*exp(-k_2.*tau1).*(exp(2.*k_1.*t) - exp(2.*k_1.*tau1)))./(4.*k_1.*k_2);
    
    x_1_t_1_0_tau1_second = @(t, psi_1, psi_2_0, tau1) k_2.*tau1 + k_1.*(t - tau1) + (alpha.*tau1)./(2.*k_2) + (psi_1.*tau1)./(2.*k_2.^2) + (alpha.*(t - tau1))./(2.*k_1) + (psi_1.*(t - tau1))./(2.*k_1.^2) + (psi_1 - psi_1.*exp(-k_1.*(t - tau1)))./(4.*k_1.^2.*k_2) - (alpha.*(exp(k_2.*tau1) - 1))./(2.*k_2.^2) - (psi_1.*(exp(k_2.*tau1) - 1))./(4.*k_2.^3) + (psi_1.*(exp(-k_2.*tau1) - 1))./(4.*k_2.^3) + (psi_2_0.*(exp(k_2.*tau1) - 1))./(4.*k_2.^2) + (psi_2_0.*(exp(-k_2.*tau1) - 1))./(4.*k_2.^2) - (alpha.*(exp(k_1.*(t - tau1)) - 1))./(2.*k_1.^2) - (psi_1.*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.^3) + (psi_1.*exp(2.*k_1.*tau1).*(exp(-k_1.*(t + tau1)) - exp(-2.*k_1.*tau1)))./(4.*k_1.^3) + (psi_2_0.*exp(-k_2.*tau1).*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.^2) + (alpha.*(exp(k_1.*(t - tau1)) - 1))./(2.*k_1.*k_2) + (psi_1.*(exp(k_1.*(t - tau1)) - 1))./(2.*k_1.*k_2.^2) - (psi_1.*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.^2.*k_2) + (psi_2_0.*exp(tau1.*(2.*k_1 - k_2)).*(exp(-k_1.*(t + tau1)) - exp(-2.*k_1.*tau1)))./(4.*k_1.^2) - (alpha.*exp(k_2.*tau1).*(exp(k_1.*(t - tau1)) - 1))./(2.*k_1.*k_2) - (psi_1.*exp(k_2.*tau1).*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.*k_2.^2) - (psi_1.*exp(-k_2.*tau1).*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.*k_2.^2) + (psi_1.*exp(-k_2.*tau1).*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.^2.*k_2) + (psi_2_0.*exp(k_2.*tau1).*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.*k_2) - (psi_2_0.*exp(-k_2.*tau1).*(exp(k_1.*(t - tau1)) - 1))./(4.*k_1.*k_2) + (psi_1.*exp(tau1.*(2.*k_1 - k_2)).*(exp(-k_1.*(t + tau1)) - exp(-2.*k_1.*tau1)))./(4.*k_1.^2.*k_2);
    
    % psi_2(T) < 0
    
    psi_1_cur_tau1 = @(tau1) (2.*k_1.*k_2.*(4.*alpha.*k_2.^3.*exp(-k_2.*tau1) + 4.*alpha.*k_1.^3.*exp(k_1.*(T - tau1)) - 4.*S.*k_1.^3.*k_2 - 2.*alpha.*k_1.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*alpha.*k_1.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*alpha.*k_1.*k_2.^2 + 2.*alpha.*k_1.^2.*k_2 - 4.*epsilon.*k_1.^3.*k_2 + 2.*S.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - alpha.*k_1.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.^2.*k_2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*alpha.*k_1.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*alpha.*k_1.^2.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.*exp(k_2.*tau1) - 4.*S.*k_1.*k_2.^3.*exp(-k_2.*tau1) + 2.*S.*k_1.^3.*k_2.*exp(-k_2.*tau1) + alpha.*k_1.*k_2.^2.*exp(k_2.*tau1) - alpha.*k_1.^2.*k_2.*exp(k_2.*tau1) - 3.*alpha.*k_1.*k_2.^2.*exp(-k_2.*tau1) - alpha.*k_1.^2.*k_2.*exp(-k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.*exp(k_2.*tau1) - 4.*epsilon.*k_1.*k_2.^3.*exp(-k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.*exp(-k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*L.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*L.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^3.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*T.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.*k_2.^2.*exp(k_1.*(T - tau1)) - 3.*alpha.*k_1.^2.*k_2.*exp(k_1.*(T - tau1)) - alpha.*k_1.*k_2.^2.*exp(-k_1.*(T - tau1)) + alpha.*k_1.^2.*k_2.*exp(-k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^4.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*k_1.^2.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 4.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*k_1.^4.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) + T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + T.*alpha.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - T.*alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - alpha.*k_1.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - alpha.*k_1.^3.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1)))./(4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 8.*k_2.^4.*exp(-k_2.*tau1) - 8.*k_1.^4.*exp(k_1.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.*k_2.^3.*exp(k_2.*tau1) + 4.*k_1.*k_2.^3.*exp(-k_2.*tau1) + 4.*k_1.^3.*k_2.*exp(k_1.*(T - tau1)) - 4.*k_1.^3.*k_2.*exp(-k_1.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*(T - tau1)) + 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 3.*k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*T.*k_1.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1));
    psi_2_0_cur_tau1 = @(tau1) (2.*(4.*L.*k_1.^2.*k_2.^4 - 4.*alpha.*k_2.^4 - 4.*alpha.*k_1.^4.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_2.^4.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_2.^4.*exp(-k_1.*(T - tau1)) - 4.*S.*k_1.^2.*k_2.^3 + 4.*S.*k_1.^3.*k_2.^2 - 4.*T.*k_1.^3.*k_2.^4 - 8.*alpha.*k_1.^2.*k_2.^2 - 4.*epsilon.*k_1.^2.*k_2.^3 + 4.*epsilon.*k_1.^3.*k_2.^2 - 4.*epsilon.*k_1.^2.*k_2.^4 - 4.*k_1.^2.*k_2.^5.*tau1 + 4.*k_1.^3.*k_2.^4.*tau1 + 2.*alpha.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*alpha.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 8.*alpha.*k_1.*k_2.^3 + 2.*T.*k_1.^3.*k_2.^4.*exp(k_1.*(T - tau1)) + 2.*T.*k_1.^4.*k_2.^3.*exp(k_1.*(T - tau1)) - 4.*T.*k_1.^5.*k_2.^2.*exp(k_1.*(T - tau1)) + 2.*T.*k_1.^3.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*T.*k_1.^4.*k_2.^3.*exp(-k_1.*(T - tau1)) + 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*(T - tau1)) + 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*(T - tau1)) - 4.*epsilon.*k_1.^3.*k_2.^2.*exp(k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(-k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^4.*exp(k_1.*(T - tau1)) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_1.*(T - tau1)) - 4.*epsilon.*k_1.^4.*k_2.^2.*exp(k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(-k_1.*(T - tau1)) + 2.*k_1.^2.*k_2.^5.*tau1.*exp(k_1.*(T - tau1)) - 6.*k_1.^4.*k_2.^3.*tau1.*exp(k_1.*(T - tau1)) + 4.*k_1.^5.*k_2.^2.*tau1.*exp(k_1.*(T - tau1)) + 2.*k_1.^2.*k_2.^5.*tau1.*exp(-k_1.*(T - tau1)) - 4.*k_1.^3.*k_2.^4.*tau1.*exp(-k_1.*(T - tau1)) + 2.*k_1.^4.*k_2.^3.*tau1.*exp(-k_1.*(T - tau1)) + 3.*alpha.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) + alpha.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) + 2.*alpha.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*alpha.*k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*S.*k_1.^4.*k_2.*exp(k_2.*tau1) - 2.*S.*k_1.^4.*k_2.*exp(-k_2.*tau1) - 4.*alpha.*k_1.*k_2.^3.*exp(k_2.*tau1) - alpha.*k_1.^3.*k_2.*exp(k_2.*tau1) - 4.*alpha.*k_1.*k_2.^3.*exp(-k_2.*tau1) + alpha.*k_1.^3.*k_2.*exp(-k_2.*tau1) - 4.*S.*T.*k_1.^2.*k_2.^4 + 2.*epsilon.*k_1.^4.*k_2.*exp(k_2.*tau1) - 2.*epsilon.*k_1.^4.*k_2.*exp(-k_2.*tau1) + 2.*S.*k_1.*k_2.^4.*exp(k_1.*(T - tau1)) - 2.*S.*k_1.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*L.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*L.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^5.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*T.*k_1.^4.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*T.*k_1.^5.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 4.*T.*epsilon.*k_1.^2.*k_2.^4 - 5.*alpha.*k_1.*k_2.^3.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_1.^3.*k_2.*exp(k_1.*(T - tau1)) - 3.*alpha.*k_1.*k_2.^3.*exp(-k_1.*(T - tau1)) - 2.*alpha.*k_1.^3.*k_2.*exp(-k_1.*(T - tau1)) + 2.*epsilon.*k_1.*k_2.^4.*exp(k_1.*(T - tau1)) - 2.*epsilon.*k_1.*k_2.^4.*exp(-k_1.*(T - tau1)) - 3.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*S.*k_1.^2.*k_2.^4.*tau1 - 4.*S.*k_1.^4.*k_2.^2.*tau1 + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*epsilon.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^3.*tau1 + 2.*alpha.*k_1.^3.*k_2.^2.*tau1 + 2.*k_1.^4.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^5.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^3.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*k_1.^4.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*k_1.^5.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*epsilon.*k_1.^2.*k_2.^4.*tau1 - 4.*epsilon.*k_1.^4.*k_2.^2.*tau1 - 2.*S.*k_1.^3.*k_2.^2.*exp(k_2.*tau1) + 4.*S.*k_1.^2.*k_2.^3.*exp(-k_2.*tau1) - 2.*S.*k_1.^3.*k_2.^2.*exp(-k_2.*tau1) - 2.*L.*k_1.^2.*k_2.^4.*exp(k_1.*(T - tau1)) - 2.*L.*k_1.^3.*k_2.^3.*exp(k_1.*(T - tau1)) + 4.*L.*k_1.^4.*k_2.^2.*exp(k_1.*(T - tau1)) - 2.*L.*k_1.^2.*k_2.^4.*exp(-k_1.*(T - tau1)) + 2.*L.*k_1.^3.*k_2.^3.*exp(-k_1.*(T - tau1)) + 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*tau1) + 4.*epsilon.*k_1.^2.*k_2.^3.*exp(-k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(-k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^3.*exp(k_1.*(T - tau1)) - 4.*S.*k_1.^3.*k_2.^2.*exp(k_1.*(T - tau1)) + 2.*S.*k_1.^2.*k_2.^3.*exp(-k_1.*(T - tau1)) - T.*alpha.*k_1.*k_2.^4.*exp(k_1.*(T - tau1)) + T.*alpha.*k_1.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^4.*tau1.*exp(k_1.*(T - tau1)) - alpha.*k_1.*k_2.^4.*tau1.*exp(-k_1.*(T - tau1)) + 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 3.*T.*alpha.*k_1.^2.*k_2.^3.*exp(k_1.*(T - tau1)) - 2.*T.*alpha.*k_1.^3.*k_2.^2.*exp(k_1.*(T - tau1)) - T.*alpha.*k_1.^2.*k_2.^3.*exp(-k_1.*(T - tau1)) - 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_1.*(T - tau1)) + alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(-k_1.*(T - tau1)) - alpha.*k_1.^3.*k_2.^2.*tau1.*exp(-k_1.*(T - tau1)) - alpha.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1)))./(4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 8.*k_2.^4.*exp(-k_2.*tau1) - 8.*k_1.^4.*exp(k_1.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.*k_2.^3.*exp(k_2.*tau1) + 4.*k_1.*k_2.^3.*exp(-k_2.*tau1) + 4.*k_1.^3.*k_2.*exp(k_1.*(T - tau1)) - 4.*k_1.^3.*k_2.*exp(-k_1.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*(T - tau1)) + 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 3.*k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*T.*k_1.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1));
    
    psi_2_t_tau1_second = @(t, tau1) psi_2_t_1_0_tau1_second(t, psi_1_cur_tau1(tau1), psi_2_0_cur_tau1(tau1), tau1);
    x_2_t_tau1_second = @(t, tau1) x_2_t_1_0_tau1_second(t, psi_1_cur_tau1(tau1), psi_2_0_cur_tau1(tau1), tau1);
    
    tau1_cur = fsolve(@(tau1) psi_1_cur_tau1(tau1) + psi_2_t_tau1_second(tau1, tau1) .* x_2_t_tau1_second(tau1, tau1), T ./ 2, fsolve_opts);

    if (0 < tau1_cur && tau1_cur < T)
        psi_1_cur = psi_1_cur_tau1(tau1_cur);
        psi_2_0_cur = psi_2_0_cur_tau1(tau1_cur);
        
        psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);

        if (psi_1_cur > 0 && ...
            psi_2_t_second(T) < 0)

            u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
            u_1_t_second = @(t) u_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);
            
            psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);
            x_2_t_second = @(t) x_2_t_tau1_second(t, tau1_cur);
            
            index = find(t_opt - tau1_cur >= 0, 1);
            
            if (prod(psi_1_cur + psi_2_t_second(t_opt(index : end)) .* x_2_t_second(t_opt(index : end)) <= 0))

                functional_cur = integral(u_1_t_first, 0, tau1_cur) + integral(u_1_t_second, tau1_cur, T);

                if (functional_cur < functional_min)

                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    index = find(t_opt - tau1_cur >= 0, 1);

                    psi_2_t_first =  @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
                    x_1_t_second = @(t) x_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);

                    x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_final .* ones(1, number_of_points_for_splitting);
                    psi2_opt = [psi_2_t_first(t_opt(1 : (index - 1))), psi_2_t_second(t_opt(index : end))];

                    x1_opt = [x_1_t_first(t_opt(1 : (index - 1))), x_1_t_second(t_opt(index : end))];
                    x2_opt = [x_2_t_first(t_opt(1 : (index - 1))), x_2_t_second(t_opt(index : end))];

                    u1_opt = [u_1_t_first(t_opt(1 : (index - 1))), u_1_t_second(t_opt(index : end))];
                    u2_opt = [u2_first .* ones(1, index - 1), u2_second .* ones(1, number_of_points_for_splitting - index + 1)];

                    functional_min = functional_cur;

                    switches = [tau1_cur];

                end
            end
        end
    end
    
    % psi_2(T) > 0
    
    psi_1_cur_tau1 = @(tau1) -(2.*k_1.*k_2.*(4.*S.*k_1.^3.*k_2 - 4.*alpha.*k_1.^3.*exp(k_1.*(T - tau1)) - 4.*alpha.*k_2.^3.*exp(-k_2.*tau1) + 2.*alpha.*k_1.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*alpha.*k_1.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*alpha.*k_1.*k_2.^2 - 2.*alpha.*k_1.^2.*k_2 - 4.*epsilon.*k_1.^3.*k_2 - 2.*S.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + alpha.*k_1.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - alpha.*k_1.^2.*k_2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*alpha.*k_1.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*S.*k_1.^3.*k_2.*exp(k_2.*tau1) + 4.*S.*k_1.*k_2.^3.*exp(-k_2.*tau1) - 2.*S.*k_1.^3.*k_2.*exp(-k_2.*tau1) - alpha.*k_1.*k_2.^2.*exp(k_2.*tau1) + alpha.*k_1.^2.*k_2.*exp(k_2.*tau1) + 3.*alpha.*k_1.*k_2.^2.*exp(-k_2.*tau1) + alpha.*k_1.^2.*k_2.*exp(-k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.*exp(k_2.*tau1) - 4.*epsilon.*k_1.*k_2.^3.*exp(-k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.*exp(-k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*L.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*T.*k_1.^3.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^3.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^2.*exp(k_1.*(T - tau1)) + 3.*alpha.*k_1.^2.*k_2.*exp(k_1.*(T - tau1)) + alpha.*k_1.*k_2.^2.*exp(-k_1.*(T - tau1)) - alpha.*k_1.^2.*k_2.*exp(-k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*k_1.^4.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*k_1.^3.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*k_1.^4.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) - 2.*epsilon.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) - T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + T.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^2.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - T.*alpha.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - alpha.*k_1.^3.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1)))./(4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 8.*k_2.^4.*exp(-k_2.*tau1) - 8.*k_1.^4.*exp(k_1.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.*k_2.^3.*exp(k_2.*tau1) + 4.*k_1.*k_2.^3.*exp(-k_2.*tau1) + 4.*k_1.^3.*k_2.*exp(k_1.*(T - tau1)) - 4.*k_1.^3.*k_2.*exp(-k_1.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*(T - tau1)) + 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 3.*k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*T.*k_1.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1));
    psi_2_0_cur_tau1 = @(tau1) (2.*(4.*L.*k_1.^2.*k_2.^4 - 4.*alpha.*k_2.^4 - 4.*alpha.*k_1.^4.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_2.^4.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_2.^4.*exp(-k_1.*(T - tau1)) - 4.*S.*k_1.^2.*k_2.^3 + 4.*S.*k_1.^3.*k_2.^2 - 4.*T.*k_1.^3.*k_2.^4 - 8.*alpha.*k_1.^2.*k_2.^2 + 4.*epsilon.*k_1.^2.*k_2.^3 - 4.*epsilon.*k_1.^3.*k_2.^2 - 4.*epsilon.*k_1.^2.*k_2.^4 - 4.*k_1.^2.*k_2.^5.*tau1 + 4.*k_1.^3.*k_2.^4.*tau1 + 2.*alpha.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*alpha.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 8.*alpha.*k_1.*k_2.^3 + 2.*T.*k_1.^3.*k_2.^4.*exp(k_1.*(T - tau1)) + 2.*T.*k_1.^4.*k_2.^3.*exp(k_1.*(T - tau1)) - 4.*T.*k_1.^5.*k_2.^2.*exp(k_1.*(T - tau1)) + 2.*T.*k_1.^3.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*T.*k_1.^4.*k_2.^3.*exp(-k_1.*(T - tau1)) + 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*(T - tau1)) + 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_1.*(T - tau1)) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*(T - tau1)) + 4.*epsilon.*k_1.^3.*k_2.^2.*exp(k_1.*(T - tau1)) - 2.*epsilon.*k_1.^2.*k_2.^3.*exp(-k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^4.*exp(k_1.*(T - tau1)) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_1.*(T - tau1)) - 4.*epsilon.*k_1.^4.*k_2.^2.*exp(k_1.*(T - tau1)) + 2.*epsilon.*k_1.^2.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(-k_1.*(T - tau1)) + 2.*k_1.^2.*k_2.^5.*tau1.*exp(k_1.*(T - tau1)) - 6.*k_1.^4.*k_2.^3.*tau1.*exp(k_1.*(T - tau1)) + 4.*k_1.^5.*k_2.^2.*tau1.*exp(k_1.*(T - tau1)) + 2.*k_1.^2.*k_2.^5.*tau1.*exp(-k_1.*(T - tau1)) - 4.*k_1.^3.*k_2.^4.*tau1.*exp(-k_1.*(T - tau1)) + 2.*k_1.^4.*k_2.^3.*tau1.*exp(-k_1.*(T - tau1)) + 3.*alpha.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) + alpha.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) + 2.*alpha.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*alpha.*k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*alpha.*k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*S.*k_1.^4.*k_2.*exp(k_2.*tau1) - 2.*S.*k_1.^4.*k_2.*exp(-k_2.*tau1) - 4.*alpha.*k_1.*k_2.^3.*exp(k_2.*tau1) - alpha.*k_1.^3.*k_2.*exp(k_2.*tau1) - 4.*alpha.*k_1.*k_2.^3.*exp(-k_2.*tau1) + alpha.*k_1.^3.*k_2.*exp(-k_2.*tau1) - 4.*S.*T.*k_1.^2.*k_2.^4 - 2.*epsilon.*k_1.^4.*k_2.*exp(k_2.*tau1) + 2.*epsilon.*k_1.^4.*k_2.*exp(-k_2.*tau1) + 2.*S.*k_1.*k_2.^4.*exp(k_1.*(T - tau1)) - 2.*S.*k_1.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*L.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*L.*k_1.^3.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*L.*k_1.^3.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*L.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*S.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*S.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^5.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*T.*k_1.^4.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.^4.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*T.*k_1.^5.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*T.*epsilon.*k_1.^2.*k_2.^4 - 5.*alpha.*k_1.*k_2.^3.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_1.^3.*k_2.*exp(k_1.*(T - tau1)) - 3.*alpha.*k_1.*k_2.^3.*exp(-k_1.*(T - tau1)) - 2.*alpha.*k_1.^3.*k_2.*exp(-k_1.*(T - tau1)) - 2.*epsilon.*k_1.*k_2.^4.*exp(k_1.*(T - tau1)) + 2.*epsilon.*k_1.*k_2.^4.*exp(-k_1.*(T - tau1)) - 3.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*S.*k_1.^2.*k_2.^4.*tau1 - 4.*S.*k_1.^4.*k_2.^2.*tau1 - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*epsilon.*k_1.^3.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*epsilon.*k_1.^4.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*alpha.*k_1.^2.*k_2.^3.*tau1 + 2.*alpha.*k_1.^3.*k_2.^2.*tau1 + 2.*k_1.^4.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^5.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^3.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*k_1.^4.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^4.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*k_1.^5.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 4.*epsilon.*k_1.^2.*k_2.^4.*tau1 + 4.*epsilon.*k_1.^4.*k_2.^2.*tau1 - 2.*S.*k_1.^3.*k_2.^2.*exp(k_2.*tau1) + 4.*S.*k_1.^2.*k_2.^3.*exp(-k_2.*tau1) - 2.*S.*k_1.^3.*k_2.^2.*exp(-k_2.*tau1) - 2.*L.*k_1.^2.*k_2.^4.*exp(k_1.*(T - tau1)) - 2.*L.*k_1.^3.*k_2.^3.*exp(k_1.*(T - tau1)) + 4.*L.*k_1.^4.*k_2.^2.*exp(k_1.*(T - tau1)) - 2.*L.*k_1.^2.*k_2.^4.*exp(-k_1.*(T - tau1)) + 2.*L.*k_1.^3.*k_2.^3.*exp(-k_1.*(T - tau1)) + 5.*alpha.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 3.*alpha.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(k_2.*tau1) - 4.*epsilon.*k_1.^2.*k_2.^3.*exp(-k_2.*tau1) + 2.*epsilon.*k_1.^3.*k_2.^2.*exp(-k_2.*tau1) + 2.*S.*k_1.^2.*k_2.^3.*exp(k_1.*(T - tau1)) - 4.*S.*k_1.^3.*k_2.^2.*exp(k_1.*(T - tau1)) + 2.*S.*k_1.^2.*k_2.^3.*exp(-k_1.*(T - tau1)) - T.*alpha.*k_1.*k_2.^4.*exp(k_1.*(T - tau1)) + T.*alpha.*k_1.*k_2.^4.*exp(-k_1.*(T - tau1)) - 2.*T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - T.*alpha.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + T.*alpha.*k_1.^2.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + T.*alpha.*k_1.^3.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.*k_2.^4.*tau1.*exp(k_1.*(T - tau1)) - alpha.*k_1.*k_2.^4.*tau1.*exp(-k_1.*(T - tau1)) + 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*alpha.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 3.*T.*alpha.*k_1.^2.*k_2.^3.*exp(k_1.*(T - tau1)) - 2.*T.*alpha.*k_1.^3.*k_2.^2.*exp(k_1.*(T - tau1)) - T.*alpha.*k_1.^2.*k_2.^3.*exp(-k_1.*(T - tau1)) - 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(k_1.*(T - tau1)) + alpha.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*(T - tau1)) + 2.*alpha.*k_1.^2.*k_2.^3.*tau1.*exp(-k_1.*(T - tau1)) - alpha.*k_1.^3.*k_2.^2.*tau1.*exp(-k_1.*(T - tau1)) - alpha.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + alpha.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1)))./(4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 4.*k_1.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 4.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 8.*k_2.^4.*exp(-k_2.*tau1) - 8.*k_1.^4.*exp(k_1.*(T - tau1)) - 8.*k_1.^2.*k_2.^2 - 4.*k_1.*k_2.^3.*exp(k_2.*tau1) + 4.*k_1.*k_2.^3.*exp(-k_2.*tau1) + 4.*k_1.^3.*k_2.*exp(k_1.*(T - tau1)) - 4.*k_1.^3.*k_2.*exp(-k_1.*(T - tau1)) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^2.*k_2.^2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(-k_2.*tau1) + 4.*k_1.^2.*k_2.^2.*exp(k_1.*(T - tau1)) + 4.*k_1.^2.*k_2.^2.*exp(-k_1.*(T - tau1)) + 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 + k_2.*tau1) - 3.*k_1.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - k_1.*k_2.^3.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 3.*k_1.^3.*k_2.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 3.*k_1.^3.*k_2.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*T.*k_1.*k_2.^4.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*T.*k_1.*k_2.^4.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) - 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*k_1.*k_2.^4.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.*k_2.^4.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1) + 2.*k_1.^4.*k_2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) + 2.*T.*k_1.^2.*k_2.^3.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 + k_2.*tau1) - 2.*k_1.^2.*k_2.^3.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) - 2.*k_1.^3.*k_2.^2.*tau1.*exp(T.*k_1 - k_1.*tau1 - k_2.*tau1) + 2.*k_1.^3.*k_2.^2.*tau1.*exp(k_1.*tau1 - T.*k_1 - k_2.*tau1));
    
    tau1_cur = fsolve(@(tau1) psi_1_cur_tau1(tau1) + psi_2_t_tau1_second(tau1, tau1) .* x_2_t_tau1_second(tau1, tau1), T ./ 2, fsolve_opts);

    if (0 < tau1_cur && tau1_cur < T)
        psi_1_cur = psi_1_cur_tau1(tau1_cur);
        psi_2_0_cur = psi_2_0_cur_tau1(tau1_cur);
        
        psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);

        if (psi_1_cur > 0 && ...
            psi_2_t_second(T) > 0)

            u_1_t_first = @(t) u_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
            u_1_t_second = @(t) u_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);
            
            psi_2_t_second = @(t) psi_2_t_tau1_second(t, tau1_cur);
            x_2_t_second = @(t) x_2_t_tau1_second(t, tau1_cur);
            
            index = find(t_opt - tau1_cur >= 0, 1);
            
            if (prod(psi_1_cur + psi_2_t_second(t_opt(index : end)) .* x_2_t_second(t_opt(index : end)) <= 0))

                functional_cur = integral(u_1_t_first, 0, tau1_cur) + integral(u_1_t_second, tau1_cur, T);

                if (functional_cur < functional_min)

                    t_opt = linspace(0, T, number_of_points_for_splitting);

                    index = find(t_opt - tau1_cur >= 0, 1);

                    psi_2_t_first =  @(t) psi_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    x_1_t_first = @(t) x_1_t_1_0_first(t, psi_1_cur, psi_2_0_cur);
                    x_1_t_second = @(t) x_1_t_1_0_tau1_second(t, psi_1_cur, psi_2_0_cur, tau1_cur);

                    x_2_t_first = @(t) x_2_t_1_0_first(t, psi_1_cur, psi_2_0_cur);

                    psi1_opt = psi_1_final .* ones(1, number_of_points_for_splitting);
                    psi2_opt = [psi_2_t_first(t_opt(1 : (index - 1))), psi_2_t_second(t_opt(index : end))];

                    x1_opt = [x_1_t_first(t_opt(1 : (index - 1))), x_1_t_second(t_opt(index : end))];
                    x2_opt = [x_2_t_first(t_opt(1 : (index - 1))), x_2_t_second(t_opt(index : end))];

                    u1_opt = [u_1_t_first(t_opt(1 : (index - 1))), u_1_t_second(t_opt(index : end))];
                    u2_opt = [u2_first .* ones(1, index - 1), u2_second .* ones(1, number_of_points_for_splitting - index + 1)];

                    functional_min = functional_cur;

                    switches = [tau1_cur];

                end
            end
        end
    end
    
    % two or more switches
    
    psi_1_cur_1_tau1_tau2_0 = @(tau1, tau2, psi_2_0) (k_2.*(2.*psi_2_0 + exp(k_2.*tau1).*(4.*alpha.^2 - 16.*k_2.^2.*psi_2_0 - 4.*alpha.*psi_2_0 - 16.*alpha.^2.*exp(k_2.*tau1) + 24.*alpha.^2.*exp(2.*k_2.*tau1) - 16.*alpha.^2.*exp(3.*k_2.*tau1) + 4.*alpha.^2.*exp(4.*k_2.*tau1) + 16.*k_2.^4.*exp(2.*k_2.*tau1) + psi_2_0.^2 - 4.*psi_2_0.^2.*exp(k_2.*tau1) + 6.*psi_2_0.^2.*exp(2.*k_2.*tau1) - 4.*psi_2_0.^2.*exp(3.*k_2.*tau1) + psi_2_0.^2.*exp(4.*k_2.*tau1) + 16.*alpha.*k_2.^2.*exp(k_2.*tau1) - 32.*alpha.*k_2.^2.*exp(2.*k_2.*tau1) + 16.*alpha.*k_2.^2.*exp(3.*k_2.*tau1) + 24.*k_2.^2.*psi_2_0.*exp(k_2.*tau1) - 8.*k_2.^2.*psi_2_0.*exp(3.*k_2.*tau1) + 16.*alpha.*psi_2_0.*exp(k_2.*tau1) - 24.*alpha.*psi_2_0.*exp(2.*k_2.*tau1) + 16.*alpha.*psi_2_0.*exp(3.*k_2.*tau1) - 4.*alpha.*psi_2_0.*exp(4.*k_2.*tau1)).^(1./2) - 2.*alpha.*exp(k_2.*tau1) + 4.*alpha.*exp(2.*k_2.*tau1) - 2.*alpha.*exp(3.*k_2.*tau1) - 3.*psi_2_0.*exp(k_2.*tau1) + psi_2_0.*exp(3.*k_2.*tau1) - 4.*k_2.^2.*exp(2.*k_2.*tau1)))./(2.*(exp(k_2.*tau1) - 1).^3);
    psi_1_cur_2_tau1_tau2_0 = @(tau1, tau2, psi_2_0) -(k_2.*(exp(k_2.*tau1).*(4.*alpha.^2 - 16.*k_2.^2.*psi_2_0 - 4.*alpha.*psi_2_0 - 16.*alpha.^2.*exp(k_2.*tau1) + 24.*alpha.^2.*exp(2.*k_2.*tau1) - 16.*alpha.^2.*exp(3.*k_2.*tau1) + 4.*alpha.^2.*exp(4.*k_2.*tau1) + 16.*k_2.^4.*exp(2.*k_2.*tau1) + psi_2_0.^2 - 4.*psi_2_0.^2.*exp(k_2.*tau1) + 6.*psi_2_0.^2.*exp(2.*k_2.*tau1) - 4.*psi_2_0.^2.*exp(3.*k_2.*tau1) + psi_2_0.^2.*exp(4.*k_2.*tau1) + 16.*alpha.*k_2.^2.*exp(k_2.*tau1) - 32.*alpha.*k_2.^2.*exp(2.*k_2.*tau1) + 16.*alpha.*k_2.^2.*exp(3.*k_2.*tau1) + 24.*k_2.^2.*psi_2_0.*exp(k_2.*tau1) - 8.*k_2.^2.*psi_2_0.*exp(3.*k_2.*tau1) + 16.*alpha.*psi_2_0.*exp(k_2.*tau1) - 24.*alpha.*psi_2_0.*exp(2.*k_2.*tau1) + 16.*alpha.*psi_2_0.*exp(3.*k_2.*tau1) - 4.*alpha.*psi_2_0.*exp(4.*k_2.*tau1)).^(1./2) - 2.*psi_2_0 + 2.*alpha.*exp(k_2.*tau1) - 4.*alpha.*exp(2.*k_2.*tau1) + 2.*alpha.*exp(3.*k_2.*tau1) + 3.*psi_2_0.*exp(k_2.*tau1) - psi_2_0.*exp(3.*k_2.*tau1) + 4.*k_2.^2.*exp(2.*k_2.*tau1)))./(2.*(exp(k_2.*tau1) - 1).^3);

    solving_equat = @(tau1, tau2, psi_1, psi_2_0) psi_1 - (exp(k_1.*(tau1 - tau2)).*(psi_2_0.*exp(-k_2.*tau1) + (psi_1.*(exp(-k_2.*tau1) - 1))./k_2) + (psi_1.*(exp(k_1.*(tau1 - tau2)) - 1))./k_1).*((alpha.*(exp(-k_1.*(tau1 - tau2)) - 1))./(2.*k_1) + (psi_1.*(exp(-k_1.*(tau1 - tau2)) - 1))./(2.*k_1.^2) + (psi_1.*exp(- k_1.*tau1 - k_1.*tau2).*(exp(2.*k_1.*tau1) - exp(2.*k_1.*tau2)))./(4.*k_1.^2) + (psi_2_0.*exp(- k_1.*tau1 - k_1.*tau2).*exp(-k_2.*tau1).*(exp(2.*k_1.*tau1) - exp(2.*k_1.*tau2)))./(4.*k_1) - (psi_1.*exp(- k_1.*tau1 - k_1.*tau2).*(exp(2.*k_1.*tau1) - exp(2.*k_1.*tau2)))./(4.*k_1.*k_2) - (exp(-k_2.*tau1).*exp(-k_1.*(tau1 - tau2)).*(exp(k_2.*tau1) - 1).*(psi_1 + k_2.*psi_2_0 - psi_1.*exp(k_2.*tau1) - 2.*alpha.*k_2.*exp(k_2.*tau1) + k_2.*psi_2_0.*exp(k_2.*tau1)))./(4.*k_2.^2) + (psi_1.*exp(- k_1.*tau1 - k_1.*tau2).*exp(-k_2.*tau1).*(exp(2.*k_1.*tau1) - exp(2.*k_1.*tau2)))./(4.*k_1.*k_2));
    
    psi_2_0_cur_11_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_1_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), alpha, fsolve_opts);
    
    psi_2_0_cur_12_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_1_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), -alpha, fsolve_opts);
    
    psi_2_0_cur_21_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_2_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), alpha, fsolve_opts);
    
    psi_2_0_cur_22_tau1_tau2 = @(tau1, tau2) fsolve(@(psi_2_0) solving_equat(...
        tau1, tau2, psi_1_cur_2_tau1_tau2_0(tau1, tau2, psi_2_0), psi_2_0), -alpha, fsolve_opts);

    tau1_splitting = linspace(0 + trimmer, T - trimmer, number_of_points_for_enum);
        
    counter = 0;
    for tau1_cur = tau1_splitting
        counter = counter + 1;
        tau2_splitting = tau1_splitting((counter + 1) : end);
        for tau2_cur = tau2_splitting

            psi_2_0_cur_11 = psi_2_0_cur_11_tau1_tau2(tau1_cur, tau2_cur);
            psi_2_0_cur_12 = psi_2_0_cur_12_tau1_tau2(tau1_cur, tau2_cur);
            psi_2_0_cur_21 = psi_2_0_cur_21_tau1_tau2(tau1_cur, tau2_cur);
            psi_2_0_cur_22 = psi_2_0_cur_22_tau1_tau2(tau1_cur, tau2_cur);
            
            psi_1_cur_11 = psi_1_cur_1_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_11);
            psi_1_cur_12 = psi_1_cur_1_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_12);
            psi_1_cur_21 = psi_1_cur_2_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_21);
            psi_1_cur_22 = psi_1_cur_2_tau1_tau2_0(tau1_cur, tau2_cur, psi_2_0_cur_22);

            if (imag(psi_1_cur_11) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_2, k_1, alpha, psi_1_cur_11), [0, T], [0, 0, psi_2_0_cur_11]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_11 > 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_11 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_11 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_11 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;

                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);
                    end
                end
            end    

            if (imag(psi_1_cur_12) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_2, k_1, alpha, psi_1_cur_12), [0, T], [0, 0, psi_2_0_cur_12]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_12 > 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_12 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_12 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_12 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;

                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);
                    end
                end
            end    
            
            if (imag(psi_1_cur_21) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_2, k_1, alpha, psi_1_cur_11), [0, T], [0, 0, psi_2_0_cur_21]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_21 > 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_21 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_21 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_21 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;
                        
                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);

                    end
                end
            end    
            
            if (imag(psi_1_cur_22) == 0)
                [t_cur, x_cur] = ode45(@(t, x) odefun_first(t, x, k_2, k_1, alpha, psi_1_cur_22), [0, T], [0, 0, psi_2_0_cur_22]);
                x_1_vect = x_cur(:, 1);
                x_2_vect = x_cur(:, 2);
                psi_2_vect = x_cur(:, 3);
                if (psi_1_cur_22 > 0 && ...
                    abs(L - x_1_vect(end)) <= epsilon + delta && ...
                    abs(S - x_2_vect(end)) <= epsilon + delta)

                    u_1_vect = (psi_2_vect - alpha) ./ 2;

                    functional_cur = trapz(t_cur, u_1_vect .^ 2 + alpha .* u_1_vect);

                    if (functional_cur < functional_min)

                        u_2_vect = (psi_1_cur_22 + psi_2_vect .* x_2_vect > 0) .* k_2 + ...
                                   (psi_1_cur_22 + psi_2_vect .* x_2_vect < 0) .* k_1;

                        t_opt = t_cur.';

                        psi1_opt = psi_1_cur_22 .* ones(1, length(t_cur));
                        psi2_opt = psi_2_vect.';

                        x1_opt = x_1_vect.';
                        x2_opt = x_2_vect.';

                        u1_opt = u_1_vect.';
                        u2_opt = u_2_vect.';

                        functional_min = functional_cur;

                        switches = t_opt(u_2_vect(1 : (end - 1)) .* u_2_vect(2 : end) < 0);
                    end
                end
            end    
        end
    end
end

try
    u_opt = [u1_opt; u2_opt; t_opt];
    x_opt = [x1_opt; x2_opt; t_opt];
    psi_opt = [psi1_opt; psi2_opt; t_opt];
catch
    u_opt = [];
    x_opt = [];
    psi_opt = [];
    disp('No Solution!');
end
end