function [value, isterminal, direction] = event_func_psi_2(t, x)
value = x(4) + 0 .* t;
isterminal = 1;
direction = 0;
end