function [value, isterminal, direction] = event_func_x_2(t, x)
value = x(2) + 0 .* t;
isterminal = 1;
direction = 0;
end