function [value, isterminal, direction] = event_X1(t, x, d, k, eps)
value = abs(x(2) - d) + x(1) .^ 2 - k - eps;
isterminal = 1;
direction = -1;
end