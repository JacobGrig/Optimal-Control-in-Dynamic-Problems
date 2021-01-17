function [func, argmax] = support_X1(l_vect, d, k)
size_l = size(l_vect);
len_l = size_l(2);

% func = zeros(1, len_l);
argmax = zeros(2, len_l);

cond_1 = logical(abs(l_vect(2, :)) < abs(l_vect(1, :)) ./ (2 .* sqrt(k)));
cond_2 = logical(~cond_1);

l1_1 = l_vect(1, cond_1);
l2_1 = l_vect(2, cond_1);
l1_2 = l_vect(1, cond_2);
l2_2 = l_vect(2, cond_2);

% func(cond_1) = d .* l1_1 + sqrt(k) .* abs(l2_1);
% func(cond_2) = abs(l1_2) .* k + (l2_2 .^ 2) ./ (4 .* abs(l1_2)) + d .* l1_2;

argmax(:, cond_1) = [sqrt(k) .* sign(l1_1); d];
argmax(:, cond_2) = [l1_2 ./ (2 .* abs(l2_2)); ...
    d + sign(l2_2) .* (k - (l1_2 .^ 2) ./ (4 .* (l2_2 .^ 2)))];

func = l_vect(1, :) .* argmax(1, :) + l_vect(2, :) .* argmax(2, :);

end