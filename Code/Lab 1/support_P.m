function [func, argmax] = support_P (l_vect, a, b, c)
size_l = size(l_vect);
len_l = size_l(2);

% func = zeros(1, len_l);
argmax = zeros(2, len_l);

case_l2 = (l_vect(2, :) >= 0);

case_l1_1 = (abs(l_vect(1, :)) <= 2 .* c .* l_vect(2, :) .* sqrt(b ./ (a + c)));
case_l1_2 = (abs(l_vect(1, :)) > 2 .* a .* abs(l_vect(2, :)) .* sqrt(b ./ (a + c)));

cond_1 = logical(case_l2 .* case_l1_1);
cond_2 = logical((case_l2 .* (~case_l1_1)) + ((~case_l2) .* case_l1_2));
cond_3 = logical((~case_l2) .* (~case_l1_2));

l1_1 = l_vect(1, cond_1);
l2_1 = l_vect(2, cond_1);
l1_2 = l_vect(1, cond_2);
% l2_2 = l_vect(2, cond_2);
l1_3 = l_vect(1, cond_3);
l2_3 = l_vect(2, cond_3);

% func(cond_1) = b .* l2_1 + (l1_1 .^ 2) ./ (4 .* c .* l2_1);
% func(cond_2) = abs(l1_2) .* sqrt(b ./ (a + c)) + l2_2 .* (a .* b) .* (a + c);
% func(cond_3) = (l1_3 .^ 2) ./ (4 .* a .* (-l2_3));

argmax(:, cond_1) = [l1_1 ./ (2 .* c .* l2_1); b - (l1_1 .^ 2) ./ (4 .* c .* (l2_1 .^ 2))];
argmax(:, cond_2) = [sign(l1_2) .* sqrt(b ./ (a + c)); a .* b ./ (a + c)];
argmax(:, cond_3) = [-l1_3 ./ (2 .* a .* l2_3); (l1_3 .^ 2) ./ (4 .* a .* (l2_3 .^ 2))];

func = l_vect(1, :) .* argmax(1, :) + l_vect(2, :) .* argmax(2, :);

end