function [func, argmax] = support_X0(l_vect, radius)
norm_l = sqrt(l_vect(1, :) .^ 2 + l_vect(2, :) .^ 2);

argmax = radius .* l_vect ./ norm_l;

func = radius .* norm_l ;
% func = l_vect(1, :) .* argmax(1, :) + l_vect(2, :) .* argmax(2, :);

end