function [plot_cur, plot_opt] = draw_phase(fig, func_opt, func_1, func_2, star)
figure(fig);
hold on;
plot_cur = [];
if (star)
   % plot_cur = cellfun(@ (x, y) plot(x, y, 'b*'), func_1(2, :), func_2(2, :));
     for i = 1 : length(func_1)
         plot_cur = plot(func_1{i}(2, :), func_2{i}(2, :), 'b*');
     end
else
   % plot_cur = cellfun(@ (x, y) plot(x, y, 'b'), func_1(2, :), func_2(2, :));
     for i = 1 : length(func_1)
         plot_cur = plot(func_1{i}(2, :), func_2{i}(2, :), 'b');
     end
end
plot_opt = [];
if (star)
    if (~isempty(func_opt))
        plot_opt = plot(func_opt(:, 1), func_opt(:, 2), 'r*');
    end
else
    if (~isempty(func_opt))
        plot_opt = plot(func_opt(:, 1), func_opt(:, 2), 'r');
    end
end
end