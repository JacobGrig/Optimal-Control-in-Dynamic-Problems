function [plot_cur, plot_opt] = draw_time(fig, func_opt, time_opt, func)
figure(fig);
hold on;
plot_cur = [];
for i = 1 : length(func)
    plot_cur = plot(func{i}(1, :), func{i}(2, :), 'b');
end
plot_opt = [];
if (~isempty(func_opt))
    plot_opt = plot(time_opt, func_opt, 'r');
end
end