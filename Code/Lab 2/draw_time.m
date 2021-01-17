function [plot_opt, plot_switches] = draw_time(fig, func, i, switches)
figure(fig);
hold on;
func_i = func(i, :);
time = func(3, :);
plot_opt = plot(time, func_i, 'r');
plot_switches = [];
for j = 1 : length(switches)
    index = find(time - switches(j) > 0, 1);
    plot_switches = plot(time(index), func_i(index), '*k');
end
end