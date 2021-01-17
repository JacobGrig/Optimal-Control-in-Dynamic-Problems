function [plot_opt, plot_switches] = draw_phase(fig, func, switches)
figure(fig);
grid;
hold on;
func1 = func(1, :);
func2 = func(2, :);
time = func(3, :);
plot_opt = plot(func1, func2, 'b');
plot_switches = [];
for i = 1 : length(switches)
    index = find(time - switches(i) > 0, 1);
    plot_switches = plot(func1(index), func2(index), '*k');
end
end