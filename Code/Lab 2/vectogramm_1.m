psi_1 = 1;
alpha = 2;

k_1 = -1;
k_2 = 1;

left_point = -10;
right_point = 10;

step = 0.05;
big_step = 1;

trimmer = 0.5;

psi_2_u_1 = @(x_2) 0 .* x_2 + alpha;
psi_2_u_2 = @(x_2) -psi_1 ./ x_2;
psi_2_cond = @(x_2) 0 .* x_2 -psi_1 ./ k_2;

x_2_start = @(psi_2) 0 .* psi_2;

split = left_point : step : right_point;
new_split = left_point : step : k_2;
plot_1 = plot(split, psi_2_u_1(split), 'b', 'LineWidth', 3);
hold on;
plot_2 = plot(split, psi_2_u_2(split), 'm', 'LineWidth', 3);
xlabel('x_2');
ylabel('\psi_2');
grid;
plot_3 = plot(x_2_start(split), split, 'r', 'LineWidth', 4);
%plot_4 = plot(new_split, psi_2_cond(new_split), 'c', 'LineWidth', 2);
%plot_5 = plot(0, -psi_1/k_1, '*k', 'LineWidth', 3);

axis([left_point - trimmer, right_point + trimmer, ...
     left_point - trimmer, right_point + trimmer]);

big_split = left_point : big_step : right_point;

for i = big_split
    for j = big_split
        if (j >= alpha) && (psi_1 + i * j > 0)
            quiver(i, j, (j - alpha) / 2 + k_2 * i, -psi_1 - j * k_2, 0.2, 'k');
        end
        if (j < alpha) && (psi_1 + i * j > 0)
            quiver(i, j, k_2 * i, -psi_1 - j * k_2, 0.2, 'k');
        end
        if (j >= alpha) && (psi_1 + i * j < 0)
            quiver(i, j, (j - alpha) / 2 + k_1 * i, -psi_1 - j * k_1, 0.2, 'k');
        end
        if (j < alpha) && (psi_1 + i * j < 0)
            quiver(i, j, k_1 * i, -psi_1 - j * k_1, 0.2, 'k');
        end
    end
end
hold off;

legend([plot_1(1), plot_2(1), plot_3(1)], ... %plot_4(1), plot_5(1)], ...
       '\psi_2=\alpha', ...
       '\psi_1+\psi_2x_2=0', ...
       'x_2(0)=0'); %...
       % 'k_2\psi_2=-\psi_1', ...
       % 'k_1\psi_2=-\psi_1');
