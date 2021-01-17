psi_1 = 1;
alpha = 1;

k_1 = -0.5;
k_2 = 1;

left_point = -10;
right_point = 10;

step = 0.05;
big_step = 1;

trimmer = 0.5;

psi_2_u_1 = @(x_2) 0 .* x_2 + alpha;
psi_2_u_2 = @(x_2) -psi_1 ./ x_2;

new_left_point = (alpha - sqrt(alpha .^ 2 + 8 * psi_1 * k_1)) / 2;
new_right_point = (alpha + sqrt(alpha .^ 2 + 8 * psi_1 * k_1)) / 2;

old_left_point = (alpha - sqrt(alpha .^ 2 + 8 * psi_1 * k_2)) / 2;
old_right_point = (alpha + sqrt(alpha .^ 2 + 8 * psi_1 * k_2)) / 2;

x_2_start = @(psi_2) 0 .* psi_2;
x_2_cond_1 = @(psi_2) (alpha - psi_2) ./ (2 .* k_1);
x_2_cond_2 = @(psi_2) (alpha - psi_2) ./ (2 .* k_2);

split = left_point : step : right_point;
new_split = new_left_point : step : new_right_point;
if (~isnan(old_left_point))
    old_split = old_left_point : step : old_right_point;
end
%plot_1 = plot(split, psi_2_u_1(split), 'b', 'LineWidth', 3);
%hold on;
plot_2 = plot(split, psi_2_u_2(split), 'm', 'LineWidth', 3);
hold on;
xlabel('x_2');
ylabel('\psi_2');
grid;
plot_3 = plot(x_2_start(split), split, 'r', 'LineWidth', 4);
plot_4 = plot(x_2_cond_1(new_split), new_split, 'c', 'LineWidth', 2);
plot_5 = plot(x_2_cond_2(old_split), old_split, 'b', 'LineWidth', 2);
plot_6 = plot(0, -psi_1 / k_1, '*k', 'LineWidth', 3);
plot_7 = plot(0, -psi_1 / k_2, 'ok', 'LineWidth', 2);
axis([left_point - trimmer, right_point + trimmer, ...
     left_point - trimmer, right_point + trimmer]);

big_split = left_point : big_step : right_point;

for i = big_split
    for j = big_split
        if (psi_1 + i * j > 0)
            quiver(i, j, (j - alpha) / 2 + k_2 * i, -psi_1 - j * k_2, 0.2, 'k');
        end
        if (psi_1 + i * j < 0)
            quiver(i, j, (j - alpha) / 2 + k_1 * i, -psi_1 - j * k_1, 0.2, 'k');
        end
    end
end
hold off;

legend([plot_2(1), plot_3(1), plot_4(1), plot_5(1), plot_6(1), plot_7(1)], ...
       '\psi_1+\psi_2x_2=0', ...
       'x_2(0)=0', ...
       '2k_1x_2=\alpha-\psi_2', ...
       '2k_2x_2=\alpha-\psi_2', ...
       'k_1\psi_2=-\psi_1', ...
       'k_2\psi_2=-\psi_1');