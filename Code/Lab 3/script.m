% Laboratory Work #3

%% Function #1

alpha = 0.1;
t = 6;

dist = 0.5;

N = 100;

K = 3;

plot(0, 0, 'o');
grid on;
hold on;

plot_traj = false;

[X, Y, switches_line] = reachset(alpha, t, N, plot_traj, dist, K);

plot(X, Y, 'b', 'LineWidth', 3);
plot(switches_line(:, 1), switches_line(:, 2), '.r');

%% Function #2

alpha = 0.1;

t1 = 1;
t2 = 9;

N = 10;

K = 10;

K_reach = 3;

dist = 0.5;

filename = 'First_Exp';

time_between_frames = 0.01;

reachsetdyn(alpha, t1, t2, N, filename, time_between_frames, K, K_reach, dist, plot_traj);

%reachsetdyn(alpha, t1, t2, N, [], time_between_frames, K, dist);
