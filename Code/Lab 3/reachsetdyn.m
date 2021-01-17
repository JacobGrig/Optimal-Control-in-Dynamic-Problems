function reachsetdyn(alpha, t1, t2, N, filename, time_between_frames, K, K_reach, dist, plot_traj)

mov(N) = struct('cdata', [], 'colormap', []);

plot(0, 0, 'ob');
hold on;
grid on;

trimmer = 0.1;

[X_max, Y_max] = reachset(alpha, t2, K, plot_traj, dist, K_reach);

axis([min(X_max) - trimmer, max(X_max) + trimmer, min(Y_max) - trimmer, max(Y_max) + trimmer]);

for i = 0 : N
    
    t_cur = t1 + (t2 - t1) * i / N;
    [X, Y] = reachset(alpha, t_cur, K, plot_traj, dist * (i + 1) / N, K_reach);
    
    plot(X, Y, 'b', 'LineWidth', 1);
    
    mov(i + 1) = getframe;
    
    pause(time_between_frames);
    
end

hold off;

if (isempty(filename))
    times = 1;                              % number of repeats
    frames_per_sec = 60;                    % fps
    movie(mov, times, frames_per_sec);
else
    video_object = VideoWriter(strcat(filename, '.avi'));
    open(video_object);
    writeVideo(video_object, mov);
    close(video_object);
end

end